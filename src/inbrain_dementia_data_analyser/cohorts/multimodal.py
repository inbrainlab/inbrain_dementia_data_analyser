"""Multimodal cohort creator: subjects with a set of 2, 3, 4 ... modalities close in time.

How a set is chosen, per subject:

* every acquisition of every listed modality is placed on the subject's time
  axis (calendar dates, or days since entry for OASIS / A4);
* a set is valid when all its acquisitions fit inside ``window_months``
  (first-to-last span), it contains every ``required`` modality and it has at
  least ``min_modalities`` modalities in total;
* among the valid sets, ``prefer`` decides: ``most_modalities`` -> the set
  with more modalities, ties -> smallest span; ``tightest`` -> the smallest
  span between the required modalities, ties -> more optional modalities,
  then smallest span;
* one set (the best) per subject.

The search tries each acquisition as the start of the window and takes each
modality's first acquisition inside it; it was checked against exhaustive
search and gives the same subjects and spans.

Modalities are exact stems (``MRI_T2w`` does not pick up ``MRI_T2w_FLAIR``)
or groups (``PET_ABETA``, ``PET_TAU``, ``PLASMA_PTAU217`` ...). Scored
measures (CDR, plasma, CSF ...) can be modalities; their value is returned in
``match_<label>_value``.
"""
from __future__ import annotations

import numpy as np
import os
import pandas as pd
import re
from bisect import bisect_left

from ..config import ProjectConfig, load_config
from ..integration.common import configure as configure_common
from ..provenance import write_provenance
from .outputs import CohortResult

from ..integration.common import VRS_DAYS_PER_MONTH, vrs_events, vrs_measure_stems, vrs_read_expanded

import logging as _logging
from ..logging_utils import make_print as _make_print
from ..logging_utils import quiet_pandas_warnings

log = _logging.getLogger(__name__)
_print = _make_print(log)

# ============================================================
#  DEFAULTS (override them through the function arguments)
# ============================================================
INPUT_PATH: str | None = None   # set by create_multimodal_cohort / the input_path argument

SELECT_DATABANKS = "ALL"     # 'ALL' or ['ADNI', 'OASIS', ...]

# -- Time window ---------------------------------------------
WINDOW_MONTHS = 12          # max span between the first and the last acquisition
                            # of the set (e.g. 6, 12, 18, 24)

# -- Modalities ----------------------------------------------
# Each entry:
#   "label"   : short name used in the output columns  (e.g. "MRI_T1w")
#   "prefixes": measures that satisfy it — exact stems ("MRI_T1w", "CDR",
#               "Plasma_pTau217_simoa") and/or groups ("PET_ABETA", "PET_TAU",
#               "PLASMA_PTAU217", "PLASMA_NFL", "CSF_AB42" …). With several,
#               the one actually used is written in match_<label>_tracer.
#   "required": True  -> the subject is dropped if it is not in the set
#               False -> counted / returned when it fits in the window
#
# To ADD a modality   -> append a dict.     To REMOVE one -> delete its dict.
DEFAULT_MODALITIES = MODALITIES = [
    {
        "label":    "MRI_T1w",
        "prefixes": ["MRI_T1w"],
        "required": True,
    },
    {
        "label":    "MRI_fMRI",
        "prefixes": ["MRI_fMRI"],
        "required": False,
    },
    {
        "label":    "PET_FDG",
        "prefixes": ["PET_FDG"],
        "required": True,
    },
    {
        "label":    "PET_ABETA",
        "prefixes": [
            "ABETA_AV45",
            "ABETA_FBB",
            "ABETA_Flutemetamol",
            "ABETA_NAV4694",
            "ABETA_PIB",
        ],
        "required": False,
    },
    {
        "label":    "PET_TAU",
        "prefixes": [
            "TAU_Flortaucipir",
            "TAU_MK6240",
            "TAU_PI2620",
        ],
        "required": False,
    },
]

# -- How many modalities --------------------------------------
# None -> every "required" modality (bi-modal with 2 required, tri-modal with 3 …)
# k    -> at least k modalities in the set, of which all the required ones.
#         With no required modality at all: ANY k of the list (2, 3, 4 …).
MIN_MODALITIES = None

# -- Which set to keep when several are possible ---------------
PREFER = "most_modalities"      # "most_modalities" | "tightest"

# Subject-level columns copied to the output (when present in the master)
KEEP_COLUMNS = [
    "Databank", "ID",
    "1st_diag", "Lst_diag", "1st_diag_date",
    "1st_stage", "Lst_stage", "max_stage", "1st_etiology", "Lst_etiology", "epilepsy_any",
    "converted", "conversion_to", "conversion_to_stage", "conversion_direction",
    "months_to_conversion",
    "gender", "education", "race", "age_at_baseline",
    "APOE_genotype", "APOE_e4_count",
]






_IMAGING_STEM = re.compile(r"^(MRI_|PET_|ABETA_|TAU_)")


def _best_set(events: dict, required: list, window_days: float,
              min_mod: int, prefer: str):
    """
    events : {label: [(t_days, row), ...] sorted by t}
    Tries every acquisition as the FIRST one of the set; for that start, each
    modality contributes its first acquisition inside the window (this gives
    the smallest span for that start). Returns (key, {label: (t, row)}) of the
    best valid set, or None.
    """
    times = {lab: [e[0] for e in ev] for lab, ev in events.items()}
    best_key, best = None, None
    for t0 in sorted({t for ts in times.values() for t in ts}):
        pick = {}
        for lab, ts in times.items():
            i = bisect_left(ts, t0)
            if i < len(ts) and ts[i] - t0 <= window_days:
                pick[lab] = events[lab][i]
        if len(pick) < min_mod or any(r not in pick for r in required):
            continue
        span = max(e[0] for e in pick.values()) - t0
        if prefer == "tightest":
            core = required or sorted(pick, key=lambda k: pick[k][0])[:min_mod]
            ct = [pick[k][0] for k in core]
            key = (max(ct) - min(ct), -len(pick), span, t0)
        else:
            key = (-len(pick), span, t0)
        if best_key is None or key < best_key:
            best_key, best = key, pick
    return (best_key, best) if best is not None else None


def filter_multimodal_patients(
    df: pd.DataFrame = None,
    window_months: float = None,
    modalities: list = None,
    add_match_cols: bool = True,
    min_modalities: int = None,
    prefer: str = None,
    input_path: str = None,
    keep_columns: list = None,
    databanks=None,
) -> pd.DataFrame:
    """
    One row per subject that has a valid set of modalities (see the header).

    df : expanded master already in memory (all its columns are returned), or
         None -> the needed columns are read from input_path / INPUT_PATH and
         only the identification + keep_columns are returned.

    Columns added (add_match_cols=True), for each modality <L> of the set:
      match_<L>_date        acquisition date (or day offset: OASIS / A4)
      match_<L>_id          image ID (imaging)      match_<L>_image_ids  all IDs of the session
      match_<L>_value       value (scored measures: CDR, plasma, CSF …)
      match_<L>_tracer      stem used, when the modality lists several
      match_<L>_manufacturer / _model / _suvr / _centiloid   when available
      match_<L>_offset_days days after the first acquisition of the set
    plus match_span_days, match_n_modalities, match_modalities, match_time_axis.
    """
    window_months = WINDOW_MONTHS if window_months is None else window_months
    modalities = MODALITIES if modalities is None else modalities
    min_modalities = MIN_MODALITIES if min_modalities is None else min_modalities
    prefer = PREFER if prefer is None else prefer
    keep_columns = KEEP_COLUMNS if keep_columns is None else keep_columns
    databanks = SELECT_DATABANKS if databanks is None else databanks
    if prefer not in ("most_modalities", "tightest"):
        raise ValueError("PREFER must be 'most_modalities' or 'tightest'.")

    labels = [m["label"] for m in modalities]
    if len(set(labels)) != len(labels):
        raise ValueError(f"Duplicated modality labels: {labels}")
    required = [m["label"] for m in modalities if m.get("required")]
    optional = [lab for lab in labels if lab not in required]
    min_mod = max(len(required), 1) if min_modalities is None else int(min_modalities)
    if min_mod > len(labels):
        raise ValueError(f"MIN_MODALITIES={min_mod} but only {len(labels)} modalities listed.")
    window_days = window_months * VRS_DAYS_PER_MONTH

    _print(f"Window     : {window_months} months ({window_days:.0f} days) first-to-last")
    _print(f"Required   : {required}")
    _print(f"Optional   : {optional}")
    _print(f"Minimum    : {min_mod} modalities in the set   |   prefer = '{prefer}'")

    lean = df is None
    if lean:
        df = vrs_read_expanded(input_path or INPUT_PATH,
                               [p for m in modalities for p in m["prefixes"]],
                               keep_cols=keep_columns, databanks=databanks)
    elif "uid" not in df.columns:
        df = df.assign(uid=df["Databank"].astype(str) + "_" + df["ID"].astype(str))

    # ── acquisitions of every modality ───────────────────────────────────────
    parts, multi = [], {}
    for m in modalities:
        ev = vrs_events(df, m["prefixes"], label=m["label"])
        multi[m["label"]] = len(vrs_measure_stems(m["prefixes"], df.columns)) > 1
        _print(f"  {m['label']:<14} {len(ev):>7,} acquisitions · {ev['uid'].nunique():>6,} subjects"
              + ("" if len(ev) else "   <- no column in the master for these prefixes"))
        parts.append(ev)
    E = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    missing_req = [lab for lab, p in zip(labels, parts) if lab in required and p.empty]
    if missing_req:
        raise ValueError(f"Required modalities without any acquisition: {missing_req}")
    E = E.sort_values(["uid", "kind", "t", "stem"]).reset_index(drop=True)

    # ── best set per subject (dates with dates, day offsets with day offsets) ──
    recs = E.to_dict("records")
    chosen = {}
    for (uid, kind), idx in E.groupby(["uid", "kind"], sort=False).indices.items():
        events = {}
        for i in idx:
            r = recs[i]
            lst = events.setdefault(r["label"], [])
            if not lst or lst[-1][0] != r["t"]:          # one per label and date
                lst.append((r["t"], r))
        got = _best_set(events, required, window_days, min_mod, prefer)
        if got and (uid not in chosen or got[0] < chosen[uid][0]):
            chosen[uid] = (got[0], got[1], kind)

    # ── output ───────────────────────────────────────────────────────────────
    rows = []
    for uid, (_, pick, kind) in chosen.items():
        t0 = min(e[0] for e in pick.values())
        row = {"uid": uid,
               "match_n_modalities": len(pick),
               "match_modalities": "+".join(lab for lab in labels if lab in pick),
               "match_span_days": int(round(max(e[0] for e in pick.values()) - t0)),
               "match_time_axis": kind}
        for lab in labels:
            if lab not in pick:
                continue
            t, r = pick[lab]
            row[f"match_{lab}_date"] = r["acq"]
            row[f"match_{lab}_id"] = r["img_id"]
            if multi[lab]:
                row[f"match_{lab}_tracer"] = r["stem"]
            if not _IMAGING_STEM.match(r["stem"]):
                row[f"match_{lab}_value"] = r["value"]
            for ch in ("image_ids", "manufacturer", "model", "suvr", "centiloid"):
                if r.get(ch) is not None:
                    row[f"match_{lab}_{ch}"] = r[ch]
            row[f"match_{lab}_offset_days"] = int(round(t - t0))
        rows.append(row)

    base = df if not lean else df[[c for c in df.columns
                                   if c in {"uid", "ID", "Databank", "dataset"} | set(keep_columns)]]
    result = base[base["uid"].isin(chosen)].copy()
    if add_match_cols and rows:
        M = pd.DataFrame(rows).set_index("uid")
        front = ["match_n_modalities", "match_modalities", "match_span_days", "match_time_axis"]
        order = front + [c for lab in labels for c in M.columns if c.startswith(f"match_{lab}_")]
        M = M[list(dict.fromkeys(order))].dropna(axis=1, how="all")   # e.g. no ID for plasma
        result = result.join(M, on="uid")

    n_total = df["uid"].nunique()
    _print(f"\n{'=' * 64}")
    _print(f"  Subjects in the master      : {n_total:,}")
    _print(f"  With a valid set            : {len(result):,} ({100 * len(result) / max(n_total, 1):.1f}%)")
    if len(result) and add_match_cols:
        bank = "Databank" if "Databank" in result.columns else None
        if bank:
            _print("\n  By databank:")
            _print(result[bank].value_counts().to_string())
        _print("\n  By combination:")
        _print(result["match_modalities"].value_counts().to_string())
        _print(f"\n  Span (days): median {result['match_span_days'].median():.0f} · "
              f"max {result['match_span_days'].max():.0f}")
        for lab in labels:
            if f"match_{lab}_tracer" in result.columns:
                _print(f"  {lab} stems: {result[f'match_{lab}_tracer'].value_counts().to_dict()}")
    _print(f"{'=' * 64}")
    return result


def multimodal_image_list(result: pd.DataFrame, modalities: list = None) -> pd.DataFrame:
    """Long list of the images of the sets: uid, Databank, ID, modality, stem,
    loniID, date, offset_days  (one row per image)."""
    modalities = MODALITIES if modalities is None else modalities
    out = []
    for m in modalities:
        lab = m["label"]
        idc = f"match_{lab}_id"
        if idc not in result.columns:
            continue
        sub = result[result[idc].notna()]
        out.append(pd.DataFrame({
            "uid": sub["uid"],
            "Databank": sub["Databank"] if "Databank" in sub.columns else np.nan,
            "ID": sub["ID"] if "ID" in sub.columns else np.nan,
            "modality": lab,
            "stem": sub[f"match_{lab}_tracer"] if f"match_{lab}_tracer" in sub.columns else lab,
            "loniID": sub[idc].astype(str),
            "date": sub[f"match_{lab}_date"],
            "offset_days": sub[f"match_{lab}_offset_days"],
        }))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# name kept so older code that calls the tri-modal function keeps working
filter_trimodal_patients = filter_multimodal_patients





@quiet_pandas_warnings
def create_multimodal_cohort(input_path: str | os.PathLike, *,
                             modalities: list | None = None,
                             window_months: float = 12,
                             min_modalities: int | None = None,
                             prefer: str = "most_modalities",
                             keep_columns: list | None = None,
                             databanks="ALL",
                             output_prefix: str = "multimodal",
                             out_dir: str | os.PathLike | None = None,
                             config: ProjectConfig | None = None) -> CohortResult:
    """Cross-sectional multimodal cohort: one set of modalities per subject.

    Parameters
    ----------
    input_path:
        ``MASTER_expanded_240m.csv``.
    modalities:
        List of ``{"label": ..., "prefixes": [...], "required": bool}``
        (default :data:`DEFAULT_MODALITIES`: T1w and FDG required; fMRI,
        amyloid and tau PET optional).
    window_months:
        Maximum span between the first and the last acquisition of the set.
    min_modalities:
        ``None`` -> every required modality; ``k`` -> at least ``k`` modalities,
        including all required ones (with none required: any ``k``).
    prefer:
        ``"most_modalities"`` or ``"tightest"``.
    keep_columns:
        Subject-level columns copied to the output.
    databanks:
        ``"ALL"`` or a list of datasets.
    output_prefix, out_dir:
        Writes ``<prefix>_multimodal.csv``, ``<prefix>_multimodal_loniIDs.csv``
        and ``<prefix>_multimodal_loniIDs.txt`` (default folder: that of
        ``input_path``).
    config:
        Project configuration (measure groups, time-axis units).

    Returns
    -------
    CohortResult
    """
    cfg = config or load_config()
    configure_common(cfg)
    modalities = DEFAULT_MODALITIES if modalities is None else modalities
    result = filter_multimodal_patients(window_months=window_months, modalities=modalities,
                                        min_modalities=min_modalities, prefer=prefer,
                                        input_path=str(input_path), keep_columns=keep_columns,
                                        databanks=databanks)
    out_dir = str(out_dir or os.path.dirname(str(input_path)))
    os.makedirs(out_dir, exist_ok=True)
    files = {"table": os.path.join(out_dir, f"{output_prefix}_multimodal.csv")}
    result.to_csv(files["table"], index=False)
    _print(f"Saved -> {files['table']}")
    imgs = multimodal_image_list(result, modalities)
    if not imgs.empty:
        files["images"] = os.path.join(out_dir, f"{output_prefix}_multimodal_loniIDs.csv")
        files["image_ids"] = os.path.join(out_dir, f"{output_prefix}_multimodal_loniIDs.txt")
        imgs.to_csv(files["images"], index=False)
        ids = imgs["loniID"].str.strip().drop_duplicates().tolist()
        with open(files["image_ids"], "w") as fh:
            fh.write(",".join(ids))
        _print(f"Saved -> {files['images']}")
        _print(f"Saved -> {files['image_ids']}  ({len(ids):,} unique image IDs)")
    write_provenance(files["table"], step="multimodal_cohort", config=cfg, inputs=[str(input_path)],
                     table=result, parameters={"modalities": modalities, "window_months": window_months,
                                               "min_modalities": min_modalities, "prefer": prefer,
                                               "databanks": databanks})
    return CohortResult(table=result, images=imgs, files=files)
