"""Longitudinal visits: subjects with repeated acquisitions of ONE modality at
defined intervals counted from their own earliest acquisition.

Per subject:

1. every acquisition of the modality is placed on the subject's time axis
   (calendar dates, or days since entry for OASIS / A4);
2. acquisitions close in time are ONE session / visit (``session_window_days``);
3. targets = baseline + ``visit_targets_months``; each target takes the session
   closest to it inside +/- ``tolerance_months``. The assignment is globally
   optimal (Hungarian algorithm; greedy fallback without scipy) and the
   tolerance never lets two targets overlap (``clamp_windows``).

The pipelines write ONE entry per session (date, a representative image ID and
every image ID of the session), so dates and IDs are always paired. The output
also reports the scanner of each visit and a ``<modality>_scanner_change`` flag.
"""
from __future__ import annotations

import numpy as np
import os
import pandas as pd
import re

from ..config import ProjectConfig, load_config
from ..integration.common import configure as configure_common
from ..provenance import write_provenance
from .outputs import CohortResult

from ..integration.common import VRS_DAYS_PER_MONTH, vrs_clean_id_list, vrs_events, vrs_read_expanded

import logging as _logging
from ..logging_utils import make_print as _make_print
from ..logging_utils import quiet_pandas_warnings

log = _logging.getLogger(__name__)
_print = _make_print(log)

# ============================================================
#  DEFAULTS (override them through the function arguments)
# ============================================================
INPUT_PATH: str | None = None   # set through the input_path argument

SELECT_DATABANKS = "ALL"      # 'ALL' or ['ADNI', 'OASIS', ...]
# optional: CSV with a 'uid' column restricting the subjects to an existing cohort (or None)
COHORT_FILTER_PATH = None

MODALITY_LABEL    = "MRI_T1w"
# exact stems and/or groups ("PET_ABETA", "PET_TAU" ...)
MODALITY_PREFIXES = ["MRI_T1w"]

# Targets in months from the subject's own baseline.
VISIT_TARGETS_MONTHS = [0, 12, 24]

# Tolerance (+/-) in months around each target.
TOLERANCE_MONTHS = 6

# Minimum number of matched visits to keep a subject.
MIN_VISITS_REQUIRED = 2

# -- Session grouping -----------------------------------------
# Acquisitions within this many days of the first one of a session are
# the SAME visit.
#   0  -> same exact date only
#   31 -> scans within ~1 month count as one session
SESSION_WINDOW_DAYS = 31

# -- Which series to use when a session has several -------------
#   "first"  -> representative ID of the first date of the session
#   "last"   -> representative ID of the last date of the session
#   "lowest_id" / "highest_id" -> by number, among ALL IDs of the session
SERIES_PICK = "first"

# Keep the tolerance from creating overlapping windows (recommended).
CLAMP_OVERLAPPING_WINDOWS = True

# Subject-level columns copied to the output (when present in the master)
KEEP_COLUMNS = [
    "Databank", "ID",
    "1st_diag", "Lst_diag", "1st_diag_date",
    "gender", "education", "race", "age_at_baseline",
    "APOE_genotype", "APOE_e4_count",
]






try:
    from scipy.optimize import linear_sum_assignment
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

_BIG = 10 ** 9


# ---------------------------------------------------------------------------
# Agrupamento em sessoes
# ---------------------------------------------------------------------------

def _id_num(img_id):
    """Number of an image ID like 'I24864' / '24864' (for sorting)."""
    if img_id is None:
        return float("inf")
    m = re.search(r"(\d+)", str(img_id))
    return int(m.group(1)) if m else float("inf")


def _group_sessions(events: list, window_days: float, pick: str) -> list:
    """
    events: records (dicts from vrs_events) of ONE subject, ordered by t.
    Groups into sessions the records within window_days of the session's first one.
    """
    if not events:
        return []
    groups, cur, anchor = [], [events[0]], events[0]["t"]
    for e in events[1:]:
        if e["t"] - anchor <= window_days:
            cur.append(e)
        else:
            groups.append(cur)
            cur, anchor = [e], e["t"]
    groups.append(cur)

    out = []
    for grp in groups:
        all_ids = vrs_clean_id_list(
            [x for e in grp for x in (e.get("image_ids") or e.get("img_id"),)])
        if pick == "last":
            chosen, cid = grp[-1], grp[-1].get("img_id")
        elif pick in ("lowest_id", "highest_id") and all_ids:
            cid = (min if pick == "lowest_id" else max)(all_ids, key=_id_num)
            chosen = next((e for e in grp if cid in vrs_clean_id_list(
                [e.get("image_ids") or e.get("img_id")])), grp[0])
        else:  # "first"
            chosen, cid = grp[0], grp[0].get("img_id")
        out.append({
            "t": chosen["t"], "acq": chosen["acq"], "id": cid,
            "all_ids": all_ids, "n_series": max(len(all_ids), len(grp)),
            "span_days": grp[-1]["t"] - grp[0]["t"],
            "manufacturer": chosen.get("manufacturer"), "model": chosen.get("model"),
            "stem": chosen.get("stem"),
        })
    return out


# ---------------------------------------------------------------------------
# Targets and tolerances
# ---------------------------------------------------------------------------

def _target_offsets_days(baseline, targets_months, kind):
    """Offsets in days. Dates: real calendar months; days since entry:
    30.4375 days per month."""
    if kind == "date":
        b = pd.Timestamp("1970-01-01") + pd.Timedelta(days=float(baseline["t"]))
        return [int(((b + pd.DateOffset(months=int(m))) - b).days) for m in targets_months]
    return [int(round(m * VRS_DAYS_PER_MONTH)) for m in targets_months]


def _tolerances_days(target_days, tolerance_days, clamp: bool):
    """Tolerance per target, optionally reduced so that neighbours do not overlap."""
    tols = []
    for i, t in enumerate(target_days):
        tol = tolerance_days
        if clamp:
            if i > 0:
                tol = min(tol, (t - target_days[i - 1]) // 2)
            if i < len(target_days) - 1:
                tol = min(tol, (target_days[i + 1] - t) // 2)
        tols.append(int(tol))
    return tols


# ---------------------------------------------------------------------------
# Globally optimal matching
# ---------------------------------------------------------------------------

def _match_visits(sessions, targets_months, tolerance_days, clamp, kind):
    if not sessions:
        return [], []
    base = sessions[0]
    target_days = _target_offsets_days(base, targets_months, kind)
    tols = _tolerances_days(target_days, tolerance_days, clamp)
    offsets = [s["t"] - base["t"] for s in sessions]

    n_t, n_s = len(target_days), len(sessions)
    cost = np.full((n_t, n_s), float(_BIG))
    for i in range(n_t):
        for j in range(n_s):
            d = abs(offsets[j] - target_days[i])
            if d <= tols[i]:
                cost[i, j] = d

    assign = {}
    if _HAS_SCIPY:
        rows, cols = linear_sum_assignment(cost)
        for i, j in zip(rows, cols):
            if cost[i, j] < _BIG:
                assign[i] = j
    else:   # fallback: globally cheapest pairs first
        pairs = sorted((cost[i, j], i, j) for i in range(n_t) for j in range(n_s)
                       if cost[i, j] < _BIG)
        used_t, used_s = set(), set()
        for _, i, j in pairs:
            if i in used_t or j in used_s:
                continue
            assign[i] = j
            used_t.add(i)
            used_s.add(j)

    matches = []
    for i, m in enumerate(targets_months):
        if i in assign:
            s = sessions[assign[i]]
            matches.append({
                "target_months": int(m), "date": s["acq"], "id": s["id"],
                "all_ids": ";".join(s["all_ids"]) if s["all_ids"] else None,
                "n_series": s["n_series"],
                "offset_days": int(round(s["t"] - base["t"])),
                "manufacturer": s["manufacturer"], "model": s["model"],
                "matched": True,
            })
        else:
            matches.append({"target_months": int(m), "date": None, "id": None,
                            "all_ids": None, "n_series": 0, "offset_days": None,
                            "manufacturer": None, "model": None, "matched": False})
    used = set(assign.values())
    return matches, [sessions[j] for j in range(n_s) if j not in used]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def filter_longitudinal_visits(
    df: pd.DataFrame = None,
    modality_label:    str  = None,
    modality_prefixes: list = None,
    visit_targets:     list = None,
    tolerance_months:  float = None,
    min_visits:        int  = None,
    session_window_days: float = None,
    series_pick:       str  = None,
    clamp_windows:     bool = None,
    input_path:        str  = None,
    keep_columns:      list = None,
    databanks=None,
    cohort_filter_path: str | None = None,
) -> pd.DataFrame:
    """One row per kept subject, with columns per visit.

    df : expanded master already in memory (all its columns are returned), or
         None -> only the needed columns are read from input_path and the
         identification + keep_columns are returned.

    Columns per visit n (label = {modality}_V{n}_{target}m):
        {label}_date, {label}_id, {label}_all_ids, {label}_n_series,
        {label}_offset_days, {label}_manufacturer, {label}_model, {label}_matched
    Plus:
        {modality}_baseline_date, {modality}_n_matched, {modality}_n_sessions,
        {modality}_n_series_total, {modality}_n_sessions_unused,
        {modality}_time_axis ('date' | 'visit' = days since entry),
        {modality}_scanner_change (True when vendor / model changes between the matched visits)
    """
    g = globals()
    modality_label    = g["MODALITY_LABEL"] if modality_label is None else modality_label
    modality_prefixes = g["MODALITY_PREFIXES"] if modality_prefixes is None else modality_prefixes
    visit_targets     = g["VISIT_TARGETS_MONTHS"] if visit_targets is None else visit_targets
    tolerance_months  = g["TOLERANCE_MONTHS"] if tolerance_months is None else tolerance_months
    min_visits        = g["MIN_VISITS_REQUIRED"] if min_visits is None else min_visits
    session_window_days = g["SESSION_WINDOW_DAYS"] if session_window_days is None else session_window_days
    series_pick       = g["SERIES_PICK"] if series_pick is None else series_pick
    clamp_windows     = g["CLAMP_OVERLAPPING_WINDOWS"] if clamp_windows is None else clamp_windows
    keep_columns      = g["KEEP_COLUMNS"] if keep_columns is None else keep_columns
    databanks         = g["SELECT_DATABANKS"] if databanks is None else databanks
    if series_pick not in ("first", "last", "lowest_id", "highest_id"):
        raise ValueError("series_pick must be 'first', 'last', 'lowest_id' or 'highest_id'.")

    tolerance_days = int(round(tolerance_months * VRS_DAYS_PER_MONTH))
    n_targets = len(visit_targets)
    mod = modality_label

    _print(f"Modality        : {mod}  (stems/groups: {modality_prefixes})")
    _print(f"Schedule        : {visit_targets} months  +/-{tolerance_months} m")
    _print(f"Session window  : {session_window_days} days   | series pick: {series_pick}")
    _print(f"Min visits      : {min_visits} / {n_targets}")
    if not _HAS_SCIPY:
        _print("[WARN] scipy not found -> using the global greedy fallback.")

    lean = df is None
    if lean:
        df = vrs_read_expanded(input_path or INPUT_PATH, modality_prefixes,
                               keep_cols=keep_columns, databanks=databanks)
    elif "uid" not in df.columns:
        df = df.assign(uid=df["Databank"].astype(str) + "_" + df["ID"].astype(str))
    if cohort_filter_path:
        keep = set(pd.read_csv(cohort_filter_path, usecols=["uid"], dtype=str)["uid"])
        df = df[df["uid"].isin(keep)]
        _print(f"  cohort filter: {len(df):,} subjects")

    _print("Extracting acquisitions ...")
    E = vrs_events(df, modality_prefixes, label=mod)
    if E.empty:
        raise ValueError(f"No acquisition found for {modality_prefixes} "
                         f"(check the stem names in the expanded master).")
    E = E.sort_values(["uid", "kind", "t"]).reset_index(drop=True)
    recs = E.to_dict("records")

    rows, collapsed, n_series_tot = {}, 0, 0
    for (uid, kind), idx in E.groupby(["uid", "kind"], sort=False).indices.items():
        events = [recs[i] for i in idx]
        sessions = _group_sessions(events, session_window_days, series_pick)
        matches, unused = _match_visits(sessions, visit_targets, tolerance_days,
                                        clamp_windows, kind)
        n_matched = sum(m["matched"] for m in matches)
        if n_matched < min_visits:
            continue
        if uid in rows and rows[uid][f"{mod}_n_matched"] >= n_matched:
            continue                               # another time axis already gave more visits
        collapsed += len(events) - len(sessions)
        scanners = {(m["manufacturer"], m["model"]) for m in matches
                    if m["matched"] and (m["manufacturer"] or m["model"])}
        row = {
            f"{mod}_baseline_date":     sessions[0]["acq"],
            f"{mod}_n_matched":         n_matched,
            f"{mod}_n_sessions":        len(sessions),
            f"{mod}_n_series_total":    sum(s["n_series"] for s in sessions),
            f"{mod}_n_sessions_unused": len(unused),
            f"{mod}_time_axis":         kind,
            f"{mod}_scanner_change":    (len(scanners) > 1) if scanners else np.nan,
        }
        for v_idx, m in enumerate(matches, start=1):
            lbl = f"{mod}_V{v_idx}_{m['target_months']}m"
            for k in ("date", "id", "all_ids", "n_series", "offset_days",
                      "manufacturer", "model", "matched"):
                row[f"{lbl}_{k}"] = m[k]
        rows[uid] = row

    base = df if not lean else df[[c for c in df.columns
                                   if c in {"uid", "ID", "Databank", "dataset"} | set(keep_columns)]]
    result = base[base["uid"].isin(rows)].copy()
    if rows:
        R = pd.DataFrame.from_dict(rows, orient="index")
        R = R.drop(columns=[c for c in R.columns
                            if c.endswith(("_manufacturer", "_model")) and R[c].isna().all()])
        result = result.join(R, on="uid")

    n_total, n_filtered = df["uid"].nunique(), len(result)
    _print(f"\n{'=' * 60}")
    _print(f"  Subjects in total          : {n_total:,}")
    _print(f"  With the modality          : {E['uid'].nunique():,}")
    _print(f"  Dates merged into sessions : {collapsed:,}")
    _print(f"  With >={min_visits} matched visits   : "
          f"{n_filtered:,} ({100 * n_filtered / max(n_total, 1):.1f}%)")
    if n_filtered:
        for v_idx, t in enumerate(visit_targets, start=1):
            col = f"{mod}_V{v_idx}_{t}m_matched"
            if col in result.columns:
                n = int(result[col].fillna(False).astype(bool).sum())
                _print(f"  V{v_idx} ({t:>3} m)  : {n:,} subjects ({100 * n / n_filtered:.1f}%)")
        if "Databank" in result.columns:
            _print("\n  Per dataset:")
            _print(result["Databank"].value_counts().to_string())
        sc = result[f"{mod}_scanner_change"]
        if sc.notna().any():
            _print(f"\n  Scanner change between the matched visits: "
                  f"{int((sc == True).sum()):,} of {int(sc.notna().sum()):,} subjects with a known scanner")
    _print(f"{'=' * 60}")
    return result


def longitudinal_image_list(result: pd.DataFrame, modality_label: str = None) -> pd.DataFrame:
    """Long list of the images of the matched visits: uid, Databank, ID, visit,
    loniID, date, offset_days."""
    mod = MODALITY_LABEL if modality_label is None else modality_label
    out = []
    for c in [c for c in result.columns if re.fullmatch(re.escape(mod) + r"_V\d+_\d+m_id", c)]:
        lbl = c[:-3]
        sub = result[result[c].notna()]
        out.append(pd.DataFrame({
            "uid": sub["uid"],
            "Databank": sub["Databank"] if "Databank" in sub.columns else np.nan,
            "ID": sub["ID"] if "ID" in sub.columns else np.nan,
            "visit": lbl[len(mod) + 1:],
            "loniID": sub[c].astype(str),
            "date": sub[f"{lbl}_date"],
            "offset_days": sub[f"{lbl}_offset_days"],
        }))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()





@quiet_pandas_warnings
def create_longitudinal_cohort(input_path: str | os.PathLike, *,
                               modality_label: str = "MRI_T1w",
                               modality_prefixes: list | None = None,
                               visit_targets_months: list = (0, 12, 24),
                               tolerance_months: float = 6,
                               min_visits: int = 2,
                               session_window_days: float = 31,
                               series_pick: str = "first",
                               clamp_windows: bool = True,
                               keep_columns: list | None = None,
                               databanks="ALL",
                               cohort_filter_path: str | os.PathLike | None = None,
                               output_prefix: str = "longitudinal",
                               out_dir: str | os.PathLike | None = None,
                               config: ProjectConfig | None = None) -> CohortResult:
    """Longitudinal cohort of one modality (stem or group).

    Parameters
    ----------
    input_path:
        ``MASTER_expanded_240m.csv``.
    modality_label, modality_prefixes:
        Output label and the stems / groups that count (default: the label).
    visit_targets_months:
        Targets in months from the subject's first session.
    tolerance_months:
        +/- months around each target.
    min_visits:
        Minimum number of matched targets to keep a subject.
    session_window_days:
        Acquisitions within this many days of the first one of a session are
        the same visit (0 = same date only).
    series_pick:
        ``first`` / ``last`` (representative ID of the first / last date of the
        session) or ``lowest_id`` / ``highest_id`` (among ALL IDs of the session).
    clamp_windows:
        Shrink tolerances so that neighbouring targets never overlap.
    keep_columns, databanks, cohort_filter_path:
        Subject-level columns, datasets and an optional CSV with a ``uid``
        column restricting the subjects.
    output_prefix, out_dir:
        Writes ``<prefix>_longitudinal.csv``, ``<prefix>_longitudinal_loniIDs.csv``
        and ``.txt`` (default folder: that of ``input_path``).
    config:
        Project configuration.

    Returns
    -------
    CohortResult
    """
    cfg = config or load_config()
    configure_common(cfg)
    result = filter_longitudinal_visits(
        modality_label=modality_label,
        modality_prefixes=list(modality_prefixes or [modality_label]),
        visit_targets=list(visit_targets_months), tolerance_months=tolerance_months,
        min_visits=min_visits, session_window_days=session_window_days,
        series_pick=series_pick, clamp_windows=clamp_windows, input_path=str(input_path),
        keep_columns=keep_columns, databanks=databanks, cohort_filter_path=cohort_filter_path)
    out_dir = str(out_dir or os.path.dirname(str(input_path)))
    os.makedirs(out_dir, exist_ok=True)
    files = {"table": os.path.join(out_dir, f"{output_prefix}_longitudinal.csv")}
    result.to_csv(files["table"], index=False)
    _print(f"\nSaved -> {files['table']}")
    imgs = longitudinal_image_list(result, modality_label)
    if not imgs.empty:
        files["images"] = os.path.join(out_dir, f"{output_prefix}_longitudinal_loniIDs.csv")
        files["image_ids"] = os.path.join(out_dir, f"{output_prefix}_longitudinal_loniIDs.txt")
        imgs.to_csv(files["images"], index=False)
        ids = imgs["loniID"].str.strip().drop_duplicates().tolist()
        with open(files["image_ids"], "w") as fh:
            fh.write(",".join(ids))
        _print(f"Saved -> {files['images']}")
        _print(f"Saved -> {files['image_ids']}  ({len(ids):,} unique image IDs)")
    write_provenance(files["table"], step="longitudinal_cohort", config=cfg, inputs=[str(input_path)],
                     table=result, parameters={"modality_label": modality_label,
                                               "modality_prefixes": modality_prefixes,
                                               "visit_targets_months": list(visit_targets_months),
                                               "tolerance_months": tolerance_months,
                                               "min_visits": min_visits,
                                               "session_window_days": session_window_days,
                                               "series_pick": series_pick, "databanks": databanks,
                                               "scipy": _HAS_SCIPY})
    return CohortResult(table=result, images=imgs, files=files)
