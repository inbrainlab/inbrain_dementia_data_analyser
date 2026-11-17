"""Anchor pairing on the expanded master (the "cohort creator").
"""
from __future__ import annotations

import ast
import numpy as np
import os
import pandas as pd
import re
from typing import Sequence

from ..config import ProjectConfig, load_config
from ..integration.common import configure as configure_common
from ..provenance import write_provenance
from .outputs import CohortResult

from ..config import ProjectConfig, load_config

import logging as _logging
from ..logging_utils import make_print as _make_print
from ..logging_utils import quiet_pandas_warnings

log = _logging.getLogger(__name__)
_print = _make_print(log)

DEFAULT_KEEP_COLUMNS = (
    "Databank", "ID",
    "1st_diag", "Lst_diag", "1st_diag_date", "Lst_diag_date",
    "months_to_conversion", "conversion_to", "converted",
    "gender", "education", "race", "hispanic", "age_at_baseline",
    "APOE_genotype", "APOE_e4_count", "APOE_e4_any",
    # seizures / epilepsy - subject level (history, not paired in time)
    "SEIZURE_any", "SEIZURE_source", "SEIZURES_ever", "EPILEP_ever",
    "Antiepileptic_med_ever", "SEIZURE_text_flag", "SEIZURE_text_evidence",
    "AED_text_flag", "AED_text_evidence",
)


@quiet_pandas_warnings
def create_anchor_cohort(input_path: str | os.PathLike, *,
                         anchor: str = "MRI_T1w",
                         required: Sequence[str] = ("CDR",),
                         optional: Sequence[str] = (),
                         margin_months: float = 6,
                         margin_months_by_feature: dict | None = None,
                         anchor_events: str = "all",
                         allow_image_reuse: bool = True,
                         include_pre_diagnosis: bool = True,
                         max_months: int = 240,
                         keep_columns: Sequence[str] = DEFAULT_KEEP_COLUMNS,
                         required_subject_columns: Sequence[str] = (),
                         select_databanks: str | Sequence[str] = "ALL",
                         cohort_filter_path: str | os.PathLike | None = None,
                         output_prefix: str = "anchor_cohort",
                         out_dir: str | os.PathLike | None = None,
                         summary_include_delta: bool = True,
                         save_xlsx: bool = False,
                         verbose_unparsed: bool = True,
                         merge_duplicate_uids: bool = True,
                         config: ProjectConfig | None = None) -> CohortResult:
    """Anchor pairing: each acquisition of ``anchor`` with the closest acquisition
    of every other measure within ``margin_months``.

    Any stem of the expanded master can be used (``MRI_T1w``, ``MRI_fMRI``,
    ``MMSE``, ``CDR_SB``, ``Plasma_pTau217_C2N`` ...), and also the groups of
    ``imaging.yaml`` (``PET_ABETA`` = any amyloid tracer, ``PET_TAU``) and the
    fluid groups (``PLASMA_PTAU217``, ``CSF_AB42`` ...).

    Parameters
    ----------
    input_path:
        ``MASTER_expanded_240m.csv``.
    anchor:
        Measure (stem or group) whose acquisitions define the events.
    required:
        Measures that must exist within the margin, otherwise the event is dropped.
    optional:
        Measures added when they exist within the margin.
    margin_months:
        +- months around the anchor (0 = same day / visit).
    margin_months_by_feature:
        Margin per measure, overriding ``margin_months`` (e.g. ``{"BMI": 12}``).
    anchor_events:
        Events kept per subject: ``all`` (longitudinal), ``first``, ``last`` or
        ``most_optional`` (most optional matches; ties -> earliest). Use
        ``first`` / ``most_optional`` for a cross-sectional cohort.
    allow_image_reuse:
        Allow one image to be matched to two anchor events.
    include_pre_diagnosis, max_months:
        Windows used (before the first diagnosis; horizon after it).
    keep_columns:
        Subject-level columns copied to the output.
    required_subject_columns:
        Subject-level columns that must be filled (e.g. ``["SEIZURE_any"]``).
    select_databanks:
        ``"ALL"`` or a list of datasets.
    cohort_filter_path:
        CSV with a ``uid`` column restricting the subjects.
    output_prefix, out_dir:
        Output names / folder (default: the folder of ``input_path``).
    summary_include_delta, save_xlsx, verbose_unparsed, merge_duplicate_uids:
        Output and diagnostics options.
    config:
        Project configuration (measure groups, time-axis units).

    Returns
    -------
    CohortResult
        ``table`` = one row per anchor event (``<prefix>_pairs.csv``),
        ``images`` = long image list (``<prefix>_loniIDs.csv``), ``summary``
        (``<prefix>_summary.csv``) and the written ``files``. The unique image
        IDs are written comma-separated to ``<prefix>_loniIDs.txt``.

    Notes
    -----
    Diagnosis at the event: the diagnosis of the anchor window; if empty, the
    last diagnosis recorded up to the anchor date; otherwise ``Lst_diag``
    (column ``diag_source``). Untimed acquisitions (no date / visit, part of
    NACC) cannot be paired in time and are ignored.
    """
    cfg = config or load_config()
    configure_common(cfg)
    INPUT_PATH = str(input_path)
    OUT_DIR = str(out_dir) if out_dir else None
    OUTPUT_PREFIX = output_prefix
    COHORT_FILTER_PATH = str(cohort_filter_path) if cohort_filter_path else None
    SELECT_DATABANKS = select_databanks if select_databanks == "ALL" else list(select_databanks)
    ANCHOR = anchor
    MARGIN_MONTHS = margin_months
    MARGIN_MONTHS_BY_FEATURE = dict(margin_months_by_feature or {})
    REQUIRED = list(required)
    OPTIONAL = list(optional)
    ANCHOR_EVENTS = anchor_events
    ALLOW_IMAGE_REUSE = allow_image_reuse
    INCLUDE_PRE_DIAGNOSIS = include_pre_diagnosis
    MAX_MONTHS = max_months
    KEEP_COLUMNS = list(keep_columns)
    REQUIRED_SUBJECT_COLUMNS = list(required_subject_columns)
    MEASURE_GROUPS = {k: list(v) for k, v in cfg.imaging.get("measure_groups", {}).items()}
    FLUID_GROUP_PATTERNS = dict(cfg.imaging.get("fluid_group_patterns", {}))
    SUMMARY_INCLUDE_DELTA = summary_include_delta
    SAVE_XLSX = save_xlsx
    VERBOSE_UNPARSED = verbose_unparsed
    ta = cfg.integration.get("time_axis", {})
    NUMERIC_ACQ_UNIT = dict(ta.get("numeric_acquisition_unit", {}))
    ID_DAYCODE_BANKS = list(ta.get("id_daycode_banks", ["OASIS"]))
    MERGE_DUPLICATE_UIDS = merge_duplicate_uids
    DAYS_PER_MONTH = float(ta.get("days_per_month", 30.4375))
    SCANNER_CHANNELS = ["manufacturer", "model"]
    VALUE_CHANNELS = ["suvr", "centiloid"]
    ID_CHANNELS = ["image_ids"]
    EXTRA_CHANNELS = SCANNER_CHANNELS + VALUE_CHANNELS + ID_CHANNELS
    FLUID_META = ["kit", "lab"]
    # =============================================================================
    # 2. LOAD
    # =============================================================================
    if not os.path.exists(INPUT_PATH):
        raise FileNotFoundError(f"Expanded master not found:\n  {INPUT_PATH}")

    header = pd.read_csv(INPUT_PATH, nrows=0).columns

    # build fluid groups from the stems that actually exist in the expanded master
    _stems = sorted({c[:-4].rsplit("_", 1)[0] for c in header if c.endswith("_acq")})
    for _g, _rx in FLUID_GROUP_PATTERNS.items():
        _members = [s_ for s_ in _stems if re.match(_rx, s_)]
        if _members and _g not in MEASURE_GROUPS:
            MEASURE_GROUPS[_g] = _members
    # -----------------------------------------------------------------------------
    # Reads ONLY the columns in use. The expanded master has thousands of columns
    # (each measure x ~80 windows x value/_acq/_id/...); loading everything
    # exhausts the Colab RAM. Read here:
    #   identifiers + KEEP_COLUMNS + REQUIRED_SUBJECT_COLUMNS + age / diagnosis columns
    #   diag_<window>(_acq) and <stem>_<window>(_acq|_id) of the anchor / required / optional
    # -----------------------------------------------------------------------------
    _needed_stems = set()
    for _f in [ANCHOR] + list(REQUIRED) + list(OPTIONAL):
        _needed_stems.update(MEASURE_GROUPS.get(_f, [_f]))
    _base_cols = {"uid", "ID", "Databank", "dataset", "cohort", "source",
                  "Lst_diag", "1st_diag_date", "age_at_1st_diag", "age_at_baseline",
                  f"{ANCHOR}_untimed_n"} | set(KEEP_COLUMNS) | set(REQUIRED_SUBJECT_COLUMNS)
    _win_tail = re.compile(r"_(pre\d*|pregt\d+|lt\d+|gt\d+)(_acq|_id|_"
                           + "|_".join(EXTRA_CHANNELS) + r")?$")
    # laboratory / kit of the fluid stems in play (subject-level columns)
    _base_cols |= {f"{s_}_{k_}" for s_ in _needed_stems for k_ in FLUID_META}
    usecols = []
    for c in header:
        if c in _base_cols:
            usecols.append(c)
            continue
        m_ = _win_tail.search(c)
        if m_ and (c[:m_.start()] in _needed_stems or c[:m_.start()] == "diag"):
            usecols.append(c)
    _print(f"Reading {len(usecols)} of {len(header)} columns of the expanded master "
          f"(stems: {sorted(_needed_stems)})")

    # IDs and dates as text (keeps loniID from becoming 123456.0)
    # Read in row chunks: pd.read_csv(low_memory=False) tokenises the WHOLE file
    # (every column) before applying usecols - on the expanded master that costs
    # several GB even for 50 columns. In chunks, peak RAM follows the columns used.
    str_cols = {c for c in usecols
                if c.endswith(("_id", "_acq", "_image_ids")) or c in ("uid", "ID")}
    df = pd.concat(pd.read_csv(INPUT_PATH, usecols=usecols, dtype=str, chunksize=4000),
                   ignore_index=True)
    for _c in df.columns:                 # same result as type inference on the whole file
        if _c not in str_cols:
            _num = pd.to_numeric(df[_c], errors="coerce")
            if _num.notna().sum() == df[_c].notna().sum():
                df[_c] = _num
    out_dir = OUT_DIR or os.path.dirname(INPUT_PATH)

    COHORT_COL = next((c for c in ["Databank", "dataset", "cohort", "source"]
                       if c in df.columns), None)
    if "uid" not in df.columns:
        df["uid"] = df[COHORT_COL].astype(str) + "_" + df["ID"].astype(str)

    if SELECT_DATABANKS != "ALL":
        df = df[df[COHORT_COL].isin(SELECT_DATABANKS)]
    if COHORT_FILTER_PATH:
        keep_uids = set(pd.read_csv(COHORT_FILTER_PATH, usecols=["uid"], dtype=str)["uid"])
        df = df[df["uid"].isin(keep_uids)]
    for _c in REQUIRED_SUBJECT_COLUMNS:
        if _c not in df.columns:
            raise ValueError(f"required_subject_columns: '{_c}' is not in the expanded master "
                             f"(re-run the expansion).")
        _before = df["uid"].nunique()
        df = df[df[_c].notna() & (df[_c].astype(str).str.strip() != "")]
        _print(f"  subject filter '{_c}': {_before} -> {df['uid'].nunique()} subjects")
    n_rows = len(df)
    dup = df["uid"].duplicated(keep=False)
    if dup.any():
        how = "rows merged (first non-empty value per column)" if MERGE_DUPLICATE_UIDS \
            else "first row of each kept"
        _print(f"  [WARN] {int(dup.sum())} rows with a repeated uid "
              f"({df.loc[dup, 'uid'].nunique()} uids) -> {how}. "
              f"Per dataset: {df.loc[dup, COHORT_COL].value_counts().to_dict()}")
        if MERGE_DUPLICATE_UIDS:
            merged = df[dup].groupby("uid", sort=False).first().reset_index()
            df = pd.concat([df[~dup], merged[df.columns]], ignore_index=True)
    df = df.drop_duplicates("uid").reset_index(drop=True)

    _print("=" * 70)
    _print("ANCHOR PAIRING")
    _print("=" * 70)
    _print(f"Rows in master: {n_rows}  |  subjects (unique uid): {len(df)}")
    _print(f"Anchor: {ANCHOR}  |  margin: ±{MARGIN_MONTHS} months  |  ANCHOR_EVENTS = '{ANCHOR_EVENTS}'")
    _print(f"Required: {REQUIRED}")
    _print(f"Optional: {OPTIONAL}")

    # =============================================================================
    # 3. HELPERS
    # =============================================================================
    WIN_RE = r"(pre\d*|pregt\d+|lt\d+|gt\d+)"


    def window_ok(w):
        if w.startswith("pre"):
            return INCLUDE_PRE_DIAGNOSIS
        if w.startswith("lt"):
            return int(w[2:]) <= max(MAX_MONTHS, 6)
        if w.startswith("gt"):
            return MAX_MONTHS >= 240
        return False


    def windows_of(stem):
        pat = re.compile(re.escape(stem) + "_" + WIN_RE + "$")
        return [m.group(1) for c in df.columns if (m := pat.fullmatch(c)) and window_ok(m.group(1))]


    DATE_LIKE = (r"\d{4}[-/.]\d{1,2}[-/.]\d{1,2}|\d{1,2}[-/.]\d{1,2}[-/.]\d{2,4}"
                 r"|^(?:19|20)\d{6}$|^(?:19|20)\d{2}[-/]\d{1,2}$")
    # visit codes relative to baseline (bl = 0)
    VISIT_RE = (r"^(?:(?P<bl>bl|sc|scmri|init|baseline|screening|v0+)"
                r"|m(?:onth)?[_ ]?(?P<m>\d+)"
                r"|y(?:ear)?[_ ]?(?P<y>\d+)"
                r"|w(?:eek)?[_ ]?(?P<w>\d+)"
                r"|d(?:ay)?[_ ]?(?P<d>\d+))$")        # d0123 = days since entry (OASIS)


    def split_cell(x):
        """'[a, b]' / "['a', 'b']" -> ['a', 'b'];  single value -> [value].
        The master stores a list when a window holds more than one acquisition."""
        if isinstance(x, str):
            s = x.strip()
            if s.startswith("[") and s.endswith("]"):
                try:
                    v = ast.literal_eval(s)
                    v = list(v) if isinstance(v, (list, tuple)) else [v]
                except (ValueError, SyntaxError):
                    v = [p.strip().strip("'\"") for p in s[1:-1].split(",")]
                v = [None if (p is None or str(p).strip().lower() in ("", "nan", "none"))
                     else p for p in v]
                return v if v else [None]
        return [x]


    def explode_lists(L, cols):
        """Explode rows whose cells are lists, keeping value / acq / id aligned."""
        is_list = np.zeros(len(L), dtype=bool)
        for c in cols:
            is_list |= L[c].astype("string").str.strip().str.startswith("[", na=False).to_numpy()
        if not is_list.any():
            return L
        X = L[is_list].copy()
        parts = {c: X[c].map(split_cell) for c in cols}
        n = pd.concat([p.map(len) for p in parts.values()], axis=1).max(axis=1)
        for c in cols:   # align lengths: scalar -> repeated; short list -> padded with None
            X[c] = [(p * k if len(p) == 1 else p + [None] * (k - len(p)))
                    for p, k in zip(parts[c], n)]
        X = X.explode(cols)
        return pd.concat([L[~is_list], X], ignore_index=True)


    def to_time(acq):
        """Convert the _acq column (date OR visit code) to (t_days, kind).
        Dates -> days since 1970 (kind 'date'); visits bl/mXX/yXX/wXX/dXXXX -> days
        since baseline (kind 'visit'). Dates are only compared with dates, visits with visits."""
        s = acq.astype("string").str.strip()
        s = s.str.replace(r"(\.0+)$", "", regex=True)          # 20120110.0 -> 20120110
        t = pd.Series(np.nan, index=s.index, dtype=float)
        kind = pd.Series(None, index=s.index, dtype=object)

        is_date = s.str.contains(DATE_LIKE, regex=True, na=False)
        if is_date.any():
            d = s[is_date]
            ymd8 = d.str.fullmatch(r"(?:19|20)\d{6}")
            dt = pd.Series(pd.NaT, index=d.index, dtype="datetime64[ns]")
            if ymd8.any():
                dt[ymd8] = pd.to_datetime(d[ymd8], format="%Y%m%d", errors="coerce")
            if (~ymd8).any():
                try:
                    dt[~ymd8] = pd.to_datetime(d[~ymd8], errors="coerce", format="mixed")
                except (TypeError, ValueError):
                    dt[~ymd8] = pd.to_datetime(d[~ymd8], errors="coerce")
            ok = dt.notna()
            idx = dt[ok].index
            t.loc[idx] = (dt[ok] - pd.Timestamp("1970-01-01")).dt.days.astype(float)
            kind.loc[idx] = "date"

        rest = kind.isna() & s.notna()
        if rest.any():
            ex = s[rest].str.lower().str.extract(VISIT_RE)
            days = pd.to_numeric(ex["m"], errors="coerce") * DAYS_PER_MONTH
            days = days.fillna(pd.to_numeric(ex["y"], errors="coerce") * 12 * DAYS_PER_MONTH)
            days = days.fillna(pd.to_numeric(ex["w"], errors="coerce") * 7)
            days = days.fillna(pd.to_numeric(ex["d"], errors="coerce"))
            days = days.where(ex["bl"].isna(), 0.0)
            ok = days.notna()
            t.loc[days[ok].index] = days[ok]
            kind.loc[days[ok].index] = "visit"
        return t, kind


    UID2BANK = df.set_index("uid")[COHORT_COL] if COHORT_COL else pd.Series(dtype=object)
    _UNIT_DAYS = {"days": 1.0, "weeks": 7.0, "months": DAYS_PER_MONTH, "years": 365.25}


    def resolve_time(F, id_col=None):
        """t / kind of a long table with columns uid + acq (+ optional id).
        1) date or visit code (to_time)
        2) datasets in ID_DAYCODE_BANKS: days taken from the image ID (..._d0045) -> override
        3) plain number -> unit from NUMERIC_ACQ_UNIT[dataset]"""
        t, kind = to_time(F["acq"])
        bank = F["uid"].map(UID2BANK)

        if id_col is not None and ID_DAYCODE_BANKS:
            sel = bank.isin(ID_DAYCODE_BANKS) & F[id_col].notna()
            if sel.any():
                dd = F.loc[sel, id_col].astype(str).str.extract(r"(?:^|_)d(\d{3,5})(?:_|$)")[0]
                dd = pd.to_numeric(dd, errors="coerce").dropna()
                t.loc[dd.index] = dd.astype(float)
                kind.loc[dd.index] = "visit"
                F.loc[dd.index, "acq"] = ["d%04d" % v for v in dd]   # mostra o dia real

        num = pd.to_numeric(F["acq"], errors="coerce")
        for b, unit in NUMERIC_ACQ_UNIT.items():
            sel = t.isna() & num.notna() & (bank == b)
            if not sel.any() or unit is None:
                continue
            if isinstance(unit, (int, float)):          # visit number x interval (months)
                t.loc[sel] = (num[sel] - 1).clip(lower=0) * unit * DAYS_PER_MONTH
                kind.loc[sel] = f"visitnum_{b}"          # only paired with visit numbers of the same dataset
            else:
                t.loc[sel] = num[sel] * _UNIT_DAYS[unit]
                kind.loc[sel] = "visit"
        return t, kind


    def report_unparsed(name, bad):
        """Report, per dataset, how many acquisitions have no usable date / visit and
        examples of what the _acq column holds (to adjust the parser)."""
        bad = bad.assign(bank=bad["uid"].map(UID2BANK))
        _print(f"  [INFO] {name}: {len(bad)} acquisitions without a usable date / visit (ignored)")
        if not VERBOSE_UNPARSED:
            return
        for bank, g in bad.groupby("bank", dropna=False):
            ex = g["acq"].astype("string").fillna("<empty>").value_counts().head(4)
            ex = ", ".join(f"'{k}'({v})" for k, v in ex.items())
            _print(f"         {str(bank):<10} {len(g):>6}  e.g. {ex}")


    def long_table(name):
        """Long table (uid, stem, window, value, acq, img_id, t, kind) of a measure
        or group. Returns None when no column exists."""
        stems = MEASURE_GROUPS.get(name, [name])
        parts = []
        for stem in stems:
            for w in windows_of(stem):
                v, a, i = f"{stem}_{w}", f"{stem}_{w}_acq", f"{stem}_{w}_id"
                sub = pd.DataFrame({
                    "uid": df["uid"].values,
                    "stem": stem,
                    "window": w,
                    "value": df[v].values,
                    "acq": df[a].values if a in df.columns else np.nan,
                    "img_id": df[i].values if i in df.columns else np.nan,
                    # scanner per acquisition (aligned with acq / img_id)
                    **{ch: (df[f"{stem}_{w}_{ch}"].values if f"{stem}_{w}_{ch}" in df.columns
                            else np.nan) for ch in EXTRA_CHANNELS},
                })
                parts.append(sub[sub["value"].notna() | sub["img_id"].notna()])
        if not parts:
            return None
        L = pd.concat(parts, ignore_index=True)
        L = explode_lists(L, ["value", "acq", "img_id"] + EXTRA_CHANNELS)
        L = L[L["value"].notna() | L["img_id"].notna()]
        num = pd.to_numeric(L["value"], errors="coerce")
        L["value"] = num.where(num.notna(), L["value"])
        L["t"], L["kind"] = resolve_time(L, id_col="img_id")
        bad = L[L["t"].isna()]
        if len(bad):
            report_unparsed(name, bad)
        # alert: dataset where ALL acquisitions fall on the same time (e.g. all 0)
        g = L[L["t"].notna()].assign(bank=lambda x: x["uid"].map(UID2BANK)).groupby("bank")
        for b, sub in g:
            if len(sub) >= 20 and sub["t"].nunique() == 1 and sub["kind"].iloc[0] != "date":
                _print(f"  [WARN] {name} / {b}: all {len(sub)} acquisitions have the same time "
                      f"({sub['acq'].iloc[0]!r}) -> the date was probably lost in the expansion")
        return L[L["t"].notna()].reset_index(drop=True)


    # =============================================================================
    # 4. ANCHOR
    # =============================================================================
    A = long_table(ANCHOR)
    if A is None or A.empty:
        raise ValueError(f"No acquisition found for the anchor {ANCHOR}.")

    A = A.sort_values(["uid", "t"]).reset_index(drop=True)
    A["event"] = np.arange(len(A))
    A["anchor_event_n"] = A.groupby("uid").cumcount() + 1
    _print(f"\nAnchor events {ANCHOR}: {len(A)} in {A['uid'].nunique()} subjects")
    if COHORT_COL:
        _b = A.assign(bank=A["uid"].map(UID2BANK)).groupby("bank") \
              .agg(events=("uid", "size"), subjects=("uid", "nunique"))
        _print(_b.to_string())
    _unt = f"{ANCHOR}_untimed_n"
    if _unt in df.columns:
        only_unt = (pd.to_numeric(df[_unt], errors="coerce").fillna(0) > 0) & \
                   ~df["uid"].isin(A["uid"])
        if only_unt.any():
            _print(f"  [INFO] {int(only_unt.sum())} subjects have {ANCHOR} only 'untimed' (no date) "
                  f"-> cannot be paired: {df.loc[only_unt, COHORT_COL].value_counts().to_dict()}")

    out = A.rename(columns={"stem": f"{ANCHOR}_stem", "window": f"{ANCHOR}_window",
                            "value": f"{ANCHOR}_value", "acq": f"{ANCHOR}_acq",
                            "img_id": f"{ANCHOR}_id", "t": "anchor_t",
                            **{ch: f"{ANCHOR}_{ch}" for ch in EXTRA_CHANNELS},
                            "kind": "anchor_kind"}).set_index("event")
    if not MEASURE_GROUPS.get(ANCHOR):
        out = out.drop(columns=[f"{ANCHOR}_stem"])
    if out[f"{ANCHOR}_id"].isna().all():
        out = out.drop(columns=[f"{ANCHOR}_id"])
    for ch in EXTRA_CHANNELS:
        col = f"{ANCHOR}_{ch}"
        if col in out.columns:
            out[col] = out[col].replace({"": np.nan, "NA": np.nan})
            if out[col].isna().all():
                out = out.drop(columns=[col])

    # =============================================================================
    # 5. PAIRING
    # =============================================================================
    MARGIN_DAYS = MARGIN_MONTHS * DAYS_PER_MONTH + (0.5 if MARGIN_MONTHS == 0 else 0)


    def match(name):
        L = long_table(name)
        if L is None or L.empty:
            _print(f"  [WARN] {name}: no column / acquisition found (ignored)")
            return None
        m = A[["event", "uid", "kind", "t"]].rename(columns={"t": "anchor_t"}) \
            .merge(L, on=["uid", "kind"], how="inner")
        m["delta_days"] = m["t"] - m["anchor_t"]
        mm = MARGIN_MONTHS_BY_FEATURE.get(name, MARGIN_MONTHS)
        m = m[m["delta_days"].abs() <= mm * DAYS_PER_MONTH + (0.5 if mm == 0 else 0)].copy()
        m["_abs"] = m["delta_days"].abs()
        m = m.sort_values(["event", "_abs", "t"]).drop_duplicates("event")

        if not ALLOW_IMAGE_REUSE and m["img_id"].notna().any():
            # each image stays only with its closest anchor event
            has = m["img_id"].notna()
            best = m[has].sort_values(["_abs", "event"]).drop_duplicates(["uid", "img_id"])
            m = pd.concat([best, m[~has]])

        keep = ["value", "acq", "img_id", "window", "delta_days"] + \
               [ch for ch in EXTRA_CHANNELS if ch in m.columns]
        if name in MEASURE_GROUPS:
            keep = ["stem"] + keep
        res = m.set_index("event")[keep]
        if res["img_id"].isna().all():
            res = res.drop(columns=["img_id"])
        for ch in EXTRA_CHANNELS:
            if ch in res.columns:
                res[ch] = res[ch].replace({"": np.nan, "NA": np.nan})
                if res[ch].isna().all():
                    res = res.drop(columns=[ch])
        res = res.rename(columns={"img_id": "id"}).add_prefix(f"{name}_")
        res[f"{name}_delta_days"] = res[f"{name}_delta_days"].round().astype("Int64")
        return res


    features = [f for f in dict.fromkeys(REQUIRED + OPTIONAL) if f != ANCHOR]
    opt_only = [f for f in dict.fromkeys(OPTIONAL) if f not in REQUIRED and f != ANCHOR]
    req_ok = [f for f in dict.fromkeys(REQUIRED) if f != ANCHOR]
    found = []

    _print("\nPairing:")
    for f in features:
        res = match(f)
        if res is None:
            if f in req_ok:
                raise ValueError(f"Required measure '{f}' does not exist in the expanded master.")
            continue
        out = out.join(res)
        found.append(f)
        _print(f"  {f:<16} {'REQUIRED' if f in req_ok else 'optional'}  "
              f"-> {len(res):>6} events with a match within "
              f"±{MARGIN_MONTHS_BY_FEATURE.get(f, MARGIN_MONTHS)} m")

    # filter on the required measures
    mask = pd.Series(True, index=out.index)
    _print("\nRequired-measure filter:")
    for f in req_ok:
        ok = out[f"{f}_delta_days"].notna()
        mask &= ok
        _print(f"  {f:<16} {int(ok.sum()):>6} pass  |  cumulative: {int(mask.sum()):>6} events, "
              f"{out.loc[mask, 'uid'].nunique():>5} subjects")
    out = out[mask].copy()

    opt_found = [f for f in opt_only if f in found]
    out["n_optional_matched"] = sum(out[f"{f}_delta_days"].notna().astype(int) for f in opt_found) \
        if opt_found else 0

    # selection of events per subject
    if ANCHOR_EVENTS == "first":
        out = out.sort_values(["uid", "anchor_t"]).drop_duplicates("uid", keep="first")
    elif ANCHOR_EVENTS == "last":
        out = out.sort_values(["uid", "anchor_t"]).drop_duplicates("uid", keep="last")
    elif ANCHOR_EVENTS == "most_optional":
        out = out.sort_values(["uid", "n_optional_matched", "anchor_t"],
                              ascending=[True, False, True]).drop_duplicates("uid")
    elif ANCHOR_EVENTS != "all":
        raise ValueError("anchor_events must be 'all', 'first', 'last' or 'most_optional'.")
    if ANCHOR_EVENTS != "all":
        _print(f"\nANCHOR_EVENTS='{ANCHOR_EVENTS}': kept 1 event per subject -> {len(out)} events")

    # =============================================================================
    # 6. DIAGNOSIS AT THE ANCHOR + DEMOGRAPHICS
    # =============================================================================
    anchor_win_col = f"{ANCHOR}_window"


    def win_order(w):
        """Chronological order of a window: pre (-inf) < preK (-K) < ltK (K) < gt (+inf)."""
        w = str(w)
        if w == "pre" or w.startswith("pregt"):
            return -1e9
        if w.startswith("pre"):
            return -int(w[3:])
        if w.startswith("lt"):
            return int(w[2:])
        if w.startswith("gt"):
            return 1e9
        return np.nan


    diag_cols = [c for c in df.columns if re.fullmatch("diag_" + WIN_RE, c)]
    out["diag_at_anchor"] = pd.NA
    out["diag_source"] = pd.NA
    if diag_cols:
        D = df[["uid"] + diag_cols].melt(id_vars="uid", var_name="col", value_name="diag")
        D = D[D["diag"].notna() & (D["diag"].astype(str).str.strip() != "")].copy()
        D["w"] = D["col"].str[5:]
        D["w_order"] = D["w"].map(win_order)
        acq_cols = [f"{c}_acq" for c in diag_cols if f"{c}_acq" in df.columns]
        if acq_cols:
            Dacq = df[["uid"] + acq_cols].melt(id_vars="uid", var_name="col", value_name="acq")
            Dacq["col"] = Dacq["col"].str[:-4]
            D = D.merge(Dacq, on=["uid", "col"], how="left")
        else:
            D["acq"] = np.nan
        # list cells ("['CN', 'MCI']") = several diagnoses in one window -> one per row
        D = D.reset_index(drop=True)
        D["_row"] = D.index
        D = explode_lists(D, ["diag", "acq"])
        D["pos"] = D.groupby("_row").cumcount()
        D["diag"] = D["diag"].astype("string").str.strip()
        D = D[D["diag"].notna() & ~D["diag"].str.lower().isin(["", "nan", "none"])]
        D = D.reset_index(drop=True)
        D["t"], D["kind"] = resolve_time(D)

        # (1) diagnosis in the anchor's window (if several: the closest to the
        #     anchor in time; without a comparable date: the last of the list)
        same = D[["uid", "w", "diag", "t", "kind", "pos"]] \
            .rename(columns={"w": anchor_win_col, "t": "_dt", "kind": "_dk"})
        tmp = out[["uid", anchor_win_col, "anchor_t", "anchor_kind"]].reset_index() \
            .merge(same, on=["uid", anchor_win_col], how="inner")
        tmp["_abs"] = np.where(tmp["_dk"] == tmp["anchor_kind"],
                               (tmp["_dt"] - tmp["anchor_t"]).abs(), 1e12)
        tmp = tmp.sort_values(["event", "_abs", "pos"], ascending=[True, True, False]) \
                 .drop_duplicates("event").set_index("event")
        out.loc[tmp.index, "diag_at_anchor"] = tmp["diag"]
        out.loc[out["diag_at_anchor"].notna(), "diag_source"] = "same_window"

        # (2) last diagnosis recorded up to the anchor (by date; otherwise by window order)
        need = out["diag_at_anchor"].isna()
        if need.any():
            E = out.loc[need, ["uid", "anchor_t", "anchor_kind", anchor_win_col]].reset_index()
            E["a_order"] = E[anchor_win_col].map(win_order)
            M = E.merge(D[["uid", "diag", "t", "kind", "w_order"]], on="uid", how="inner")
            by_time = (M["kind"] == M["anchor_kind"]) & (M["t"] <= M["anchor_t"])
            by_win = M["kind"].ne(M["anchor_kind"]) & (M["w_order"] <= M["a_order"])
            M = M[by_time | by_win].copy()
            # prefer diagnoses with a comparable date (same kind); among them, the most recent
            M["_same"] = (M["kind"] == M["anchor_kind"]).astype(int)
            M["_key"] = np.where(M["_same"] == 1, M["t"], M["w_order"])
            last = M.sort_values(["event", "_same", "_key"]).drop_duplicates("event", keep="last") \
                    .set_index("event")["diag"]
            out.loc[last.index, "diag_at_anchor"] = last
            out.loc[last.index, "diag_source"] = "last_before_anchor"

    # (3) no earlier diagnosis -> Lst_diag (last diagnosis of the subject)
    need = out["diag_at_anchor"].isna()
    if need.any() and "Lst_diag" in df.columns:
        lst = out.loc[need, "uid"].map(df.set_index("uid")["Lst_diag"])
        lst = lst.map(lambda x: split_cell(x)[-1] if pd.notna(x) else x)   # list -> last
        out.loc[need, "diag_at_anchor"] = lst
        out.loc[need & lst.notna(), "diag_source"] = "Lst_diag"

    # =============================================================================
    # 6b. months_from_1st_diag and age_at_anchor
    # =============================================================================
    # months_from_1st_diag = (anchor date - 1st_diag_date) in months
    #   'exact'  -> anchor and 1st_diag_date on the same time axis (date-date or visit-visit)
    #   'window' -> estimated from the centre of the anchor window (ltK -> K-3 ; preK -> -(K-3))
    # age_at_anchor = age_at_1st_diag + months_from_1st_diag/12
    #   (without age_at_1st_diag -> age_at_baseline + months_from_1st_diag/12, assuming
    #    baseline = first diagnostic assessment)


    def _first(x):
        return split_cell(x)[0] if pd.notna(x) else x


    def window_center(w):
        w = str(w)
        if w.startswith("lt"):
            return int(w[2:]) - 3.0
        if w.startswith("pre") and w[3:].isdigit():
            return -(int(w[3:]) - 3.0)
        return np.nan              # 'pre', 'pregtK', 'gt240' -> no defined centre


    u2 = df.set_index("uid")
    ev = out[["uid", "anchor_t", "anchor_kind", anchor_win_col]].copy()
    out["months_from_1st_diag"] = np.nan
    out["months_from_1st_diag_source"] = pd.NA

    if "1st_diag_date" in df.columns:
        F = pd.DataFrame({"uid": ev["uid"].values,
                          "acq": ev["uid"].map(u2["1st_diag_date"]).map(_first).values},
                         index=ev.index)
        F["t"], F["kind"] = resolve_time(F)
        same = F["t"].notna() & (F["kind"] == ev["anchor_kind"])
        out.loc[same, "months_from_1st_diag"] = \
            (ev.loc[same, "anchor_t"] - F.loc[same, "t"]) / DAYS_PER_MONTH
        out.loc[same, "months_from_1st_diag_source"] = "exact"

    need = out["months_from_1st_diag"].isna()
    wc = ev.loc[need, anchor_win_col].map(window_center)
    out.loc[wc.dropna().index, "months_from_1st_diag"] = wc.dropna()
    out.loc[wc.dropna().index, "months_from_1st_diag_source"] = "window"
    out["months_from_1st_diag"] = out["months_from_1st_diag"].astype(float).round(1)

    yrs = out["months_from_1st_diag"] / 12
    out["age_at_anchor"] = np.nan
    out["age_at_anchor_source"] = pd.NA
    for col in ["age_at_1st_diag", "age_at_baseline"]:
        if col not in df.columns:
            continue
        base_age = pd.to_numeric(out["uid"].map(u2[col]).map(_first), errors="coerce")
        fill = out["age_at_anchor"].isna() & base_age.notna() & yrs.notna()
        out.loc[fill, "age_at_anchor"] = base_age[fill] + yrs[fill]
        out.loc[fill, "age_at_anchor_source"] = \
            col + "+" + out.loc[fill, "months_from_1st_diag_source"].astype(str)
    out["age_at_anchor"] = out["age_at_anchor"].astype(float).round(1)

    _print("\nmonths_from_1st_diag:", out["months_from_1st_diag_source"]
          .value_counts(dropna=False).to_dict())
    _print("age_at_anchor:       ", out["age_at_anchor_source"]
          .value_counts(dropna=False).to_dict())

    keep_extra = [c for c in KEEP_COLUMNS if c in df.columns]
    out = out.reset_index().merge(df[["uid"] + keep_extra], on="uid", how="left")

    front = ["uid"] + keep_extra + ["anchor_event_n", "anchor_kind"]
    anchor_cols = [c for c in out.columns if c.startswith(f"{ANCHOR}_")]
    feat_cols = [c for f in found for c in out.columns if c.startswith(f"{f}_")]
    rest = ["age_at_anchor", "age_at_anchor_source",
            "months_from_1st_diag", "months_from_1st_diag_source",
            "diag_at_anchor", "diag_source", "n_optional_matched"]
    order = list(dict.fromkeys(front + anchor_cols + [c for c in rest if c in out.columns] + feat_cols))
    out = out[order].sort_values(["uid", "anchor_event_n"]).reset_index(drop=True)

    # laboratory / kit of every matched fluid biomarker (which platform the value is from)
    for _f in [ANCHOR] + found:
        _stem = out[f"{_f}_stem"] if f"{_f}_stem" in out.columns else pd.Series(_f, index=out.index)
        for _k in FLUID_META:
            _vals = pd.Series(np.nan, index=out.index, dtype=object)
            for _s in _stem.dropna().unique():
                if f"{_s}_{_k}" in df.columns:
                    _sel = (_stem == _s) & out[f"{_f}_acq"].notna()
                    _vals[_sel] = out.loc[_sel, "uid"].map(u2[f"{_s}_{_k}"]).values
            if _vals.notna().any():
                out[f"{_f}_{_k}"] = _vals

    # =============================================================================
    # 7. LIST OF loniIDs
    # =============================================================================
    img_rows = []
    for f in [ANCHOR] + found:
        idc = f"{f}_id"
        if idc not in out.columns:
            continue
        sub = out[out[idc].notna()]
        img_rows.append(pd.DataFrame({
            "uid": sub["uid"],
            COHORT_COL: sub[COHORT_COL] if COHORT_COL in sub.columns else np.nan,
            "ID": sub["ID"] if "ID" in sub.columns else np.nan,
            "anchor_event_n": sub["anchor_event_n"],
            "anchor_acq": sub[f"{ANCHOR}_acq"],
            "feature": f,
            "stem": sub[f"{f}_stem"] if f"{f}_stem" in sub.columns else f,
            "loniID": sub[idc],
            "acq": sub[f"{f}_acq"],
            "delta_days": sub[f"{f}_delta_days"] if f != ANCHOR else 0,
            "required": f in req_ok or f == ANCHOR,
        }))
    imgs = pd.concat(img_rows, ignore_index=True) if img_rows else pd.DataFrame()

    # =============================================================================
    # 7b. SUMMARY SHEET: basic data + (loniID or value) + date per measure
    # =============================================================================
    # Image (has _id)     -> <measure> = loniID
    # NPS / fluid marker  -> <measure> = value
    # Always              -> <measure>_date = date / visit of the acquisition
    summary = out[["uid"] + keep_extra + ["anchor_event_n",
                                          "age_at_anchor", "months_from_1st_diag",
                                          "months_from_1st_diag_source",
                                          "diag_at_anchor", "diag_source"]].copy()
    for f in [ANCHOR] + found:
        idc, valc = f"{f}_id", f"{f}_value"
        if f"{f}_stem" in out.columns:                       # group: which tracer / stem
            summary[f"{f}_stem"] = out[f"{f}_stem"]
        if idc in out.columns and out[idc].notna().any():
            summary[f"{f}_loniID"] = out[idc]
        else:
            summary[f] = out[valc]
        summary[f"{f}_date"] = out[f"{f}_acq"]
        for ch in EXTRA_CHANNELS + FLUID_META:
            if f"{f}_{ch}" in out.columns:
                summary[f"{f}_{ch}"] = out[f"{f}_{ch}"]
        if SUMMARY_INCLUDE_DELTA and f != ANCHOR:
            summary[f"{f}_delta_days"] = out[f"{f}_delta_days"]

    # =============================================================================
    # 8. SAVE + SUMMARY
    # =============================================================================
    pairs_path = os.path.join(out_dir, f"{OUTPUT_PREFIX}_pairs.csv")
    out.to_csv(pairs_path, index=False)

    _print("\n" + "=" * 70)
    _print(f"PAIRED EVENTS: {len(out)}  |  SUBJECTS: {out['uid'].nunique()}")
    _print("=" * 70)
    _print(f"Saved -> {pairs_path}")

    if not imgs.empty:
        imgs_path = os.path.join(out_dir, f"{OUTPUT_PREFIX}_loniIDs.csv")
        txt_path = os.path.join(out_dir, f"{OUTPUT_PREFIX}_loniIDs.txt")
        imgs.to_csv(imgs_path, index=False)
        unique_ids = imgs["loniID"].astype(str).str.strip().drop_duplicates().tolist()
        with open(txt_path, "w") as fh:
            fh.write(",".join(unique_ids))
        _print(f"Saved -> {imgs_path}")
        _print(f"Saved -> {txt_path}  ({len(unique_ids)} unique loniIDs)")
        _print("\nImages per measure:")
        _print(imgs.groupby("feature")["loniID"].nunique().to_string())

    summary_path = os.path.join(out_dir, f"{OUTPUT_PREFIX}_summary.csv")
    summary.to_csv(summary_path, index=False)
    _print(f"Saved -> {summary_path}  (summary sheet: {summary.shape[1]} columns)")

    if SAVE_XLSX:
        xlsx_path = os.path.join(out_dir, f"{OUTPUT_PREFIX}.xlsx")
        try:
            with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
                summary.to_excel(xw, sheet_name="summary", index=False)
                out.to_excel(xw, sheet_name="pairs", index=False)
                if not imgs.empty:
                    imgs.to_excel(xw, sheet_name="loniIDs", index=False)
            _print(f"Saved -> {xlsx_path}  (sheets: summary | pairs | loniIDs)")
        except ImportError:
            _print("  [WARN] openpyxl is not installed: pip install openpyxl  (xlsx not saved)")

    if COHORT_COL in out.columns:
        _print("\nSubjects per dataset:")
        _print(out.drop_duplicates("uid")[COHORT_COL].value_counts().to_string())
    if "diag_at_anchor" in out.columns:
        _print("\nDiagnosis at the anchor (events):")
        _print(out["diag_at_anchor"].map(lambda x: str(x).split("+")[0])
              .value_counts(dropna=False).to_string())
        _print("\nSource of the diagnosis:")
        _print(out["diag_source"].value_counts(dropna=False).to_string())
    if opt_found:
        _print("\nCoverage of the optional measures (events):")
        for f in opt_found:
            n = int(out[f"{f}_delta_days"].notna().sum())
            _print(f"  {f:<16} {n:>6} / {len(out)}  ({100 * n / max(len(out), 1):.1f}%)")
    _print("\nDone.")
    files = {"pairs": pairs_path, "summary": summary_path}
    if not imgs.empty:
        files.update(images=imgs_path, image_ids=txt_path)
    if SAVE_XLSX:
        files["xlsx"] = xlsx_path
    write_provenance(pairs_path, step="anchor_cohort", config=cfg, inputs=[INPUT_PATH], table=out,
                     parameters={"anchor": ANCHOR, "required": REQUIRED, "optional": OPTIONAL,
                                 "margin_months": MARGIN_MONTHS,
                                 "margin_months_by_feature": MARGIN_MONTHS_BY_FEATURE,
                                 "anchor_events": ANCHOR_EVENTS, "select_databanks": SELECT_DATABANKS,
                                 "cohort_filter_path": COHORT_FILTER_PATH})
    return CohortResult(table=out, images=imgs, summary=summary, files=files)
