"""Shared helpers of the cohort tools (they read the expanded master).

* :func:`vrs_measure_stems` - name or group (``PET_ABETA``, ``PLASMA_PTAU217``) -> stems
* :func:`vrs_read_expanded` - reads ONLY the columns of the requested measures
* :func:`vrs_events` - long table, one row per acquisition: uid, label, stem,
  window, value, acq, img_id, image_ids, manufacturer, model, suvr, centiloid,
  t (days), kind
* :func:`vrs_time` - acquisition token -> days on ONE axis per subject

Time axis: ``_acq`` holds a calendar date (ADNI, AIBL, NIFD, 4RTNI, ABVIB,
MCSA, NACC, PREVENT-AD) or a number of days since entry (OASIS, A4). Both
become ``t`` in days; ``kind`` says which axis (``date`` or ``visit``).
Acquisitions are only compared inside the same subject AND the same kind, so a
date is never subtracted from a day offset.

The module-level settings below are the defaults of ``configs/``;
:func:`configure` replaces them with the values of a loaded configuration.
"""
from __future__ import annotations

import ast
import numpy as np
import pandas as pd
import re

import logging as _logging
from ..logging_utils import make_print as _make_print

log = _logging.getLogger(__name__)
_print = _make_print(log)

VRS_DAYS_PER_MONTH = 30.4375
VRS_WIN_RX = r"(pre\d*|pregt\d+|lt\d+|gt\d+)"
VRS_CHANNELS = ("id", "image_ids", "manufacturer", "model", "suvr", "centiloid")

# pure numbers in '_acq': unit per databank (same table as the Cohort creator)
VRS_NUMERIC_ACQ_UNIT = {"OASIS": "days", "A4": "days", "NIFD": "months",
                        "4RTNI": "months", "MCSA": None, "ABVIB": None}
_VRS_UNIT_DAYS = {"days": 1.0, "weeks": 7.0, "months": VRS_DAYS_PER_MONTH, "years": 365.25}

# OR-groups: any member satisfies the measure; '<label>_stem' says which one
VRS_MEASURE_GROUPS = {
    "PET_ABETA": ["ABETA_AV45", "ABETA_FBB", "ABETA_Flutemetamol",
                  "ABETA_NAV4694", "ABETA_PIB"],
    "PET_TAU":   ["TAU_Flortaucipir", "TAU_MK6240", "TAU_PI2620"],
}
VRS_FLUID_GROUP_PATTERNS = {
    "PLASMA_PTAU217": r"^Plasma_pTau217_(?!ratio|AB42)",
    "PLASMA_PTAU181": r"^Plasma_pTau181_",
    "PLASMA_NFL":     r"^Plasma_(NfL|NFL)(_|$)",
    "PLASMA_GFAP":    r"^Plasma_GFAP(_|$)",
    "PLASMA_AB42_40": r"^Plasma_(AB42_40|AB_ratio)_",
    "CSF_PTAU":       r"^(PTAU$|CSF_pTau181_)",
    "CSF_AB42":       r"^(ABETA42(_lumipulse)?$|CSF_AB42_(?!40))",
}
_VRS_EMPTY = ("", "nan", "none", "nat", "na")


def vrs_header(path: str) -> list:
    return pd.read_csv(path, nrows=0).columns.tolist()


def vrs_read_csv_lean(path: str, usecols, text_cols=None, chunksize: int = 4000) -> pd.DataFrame:
    """
    read_csv for VERY wide files. pd.read_csv(..., low_memory=False) tokenises
    the whole file (every column, every row) before applying usecols — on the
    expanded master that costs several GB even to read 50 columns. Reading in
    row chunks keeps the peak proportional to the columns actually used.
    text_cols (dict col -> str) are read as text; the others are converted to
    numbers afterwards when every non-empty value is numeric (same result as
    the whole-file type inference).
    """
    text_cols = dict(text_cols or {})
    use = list(usecols)
    df = pd.concat(pd.read_csv(path, usecols=use, dtype=str, chunksize=chunksize),
                   ignore_index=True)
    for c in df.columns:
        if c in text_cols:
            continue
        num = pd.to_numeric(df[c], errors="coerce")
        if num.notna().sum() == df[c].notna().sum():
            df[c] = num
    return df[[c for c in use if c in df.columns]]


def vrs_all_stems(header) -> list:
    """Every measure stem of an expanded master (columns <stem>_<window>_acq)."""
    rx = re.compile(r"^(.+)_" + VRS_WIN_RX + r"_acq$")
    return sorted({m.group(1) for c in header if (m := rx.match(c))})


def vrs_measure_stems(name, header) -> list:
    """Stems behind a measure name. Accepts a stem ('MRI_T1w'), a group
    ('PET_ABETA', 'PLASMA_PTAU217' …) or a list of either. EXACT stem match:
    'MRI_T2w' never includes 'MRI_T2w_FLAIR', 'CDR' never includes 'CDR_SB'."""
    stems = vrs_all_stems(header)
    out = []
    for n in ([name] if isinstance(name, str) else list(name)):
        if n in VRS_MEASURE_GROUPS:
            out += [s for s in VRS_MEASURE_GROUPS[n] if s in stems]
        elif n in VRS_FLUID_GROUP_PATTERNS:
            out += [s for s in stems if re.match(VRS_FLUID_GROUP_PATTERNS[n], s)]
        elif n in stems:
            out.append(n)
    return list(dict.fromkeys(out))


def vrs_read_expanded(path: str, measures, keep_cols=(), databanks="ALL") -> pd.DataFrame:
    """
    Lean read of the expanded master: identification + keep_cols + the window
    columns (<stem>_<w>, _acq, _id, _image_ids, scanner and PET value channels)
    of the requested measures only. IDs / acquisition tokens are read as text.
    """
    header = vrs_header(path)
    stems = set()
    for m in measures:
        stems.update(vrs_measure_stems(m, header))
    tail = re.compile(r"_" + VRS_WIN_RX + r"(_acq|_" + "|_".join(VRS_CHANNELS) + r")?$")
    base = {"uid", "ID", "Databank", "dataset"} | set(keep_cols)
    use = [c for c in header
           if c in base or ((m := tail.search(c)) and c[:m.start()] in stems)]
    as_text = {c: str for c in use if c in ("uid", "ID") or tail.search(c)}
    df = vrs_read_csv_lean(path, use, as_text)
    bank_col = "Databank" if "Databank" in df.columns else "dataset"
    if "uid" not in df.columns:
        df["uid"] = df[bank_col].astype(str) + "_" + df["ID"].astype(str)
    if databanks != "ALL":
        df = df[df[bank_col].isin(databanks)]
    df = df.drop_duplicates("uid").reset_index(drop=True)
    _print(f"  read {len(use)} of {len(header)} columns · {len(df):,} subjects "
          f"· stems: {sorted(stems)}")
    return df


def vrs_split_cell(x) -> list:
    """'[a, b]' / "['a', 'b']" -> ['a', 'b'] ; scalar -> [scalar] ; empty -> [None]."""
    if isinstance(x, str):
        s = x.strip()
        if s.startswith("[") and s.endswith("]"):
            try:
                v = ast.literal_eval(s)
                v = list(v) if isinstance(v, (list, tuple)) else [v]
            except (ValueError, SyntaxError):
                v = [p.strip().strip("'\"") for p in s[1:-1].split(",")]
            return [None if (p is None or str(p).strip().lower() in _VRS_EMPTY) else p
                    for p in v] or [None]
        return [None if s.lower() in _VRS_EMPTY else x]
    return [None if (x is None or (isinstance(x, float) and np.isnan(x))) else x]


def vrs_time(acq: pd.Series, bank: pd.Series, units: dict = None):
    """
    Acquisition tokens -> (t in days, kind).
      'YYYY-MM-DD' (or any parseable date) -> days since 1970-01-01, kind 'date'
      pure number -> days via units[bank] ('days' | 'weeks' | 'months' | 'years'),
                     kind 'visit'. Banks with unit None are left without time.
    """
    units = VRS_NUMERIC_ACQ_UNIT if units is None else units
    s = acq.astype("string").str.strip().str.replace(r"\.0+$", "", regex=True)
    t = pd.Series(np.nan, index=s.index, dtype=float)
    kind = pd.Series(None, index=s.index, dtype=object)
    is_date = s.str.contains(r"\d{4}-\d{1,2}-\d{1,2}|\d{1,2}/\d{1,2}/\d{2,4}", regex=True, na=False)
    if is_date.any():
        dt = pd.to_datetime(s[is_date].str.slice(0, 10), errors="coerce", format="%Y-%m-%d")
        if dt.isna().any():
            dt = dt.fillna(pd.to_datetime(s[is_date][dt.isna()], errors="coerce"))
        ok = dt.notna()
        t.loc[dt[ok].index] = (dt[ok] - pd.Timestamp("1970-01-01")).dt.days.astype(float)
        kind.loc[dt[ok].index] = "date"
    num = pd.to_numeric(s.where(~is_date), errors="coerce")
    for b, unit in units.items():
        sel = t.isna() & num.notna() & (bank == b)
        if sel.any() and unit in _VRS_UNIT_DAYS:
            t.loc[sel] = num[sel] * _VRS_UNIT_DAYS[unit]
            kind.loc[sel] = "visit"
    return t, kind


def vrs_events(df: pd.DataFrame, measure, label: str = None, units: dict = None) -> pd.DataFrame:
    """
    One row per acquisition of a measure (stem, group or list of stems):
      uid, label, stem, window, value, acq, img_id, image_ids, manufacturer,
      model, suvr, centiloid, t (days), kind ('date' | 'visit')
    List cells of the expanded master (several acquisitions in one window) are
    exploded keeping value / acq / id / channels aligned by position.
    """
    stems = vrs_measure_stems(measure, df.columns)
    label = label or (measure if isinstance(measure, str) else "+".join(measure))
    chans = [c for c in VRS_CHANNELS if c != "id"]
    parts = []
    for stem in stems:
        rx = re.compile(re.escape(stem) + "_" + VRS_WIN_RX + "$")
        for c in df.columns:
            m = rx.fullmatch(c)
            if not m:
                continue
            w = m.group(1)
            acq = df.get(f"{c}_acq")
            if acq is None:
                continue
            sel = acq.notna()
            if not sel.any():
                continue
            sub = pd.DataFrame({"uid": df.loc[sel, "uid"].values, "stem": stem, "window": w,
                                "value": df.loc[sel, c].values, "acq": acq[sel].values,
                                "img_id": (df.loc[sel, f"{c}_id"].values
                                           if f"{c}_id" in df.columns else None)})
            for ch in chans:
                sub[ch] = df.loc[sel, f"{c}_{ch}"].values if f"{c}_{ch}" in df.columns else None
            parts.append(sub)
    cols = ["uid", "label", "stem", "window", "value", "acq", "img_id"] + chans + ["t", "kind"]
    if not parts:
        return pd.DataFrame(columns=cols)
    L = pd.concat(parts, ignore_index=True)
    lst_cols = ["value", "acq", "img_id"] + chans
    is_list = np.zeros(len(L), dtype=bool)
    for c in lst_cols:
        is_list |= L[c].astype("string").str.strip().str.startswith("[", na=False).to_numpy()
    if is_list.any():
        X = L[is_list].copy()
        split = {c: X[c].map(vrs_split_cell) for c in lst_cols}
        n = pd.concat([p.map(len) for p in split.values()], axis=1).max(axis=1)
        for c in lst_cols:                 # scalar -> repeated ; short list -> padded
            X[c] = [(p * k if len(p) == 1 else p + [None] * (k - len(p)))
                    for p, k in zip(split[c], n)]
        L = pd.concat([L[~is_list], X.explode(lst_cols)], ignore_index=True)
    for c in ["img_id"] + chans:           # 'NA' / '' placeholders -> missing
        s = L[c].astype("string").str.strip()
        L[c] = s.where(~s.str.lower().isin(_VRS_EMPTY)).astype(object)
        L.loc[L[c].isna(), c] = None
    L["img_id"] = L["img_id"].map(lambda x: re.sub(r"\.0$", "", x) if isinstance(x, str) else x)
    bank = L["uid"].map(df.drop_duplicates("uid").set_index("uid")[
        "Databank" if "Databank" in df.columns else "dataset"])
    L["t"], L["kind"] = vrs_time(L["acq"], bank, units)
    n_bad = int(L["t"].isna().sum())
    if n_bad:
        ex = L.loc[L["t"].isna()].assign(bank=bank).groupby("bank")["acq"] \
              .agg(lambda x: f"{len(x)} (e.g. {x.iloc[0]!r})").to_dict()
        _print(f"    [INFO] {label}: {n_bad} acquisitions without a usable time (ignored): {ex}")
    L = L[L["t"].notna()].copy()
    L["label"] = label
    L = L.drop_duplicates(["uid", "stem", "t", "img_id"])
    return L[cols].sort_values(["uid", "t", "stem"]).reset_index(drop=True)


def vrs_clean_id_list(ids) -> list:
    """Unique, non-empty image IDs keeping order."""
    out = []
    for i in ids:
        if i is None or (isinstance(i, float) and np.isnan(i)):
            continue
        for part in str(i).split("/"):
            part = part.strip()
            if part and part.lower() not in _VRS_EMPTY and part not in out:
                out.append(part)
    return out



def configure(cfg) -> None:
    """Load the time-axis units and measure groups from a project configuration."""
    global VRS_DAYS_PER_MONTH, VRS_NUMERIC_ACQ_UNIT, VRS_MEASURE_GROUPS, VRS_FLUID_GROUP_PATTERNS
    ta = cfg.integration.get("time_axis", {})
    VRS_DAYS_PER_MONTH = float(ta.get("days_per_month", VRS_DAYS_PER_MONTH))
    VRS_NUMERIC_ACQ_UNIT = dict(ta.get("numeric_acquisition_unit", VRS_NUMERIC_ACQ_UNIT))
    _VRS_UNIT_DAYS["months"] = VRS_DAYS_PER_MONTH
    VRS_MEASURE_GROUPS = {k: list(v) for k, v in cfg.imaging.get("measure_groups", {}).items()}
    VRS_FLUID_GROUP_PATTERNS = dict(cfg.imaging.get("fluid_group_patterns", {}))
