"""Shared clinical helpers: body-mass index, UDS health-history codes,
case-insensitive column lookup, and free-text seizure / antiepileptic flags.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import re

import logging as _logging
from .logging_utils import make_print as _make_print

log = _logging.getLogger(__name__)
_print = _make_print(log)

# =============================================================================
# BMI + UDS HEALTH-HISTORY HELPERS  (used by ADNI, A4, NACC, OASIS-3, PREVENT-AD)
# =============================================================================
_W_SENTINELS = {888, 888.8, 999, 999.9, 8888, 9999}
_H_SENTINELS = {88.8, 99.9, 888, 888.8, 999, 9999}


def _unit_series(unit, index, mapping):
    """unit: 'kg'/'lb'/'cm'/'in' string, or a Series of codes mapped via mapping."""
    if isinstance(unit, pd.Series):
        return pd.to_numeric(unit, errors="coerce").map(mapping)
    return pd.Series(unit, index=index)


def compute_bmi(weight, height, weight_unit="kg", height_unit="cm",
                weight_codes=None, height_codes=None):
    """
    BMI = kg / m². Units may be fixed strings or per-row code Series
    (ADNI/A4 VSWTUNIT 1=lb 2=kg, VSHTUNIT 1=in 2=cm). Rows with a missing unit
    code get it from the magnitude (height > 100 -> cm, else in; weight follows
    the height system). Sentinels (888, 9999, 88.8 …) and negatives -> NaN.
    Implausible BMI (< 10 or > 80) -> NaN.
    """
    w = pd.to_numeric(pd.Series(weight), errors="coerce")
    h = pd.to_numeric(pd.Series(height), errors="coerce").reindex(w.index)
    w = w.where(~w.isin(_W_SENTINELS) & (w > 0))
    h = h.where(~h.isin(_H_SENTINELS) & (h > 0))
    hu = _unit_series(height_unit, w.index, height_codes or {1: "in", 2: "cm"})
    wu = _unit_series(weight_unit, w.index, weight_codes or {1: "lb", 2: "kg"})
    hu = hu.where(hu.notna(), np.where(h > 100, "cm", "in"))
    wu = wu.where(wu.notna(), np.where(hu == "cm", "kg", "lb"))
    kg = np.where(wu == "lb", w * 0.45359237, w)
    m  = np.where(hu == "in", h * 0.0254, np.where(hu == "m", h, h / 100.0))
    bmi = pd.Series(kg / (m ** 2), index=w.index)
    return bmi.where(bmi.between(10, 80)).round(1)


def clean_uds_code(s, valid=(0, 1, 2)):
    """UDS health-history codes: keep 0/1/2 (absent / recent-active / remote);
    9 = unknown, -4 = not available, 8 = not applicable -> NaN."""
    v = pd.to_numeric(s, errors="coerce")
    return v.where(v.isin(valid))


def find_col(df, name):
    """Case-insensitive column lookup (OASIS mixes 'SEIZURES' / 'seizures')."""
    low = {c.lower(): c for c in df.columns}
    return low.get(name.lower())



# =============================================================================
# SEIZURE / EPILEPSY / ANTIEPILEPTIC TEXT HELPERS  (added 2026-11-08)
# =============================================================================
# Used where a cohort has no structured seizure field and only free-text
# medical-history or medication records (ADNI RECMHIST / INITHEALTH / RECCMEDS,
# ABVIB UDS A4, PREVENT-AD Med_use). Output is a per-subject evidence flag:
#     {name}_flag      1 = at least one matching record, 0 = records but no match
#     {name}_n         number of matching records
#     {name}_1st_date  date of the first matching record
#     {name}_evidence  up to 5 distinct matching text snippets
# A text flag is weaker evidence than a coded variable — keep it separate.
SEIZURE_TEXT_RX = re.compile(
    r"seiz|epilep|convuls|grand\s*mal|petit\s*mal|tonic.?clonic|status\s+epil|"
    r"absence\s+(attack|spell)", re.I)

# Drugs used mainly as antiepileptics. Deliberately EXCLUDES gabapentin,
# pregabalin, primidone and clonazepam (common non-epilepsy indications).
AED_TEXT_RX = re.compile(
    r"levetiracetam|keppra|brivaracetam|briviact|lamotrigin|lamictal|"
    r"valpro|depakot|depakene|divalproex|epival|carbamazepin|tegretol|"
    r"oxcarbazepin|trileptal|eslicarbazepin|aptiom|phenytoin|dilantin|"
    r"fosphenytoin|phenobarb|topiramat|topamax|lacosamid|vimpat|zonisamid|"
    r"zonegran|ethosuximid|zarontin|clobazam|onfi|frisium|perampanel|fycompa|"
    r"cenobamat|xcopri|rufinamid|felbamat|tiagabin|vigabatrin|stiripentol",
    re.I)


def text_flag_series(df: pd.DataFrame, id_col: str, text_cols, regex,
                     name: str, date_col: str | None = None,
                     extra_mask: pd.Series | None = None) -> pd.DataFrame:
    """Per-subject regex flag over one or more free-text columns."""
    cols = [c for c in text_cols if c in df.columns]
    if df.empty or id_col not in df.columns or not cols:
        return pd.DataFrame()
    txt = (df[cols].fillna("").astype(str)
           .apply(lambda s: s.where(~s.str.lower().isin(["nan", "none", "nat"]), ""))
           .agg(lambda r: " | ".join(x for x in r if x), axis=1))
    hit = txt.str.contains(regex, na=False)
    if extra_mask is not None:
        hit = hit | extra_mask.reindex(df.index).fillna(False).astype(bool)
    d = pd.DataFrame({"ID": df[id_col], "_hit": hit, "_txt": txt.str.strip(" |")})
    d["_t"] = (pd.to_datetime(df[date_col], errors="coerce")
               if date_col and date_col in df.columns else pd.NaT)
    d = d.dropna(subset=["ID"])
    rows = []
    for subj, g in d.groupby("ID"):
        h = g[g["_hit"]].sort_values("_t")
        ev = list(dict.fromkeys(h["_txt"].str[:80]))[:5]
        rows.append({
            "ID"               : subj,
            f"{name}_flag"     : int(len(h) > 0),
            f"{name}_n"        : int(len(h)),
            f"{name}_1st_date" : (str(h["_t"].iloc[0].date())
                                  if len(h) and pd.notna(h["_t"].iloc[0]) else np.nan),
            f"{name}_evidence" : " || ".join(ev) if ev else np.nan,
        })
    return pd.DataFrame(rows).set_index("ID") if rows else pd.DataFrame()
