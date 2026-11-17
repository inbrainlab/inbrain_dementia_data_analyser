"""Scanner manufacturer and model per acquisition, aligned with ``all_{T}_date``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import logging as _logging
from ..logging_utils import make_print as _make_print

log = _logging.getLogger(__name__)
_print = _make_print(log)

# =============================================================================
# SCANNER MANUFACTURER / MODEL — per acquisition, aligned with all_{T}_date
# =============================================================================
MFR_COL_CANDIDATES   = ("Manufacturer", "SCANNERMANUFACTURER", "mri_mfr", "pet_mfr",
                        "mfr", "manufacturer", "scanner_manufacturer")
MODEL_COL_CANDIDATES = ("ManufacturersModelName", "SCANNERMODEL", "mri_mfr_model",
                        "pet_mfr_model", "mfr_model", "model", "scanner_model")


def norm_vendor(x):
    """'SIEMENS' / 'Siemens Healthineers' -> Siemens ; 'GE MEDICAL SYSTEMS' -> GE ;
    'Philips Healthcare' -> Philips ; Toshiba / Canon -> Canon."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return np.nan
    s = str(x).strip()
    if s.lower() in ("", "nan", "none", "-4", "unknown"):
        return np.nan
    l = s.lower()
    if "siemens" in l or l.startswith(("biograph", "ecat")):
        return "Siemens"
    if l in ("ge", "g.e.") or "general electric" in l or l.startswith("ge ") or "ge medical" in l \
            or l.startswith(("discovery", "signa")):
        return "GE"
    if "philips" in l or l.startswith(("ingenuity", "vereos", "gemini", "achieva", "ingenia")):
        return "Philips"
    if "united imaging" in l or l.startswith(("uih", "umi", "umr")):
        return "United Imaging"
    if "toshiba" in l or "canon" in l:
        return "Canon"
    if "hitachi" in l or "fujifilm" in l:
        return "Hitachi/Fujifilm"
    return s


def _tok(v):
    """Same token the builders write in all_{T}_date / all_{T}_visit."""
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return None
    return str(v.date()) if hasattr(v, "date") else str(v)


def add_scanner_channels(result: pd.DataFrame, raw: pd.DataFrame,
                         subj_col: str, type_col: str, time_col: str,
                         mfr_col: str | None = None, model_col: str | None = None,
                         list_suffix: str = "date") -> pd.DataFrame:
    """
    For every acquisition type T in `result`, add (aligned token-by-token with
    all_{T}_{list_suffix}):
        all_{T}_manufacturer   vendor of the scanner used at that acquisition
        all_{T}_model          scanner model
        {T}_manufacturer       vendor(s) used across all acquisitions
    Same-day acquisitions on different scanners are joined with '/'.
    Missing -> 'NA'. Columns are auto-detected when not given.
    """
    if result is None or result.empty or raw is None or raw.empty:
        return result
    # explicit names that are absent from `raw` fall back to auto-detection
    # (avoids KeyError when an export lacks e.g. SCANNERMANUFACTURER)
    if mfr_col and mfr_col not in raw.columns:
        mfr_col = None
    if model_col and model_col not in raw.columns:
        model_col = None
    mfr_col   = mfr_col   or next((c for c in MFR_COL_CANDIDATES   if c in raw.columns), None)
    model_col = model_col or next((c for c in MODEL_COL_CANDIDATES if c in raw.columns), None)
    if not mfr_col and not model_col:
        return result
    r = pd.DataFrame({
        "S": raw[subj_col].values, "T": raw[type_col].values,
        "K": [_tok(v) for v in raw[time_col]],
        "M": raw[mfr_col].map(norm_vendor).values if mfr_col else np.nan,
        "D": raw[model_col].astype(str).str.strip().where(raw[model_col].notna()).values
             if model_col else np.nan,
    }).dropna(subset=["S", "T", "K"])
    # same-day scans on >1 scanner: keep vendor and model PAIRED (same order)
    lut_m, lut_d = {}, {}
    for key, g in r.groupby(["S", "T", "K"], sort=False):
        pairs = sorted({("NA" if pd.isna(m) else str(m), "NA" if pd.isna(d) else str(d))
                        for m, d in zip(g["M"], g["D"])})
        if len(pairs) > 1:                       # drop all-NA pair if a real one exists
            pairs = [p for p in pairs if p != ("NA", "NA")] or pairs
        lut_m[key] = "/".join(p[0] for p in pairs)
        lut_d[key] = "/".join(p[1] for p in pairs)

    new = {}
    for T in r["T"].dropna().unique():
        lc = f"all_{T}_{list_suffix}"
        if lc not in result.columns:
            continue
        man, mod, summ = [], [], []
        for sid, cell in zip(result.index, result[lc]):
            if not isinstance(cell, str) or not cell.strip():
                man.append(np.nan); mod.append(np.nan); summ.append(np.nan); continue
            toks = cell.split("|")
            m = [lut_m.get((sid, T, t), "NA") for t in toks]
            d = [lut_d.get((sid, T, t), "NA") for t in toks]
            man.append("|".join(m) if mfr_col else np.nan)
            mod.append("|".join(d) if model_col else np.nan)
            vend = sorted({x for v in m for x in v.split("/") if x != "NA"})
            summ.append("|".join(vend) if vend else np.nan)
        if mfr_col:
            new[f"all_{T}_manufacturer"] = man
            new[f"{T}_manufacturer"]     = summ
        if model_col:
            new[f"all_{T}_model"] = mod
    if new:
        add = pd.DataFrame(new, index=result.index)
        result = result.drop(columns=[c for c in add.columns if c in result.columns])
        result = result.join(add)
    return result
