"""Clinical diagnosis as two separate axes: cognitive STAGE and primary ETIOLOGY.

Why (v0.4.0)
------------
Up to v0.3.0 each pipeline built one compound label per visit ("AD+MCI",
"CN+EPY", "AD+VaD") by joining every active flag. In the cohorts that use the
UDS Form D1 (NACC, OASIS-3, ABVIB) the Alzheimer flags are *etiologies* that
clinicians mark at any level of impairment, so MCI due to AD was labelled
"AD"; impaired-not-MCI was labelled "MCI"; adding an etiology or epilepsy
counted as a change of diagnosis; and in OASIS-3 the UDS v3 fields
(``alzdis``, ``lbdis``, ``cvd`` ...) were not read at all.

From v0.4.0 every visit has

* ``stage``      CN | IMPNOMCI | MCI | DEM | Unknown   (syndromic severity)
* ``etiology``   primary etiology: AD, LBD, FTD, VaD, PSP, CBS, MSA, HD, PD,
                 Other, Epilepsy ... ; ``none`` for CN; ``unknown`` when impaired
                 without an etiology. Several primaries are joined with '+'.
* ``label``      harmonised clinical label (column ``all_diag``):
                 CN, IMPNOMCI, MCI, the dementia etiology (AD, LBD, ...),
                 Dementia-NOS or Dementia-Other; syndromic cohorts without a
                 stage (NIFD, 4RTNI) keep the syndrome (bvFTD, PSP ...).
* ``epy``        epilepsy recorded at the visit (covariate, not part of the label)

The trajectory (number and time of changes, conversions) is computed on the
STAGE only. The pre-0.4.0 compound label is kept in ``*_diag_legacy``.
"""
from __future__ import annotations

from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd

from .series import d_fmt, days_from, list_str

STAGES = ("CN", "IMPNOMCI", "MCI", "DEM")
KNOWN_STAGES = frozenset(STAGES)
STAGE_RANK = {s: i for i, s in enumerate(STAGES)}
NO_ETIOLOGY = ("", "none", "unknown")


# --------------------------------------------------------------------------- helpers
def _num(v: Any) -> float:
    """Value as float; NaN for missing / text."""
    if v is None:
        return np.nan
    try:
        x = float(v)
    except (TypeError, ValueError):
        return np.nan
    return x


def _get(row: pd.Series | dict, name: str) -> float:
    """Numeric value of ``name`` looked up as given, upper and lower case."""
    for key in (name, name.upper(), name.lower()):
        if key in row:
            v = _num(row[key])
            if not np.isnan(v):
                return v
    return np.nan


def _flag(row, name) -> bool:
    return _get(row, name) == 1


def join_etiologies(items: Iterable[str]) -> str:
    items = sorted({i for i in items if i and i not in NO_ETIOLOGY})
    return "+".join(items) if items else "unknown"


def harmonised_label(stage: str, etiology: str) -> str:
    """Clinical label of one visit from its stage and primary etiology."""
    etiology = "" if etiology is None or (isinstance(etiology, float) and np.isnan(etiology)) else str(etiology)
    if stage in ("CN", "IMPNOMCI", "MCI"):
        return stage
    known = etiology not in NO_ETIOLOGY
    if stage == "DEM":
        if not known:
            return "Dementia-NOS"
        return "Dementia-Other" if etiology in ("Other", "Epilepsy") else etiology
    # stage not recorded (syndromic cohorts): keep the syndrome
    if known:
        return "OtherDx" if etiology == "Other" else etiology
    return "Unknown"


def stage_of_label(label: Any) -> str:
    """Stage implied by a harmonised label (used only when no stage column exists)."""
    s = str(label).strip()
    if s in ("CN", "CU", "CON", "IMPNOMCI", "MCI"):
        return {"CU": "CN", "CON": "CN"}.get(s, s)
    if s.startswith("Dementia") or s in ("AD", "LBD", "VaD"):
        return "DEM"
    return "Unknown"


def add_dx_columns(dx: pd.DataFrame, stage_fn: Callable[[pd.Series], tuple],
                   legacy_col: str = "DX_LABEL") -> pd.DataFrame:
    """Add DX_STAGE, DX_ETIOLOGY, DX_EPY and the harmonised DX_LABEL.

    ``stage_fn(row) -> (stage, etiology, epy)``. The pre-0.4.0 label, when
    present in ``legacy_col``, is moved to ``DX_LEGACY``.
    """
    dx = dx.copy()
    dx["DX_LEGACY"] = dx[legacy_col] if legacy_col in dx.columns else np.nan
    if dx.empty:
        for c in ("DX_STAGE", "DX_ETIOLOGY", "DX_EPY", "DX_LABEL"):
            dx[c] = pd.Series(dtype=object)
        return dx
    res = dx.apply(lambda r: pd.Series(stage_fn(r), index=["DX_STAGE", "DX_ETIOLOGY", "DX_EPY"]), axis=1)
    dx["DX_STAGE"] = res["DX_STAGE"].fillna("Unknown").astype(str)
    dx["DX_ETIOLOGY"] = res["DX_ETIOLOGY"].fillna("unknown").astype(str)
    dx["DX_EPY"] = res["DX_EPY"]
    dx["DX_LABEL"] = [harmonised_label(s, e) for s, e in zip(dx["DX_STAGE"], dx["DX_ETIOLOGY"])]
    return dx


# --------------------------------------------------------------------------- trajectory
def trajectory_record(ptid, grp: pd.DataFrame, dates: list, *,
                      fmt: Callable[[Any], Any] = d_fmt,
                      delta: Callable[[Any, Any], float] = days_from) -> dict:
    """Subject-level diagnosis columns from its visits (already sorted in time).

    ``grp`` holds DX_STAGE, DX_ETIOLOGY, DX_LABEL, DX_LEGACY (and DX_EPY);
    ``dates`` are the visit times in the dataset's own axis; ``fmt`` formats a
    time for the pipe lists and ``delta(first, t)`` gives days from the first visit.
    Changes are counted on the STAGE (visits with Unknown stage are skipped).
    """
    stages = grp["DX_STAGE"].tolist()
    labels = grp["DX_LABEL"].tolist()
    etios = grp["DX_ETIOLOGY"].tolist()
    legacy = grp["DX_LEGACY"].tolist() if "DX_LEGACY" in grp.columns else [np.nan] * len(grp)
    first_date = dates[0] if dates else np.nan

    known = [i for i, s in enumerate(stages) if s in KNOWN_STAGES]
    ch_dates, ch_days, prev = [], [], None
    for i in known:
        if prev is not None and stages[i] != prev:
            ch_dates.append(dates[i])
            ch_days.append(delta(first_date, dates[i]))
        prev = stages[i]
    n_label_changes = sum(1 for a, b in zip(labels, labels[1:]) if a != b)
    i0 = known[0] if known else 0
    i1 = known[-1] if known else len(stages) - 1
    max_stage = max((stages[i] for i in known), key=STAGE_RANK.get) if known else "Unknown"

    rec = {
        "ID": ptid,
        "1st_diag": labels[0],
        "Lst_diag": labels[-1],
        "1st_diag_date": dates[0],
        "Lst_diag_date": dates[-1],
        "all_diag": list_str(labels),
        "all_diag_date": list_str([fmt(d) for d in dates]),
        "nb_of_visits": len(grp),
        "1st_time_of_change": ch_dates[0] if ch_dates else np.nan,
        "1st_time_of_change_days_to": ch_days[0] if ch_days else np.nan,
        "all_time_of_change": list_str([fmt(d) for d in ch_dates]),
        "all_time_of_change_days_to": list_str(ch_days),
        "Lst_time_of_change": ch_dates[-1] if ch_dates else np.nan,
        "nb_of_change": len(ch_dates),
        "Diag_total": list_str(sorted(set(labels))),
        # v0.4.0: stage / etiology axes
        "1st_stage": stages[i0] if known else "Unknown",
        "Lst_stage": stages[i1] if known else "Unknown",
        "max_stage": max_stage,
        "all_stage": list_str(stages),
        "1st_etiology": etios[i0],
        "Lst_etiology": etios[i1],
        "all_etiology": list_str(etios),
        "nb_of_label_change": n_label_changes,
        "1st_diag_legacy": legacy[0],
        "Lst_diag_legacy": legacy[-1],
        "all_diag_legacy": list_str(legacy),
    }
    if "DX_EPY" in grp.columns and grp["DX_EPY"].notna().any():
        epy = grp["DX_EPY"].fillna(False).astype(bool)
        rec["epilepsy_any"] = int(epy.any())
        rec["nb_of_epilepsy_visits"] = int(epy.sum())
    return rec


def static_record_extra(stage: str, etiology: str, legacy: Any) -> dict:
    """Stage/etiology columns for cohorts with one static diagnosis (A4, 4RTNI)."""
    return {"1st_stage": stage, "Lst_stage": stage, "max_stage": stage, "all_stage": stage,
            "1st_etiology": etiology, "Lst_etiology": etiology, "all_etiology": etiology,
            "nb_of_label_change": 0,
            "1st_diag_legacy": legacy, "Lst_diag_legacy": legacy, "all_diag_legacy": legacy}


# --------------------------------------------------------------------------- UDS Form D1
UDS_MCI_COLS = [
    "MCIAMEM", "MCIAPLUS", "MCIAPLAN", "MCIAPATT", "MCIAPEX", "MCIAPVIS",
    "MCINON1", "MCIN1LAN", "MCIN1ATT", "MCIN1EX", "MCIN1VIS",
    "MCINON2", "MCIN2LAN", "MCIN2ATT", "MCIN2EX", "MCIN2VIS",
]
#: etiology -> (presence flag, primary/contributing column); UDS v1.2/v2 and v3 names
UDS_ETIOLOGIES = [
    ("AD", "PROBAD", "PROBADIF"), ("AD", "POSSAD", "POSSADIF"), ("AD", "ALZDIS", "ALZDISIF"),
    ("LBD", "DLB", "DLBIF"), ("LBD", "LBDIS", "LBDIF"),
    ("VaD", "VASC", "VASCIF"), ("VaD", "VASCPS", "VASCPSIF"), ("VaD", "CVD", "CVDIF"),
    ("VaD", "STROKE", "STROKIF"),
    ("FTD", "FTD", "FTDIF"), ("FTD", "PPAPH", "PPAPHIF"),
    ("FTD", "FTLDMO", "FTLDMOIF"), ("FTD", "FTLDNOS", "FTLDNOIF"),
    ("PSP", "PSP", "PSPIF"), ("CBS", "CORT", "CORTIF"), ("MSA", "MSA", "MSAIF"),
    ("HD", "HUNT", "HUNTIF"), ("PD", "PARK", "PARKIF"),
    ("Other", "PRION", "PRIONIF"), ("Other", "DOWNS", "DOWNSIF"), ("Other", "HYCEPH", "HYCEPHIF"),
    ("Other", "BRNINJ", "BRNINJIF"), ("Other", "NEOP", "NEOPIF"), ("Other", "MEDS", "MEDSIF"),
    ("Other", "DYSILL", "DYSILLIF"), ("Other", "DEP", "DEPIF"), ("Other", "OTHPSY", "OTHPSYIF"),
    ("Other", "ALCDEM", "ALCDEMIF"), ("Other", "IMPSUB", "IMPSUBIF"), ("Other", "ESSTREM", "ESSTREIF"),
    ("Other", "COGOTH", "COGOTHIF"), ("Other", "COGOTH2", "COGOTH2F"), ("Other", "COGOTH3", "COGOTH3F"),
    ("Epilepsy", "EPILEP", "EPILEPIF"),
]


def uds_d1_stage(row) -> str:
    """Stage from the raw UDS D1 questions (any version)."""
    if _flag(row, "NORMCOG"):
        return "CN"
    if _flag(row, "DEMENTED"):
        return "DEM"
    if any(_flag(row, c) for c in UDS_MCI_COLS):
        return "MCI"
    if _flag(row, "IMPNOMCI"):
        return "IMPNOMCI"
    return "Unknown"


def uds_d1_epilepsy(row) -> bool:
    return _flag(row, "EPILEP") or _get(row, "EPILEPIF") in (1, 2, 3)


def uds_d1_stage_etiology(row) -> tuple:
    """(stage, primary etiology, epilepsy) from the raw UDS D1 questions.

    Primary = presence flag 1 and its IF column 1. When no etiology is marked
    primary, etiologies with a presence flag but no IF answer are used;
    contributing ones (IF = 2) never enter the label.
    """
    stage = uds_d1_stage(row)
    epy = uds_d1_epilepsy(row)
    if stage == "CN":
        return stage, "none", epy
    primary, undetermined = [], []
    for lab, flag, ifcol in UDS_ETIOLOGIES:
        f, lvl = _get(row, flag), _get(row, ifcol)
        if lvl == 1 and (f == 1 or np.isnan(f)):
            primary.append(lab)
        elif f == 1 and np.isnan(lvl) and lab != "Epilepsy":
            undetermined.append(lab)
    etiology = join_etiologies(primary or undetermined)
    return stage, etiology, epy


# --------------------------------------------------------------------------- NACC derived variables
NACCUDSD_STAGE = {1: "CN", 2: "IMPNOMCI", 3: "MCI", 4: "DEM"}
NACCETPR_ETIOLOGY = {1: "AD", 2: "LBD", 3: "MSA", 4: "PSP", 5: "CBS", 6: "FTD", 7: "FTD", 8: "VaD",
                     9: "Other", 10: "Other", 11: "HD", 12: "Other", 13: "Other", 14: "Other",
                     15: "Epilepsy", 16: "Other", 17: "Other", 18: "Other", 19: "Other", 20: "Other",
                     21: "Other", 22: "Other", 23: "Other", 24: "Other", 25: "Other", 26: "Other",
                     27: "Other", 28: "Other", 29: "Other", 30: "Other", 88: "none", 99: "unknown"}


def nacc_stage_etiology(row) -> tuple:
    """NACC: stage = NACCUDSD, etiology = NACCETPR (one primary etiology, UDS v1-v3).

    Falls back to the raw D1 questions when NACCUDSD is missing.
    """
    epy = uds_d1_epilepsy(row) or _get(row, "NACCETPR") == 15
    st = _get(row, "NACCUDSD")
    if np.isnan(st) or int(st) not in NACCUDSD_STAGE:
        stage, etiology, epy2 = uds_d1_stage_etiology(row)
        return stage, etiology, epy or epy2
    stage = NACCUDSD_STAGE[int(st)]
    if stage == "CN":
        return stage, "none", epy
    et = _get(row, "NACCETPR")
    etiology = NACCETPR_ETIOLOGY.get(int(et), "unknown") if not np.isnan(et) else "unknown"
    if etiology == "none":           # 88 with an impaired stage: inconsistent -> unknown
        etiology = "unknown"
    return stage, etiology, epy
