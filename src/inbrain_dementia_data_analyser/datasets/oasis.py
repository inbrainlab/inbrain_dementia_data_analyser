"""OASIS-3 LONGITUDINAL INTEGRATION PIPELINE
Applies the same domain-builder logic as ADNI v2, adapted for OASIS-3's
UDS-based multi-CSV structure with OASISID subject IDs and days_to_visit
time axis.

OASIS-3-SPECIFIC ARCHITECTURE
-------------------------------
OASIS-3 ships one CSV per UDS form, all keyed on:
  - OASISID            : Primary subject identifier → normalised to PTID
  - OASIS_session_label: Session label (e.g. "OAS3_0001_UDS_d1234")
  - days_to_visit      : Integer days from study entry (day 0 = entry date)
                         — primary time axis (same role as ADNI EXAMDATE)
  - age at visit       : Float age at each visit

AVAILABLE FILES (confirmed from PC Databank file tree,
                 all in Tables_OASIS/):
  OASIS3_UDSd1_diagnoses.csv           UDS D1 — diagnosis flags (NORMCOG, MCI, AD…)
  OASIS3_UDSb4_cdr.csv                 UDS B4 — CDR, MMSE (from CDR form), CDRSUM
  OASIS3_UDSb6_gds.csv                 UDS B6 — GDS (30 items)
  OASIS3_UDSc1_cognitive_assessments.csv  UDS C1 — full NPS battery
  OASIS3_UDSa1_participant_demo.csv    UDS A1 — living situation, marital status
  OASIS3_UDSb1_physical_eval.csv       UDS B1 — vitals (BP, weight, height)
  OASIS3_UDSb3.csv                     UDS B3 — motor/Parkinsonism (UPDRS-like)
  OASIS3_MR_json.csv                   MRI acquisitions (BIDS JSON metadata)
  OASIS3_PET_json.csv                  PET acquisitions (BIDS JSON metadata)
  OASIS3_amyloid_centiloid.csv         Amyloid PET centiloid values
  OASIS3_demographics.csv              Static demographics + APOE

NOTE: OASIS3_Clinical_and_Psychometrics_dictionary.xlsx is the data dictionary.

DIAGNOSIS CODING (UDS D1 — same as NACC/ABVIB convention):
  NORMCOG = 1 → CN     (takes sole priority)
  MCI sub-flags (MCIAMEM…MCINON2) → MCI
  PROBAD/POSSAD = 1 → AD
  DLB = 1 → LBD
  VASC/VASCPS = 1 → VaD
  FTD/PSP/CORT/PPAPH/PNAPH etc. = 1 → FTD
  PARK = 1 → PD
  DEMENTED = 1 (no subtype) → Dementia-NOS
  Compound labels via sorted '+' join (same as ADNI/AIBL/ABVIB)

TIME AXIS:
  days_to_visit = integer days from study entry (day 0000 = entry date).
  Stored directly as the "date" field in longitudinal columns.
  age at visit also available for cross-cohort temporal alignment.

MRI: OASIS3_MR_json.csv — BIDS JSON metadata per session
  Columns: OASIS_ID, session_id, filename, + scan parameters
  scan category (or filename) used to classify MRI type.

PET: OASIS3_PET_json.csv — BIDS JSON metadata per session
  Columns: subject_id, session_id, Tracer.Name / TracerName, Dosage,
           InjectedRadioactivity, etc.
  Tracer name → canonical via PET_TRACER_MAP.

AMYLOID CENTILOID: OASIS3_amyloid_centiloid.csv
  Processed PET metric — amyloid centiloid scale per session.

NO CSF table in OASIS-3 public release.

Output:
  <root>/Tables_OASIS/Integrated_OASIS/OASIS_longitudinal.csv

Author  : VRS_inbrain project
Version : 1  —  2026-06-09

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("OASIS", root_path)
"""
from __future__ import annotations

import gc
import numpy as np
import os
import pandas as pd
import re
from functools import reduce

from ..datasets.base import DatasetContext
from ..diagnosis import add_dx_columns, trajectory_record, uds_d1_stage_etiology
from ..imaging.scanner import add_scanner_channels
from ..imaging.sessions import session_id_columns
from ..series import build_scored_columns, list_str
from ..variables import (
    attach_generic_variables,
    load_generic_tables,
    materialise,
    variable_sources,
)

import logging as _logging
from ..logging_utils import make_print as _make_print

log = _logging.getLogger(__name__)
_print = _make_print(log)

DATASET = "OASIS"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
OASIS = OUT_DIR = OUT_LONG = ""


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, OASIS, OUT_DIR, OUT_LONG
    _CTX = ctx
    OASIS = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write the per-subject table of OASIS; returns it."""
    configure(ctx)
    return main()



def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns of one measure (see series.build_scored_columns)."""
    return build_scored_columns(scores_by_subj, test_name)




# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================











# MRI type map — applied to 'scan category' or filename
MRI_TYPE_MAP = {
    "t1w"        : "MRI_T1w",
    "mprage"     : "MRI_T1w",
    "t1_mprage"  : "MRI_T1w",
    "t2w"        : "MRI_T2w",
    "flair"      : "MRI_T2w_FLAIR",
    "t2_flair"   : "MRI_T2w_FLAIR",
    "swi"        : "MRI_T2star",
    "t2star"     : "MRI_T2star",
    "bold"       : "MRI_fMRI",
    "resting"    : "MRI_fMRI",
    "dwi"        : "MRI_DTI",
    "dti"        : "MRI_DTI",
    "asl"        : "MRI_ASL",
    "perfusion"  : "MRI_ASL",
}

# PET tracer map
PET_TRACER_MAP = {
    "pib"            : "ABETA_PIB",
    "c11-pib"        : "ABETA_PIB",
    "11c-pib"        : "ABETA_PIB",
    "fdg"            : "Metabolism",
    "18f-fdg"        : "Metabolism",
    "av45"           : "ABETA_AV45",
    "florbetapir"    : "ABETA_AV45",
    "18f-av45"       : "ABETA_AV45",
    "flortaucipir"   : "TAU_Flortaucipir",
    "av1451"         : "TAU_Flortaucipir",
    "18f-av1451"     : "TAU_Flortaucipir",
    "florbetaben"    : "ABETA_FBB",
    "18f-fbb"        : "ABETA_FBB",
    "mk6240"         : "TAU_MK6240",
    "18f-mk6240"     : "TAU_MK6240",
}

# UDS D1 diagnosis flag groups (same as ABVIB — UDS convention)
MCI_COLS = [
    "MCIAMEM", "MCIAPLUS", "MCIAPLAN", "MCIAPATT", "MCIAPEX", "MCIAPVIS",
    "MCINON1", "MCIN1LAN", "MCIN1ATT", "MCIN1EX",  "MCIN1VIS",
    "MCINON2", "MCIN2LAN", "MCIN2ATT", "MCIN2EX",  "MCIN2VIS", "IMPNOMCI",
]
AD_COLS  = ["PROBAD", "POSSAD"]
LBD_COLS = ["DLB"]
VAD_COLS = ["VASC", "VASCPS"]
FTD_COLS = ["FTD", "PSP", "CORT", "PPAPH", "PNAPH", "SEMDEMAN",
            "SEMDEMAG", "PPAOTHR"]
PD_COLS  = ["PARK"]

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================







def load(key: str) -> pd.DataFrame:
    path = PATHS.get(key, "")
    if not path or not os.path.exists(path):
        _print(f"  [SKIP] {key}: not found")
        return pd.DataFrame()
    for i, kwargs in enumerate([
        dict(engine="c",      on_bad_lines="error", encoding="utf-8"),
        dict(engine="c",      on_bad_lines="skip",  encoding="utf-8"),
        dict(engine="python", on_bad_lines="skip",  encoding="utf-8"),
        dict(engine="python", on_bad_lines="skip",  encoding="latin-1"),
        dict(engine="python", on_bad_lines="skip",  encoding="cp1252"),
    ]):
        try:
            df = pd.read_csv(path, low_memory=False, **kwargs)
            if i > 0:
                _print(f"  [WARN] {key}: loaded via fallback strategy {i+1}")
            _print(f"  [OK]   {key}: {len(df):>6,} rows  {df.shape[1]} cols")
            return df
        except Exception:
            continue
    _print(f"  [ERROR] {key}: all strategies failed")
    return pd.DataFrame()


def normalise_id(df: pd.DataFrame) -> pd.DataFrame:
    """Rename OASISID → PTID."""
    if "PTID" in df.columns:
        return df
    col_lower = {c.lower(): c for c in df.columns}
    for alt in ("oasisid", "oasis_id", "subject_id", "ptid"):
        if alt in col_lower:
            df = df.rename(columns={col_lower[alt]: "PTID"})
            return df
    return df


def sort_days(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure _days column for chronological sorting."""
    if "days_to_visit" in df.columns:
        df["_days"] = pd.to_numeric(df["days_to_visit"], errors="coerce").fillna(0)
    else:
        df["_days"] = 0
    return df





def _nps_from_df(df: pd.DataFrame, test_name: str,
                 score_col: str) -> pd.DataFrame:
    """Extract one test longitudinally using days_to_visit as time axis."""
    if df.empty or score_col not in df.columns:
        _print(f"    [WARN] '{score_col}' not found → {test_name} skipped")
        return pd.DataFrame()
    df = normalise_id(df.copy())
    df = sort_days(df)
    df = df.sort_values(["PTID", "_days"])
    scores_by_subj = {}
    for subj, grp in df.groupby("PTID"):
        scores_by_subj[subj] = list(zip(grp["_days"], grp[score_col]))
    return build_nps_columns(scores_by_subj, test_name)


# =============================================================================
# §1  DIAGNOSIS TRAJECTORY
# =============================================================================
# Source: OASIS3_UDSd1_diagnoses.csv — UDS D1 format
# Same flag system as NACC/ABVIB. Uses days_to_visit for ordering.
#
# Compound-label logic (identical to ADNI/AIBL/ABVIB):
#   NORMCOG=1 → sole priority "CN"
#   Collect all active impairment flags → join with '+'
# =============================================================================

def _classify_oasis(row: pd.Series) -> str:
    def active(cols):
        return any(
            pd.notna(row.get(c)) and float(row.get(c, 0) or 0) == 1
            for c in cols
        )

    # epilepsy (D1 'epilep' / 'epilepif') is filled even when cognition is normal
    epi = active(["EPILEP", "epilep"]) or any(
        pd.notna(row.get(c)) and float(row.get(c) or 0) in (1, 2)
        for c in ("EPILEPIF", "epilepif"))

    if active(["NORMCOG"]):
        return "CN+EPY" if epi else "CN"

    parts = []
    if epi:               parts.append("EPY")
    if active(MCI_COLS):  parts.append("MCI")
    if active(AD_COLS):   parts.append("AD")
    if active(LBD_COLS):  parts.append("LBD")
    if active(VAD_COLS):  parts.append("VaD")
    if active(FTD_COLS):  parts.append("FTD")
    if active(PD_COLS):   parts.append("PD")

    if not parts:
        if pd.notna(row.get("DEMENTED")) and float(row.get("DEMENTED", 0) or 0) == 1:
            return "Dementia-NOS"
        return "Unknown"

    return "+".join(sorted(parts))


def build_diagnosis(df_raw: pd.DataFrame) -> pd.DataFrame:
    _print("  §1  Diagnosis trajectory (UDS D1)...")
    dx = normalise_id(df_raw.copy())
    dx = sort_days(dx)
    dx["DX_LABEL"] = dx.apply(_classify_oasis, axis=1)   # pre-0.4.0 label -> *_diag_legacy
    dx = add_dx_columns(dx, uds_d1_stage_etiology)  # UDS D1 v1.2-v3 (upper/lower case)
    dx = dx.sort_values(["PTID", "_days"])

    records = []
    for ptid, grp in dx.groupby("PTID"):
        days = grp["_days"].tolist()
        records.append(trajectory_record(
            ptid, grp, days, fmt=lambda x: x,
            delta=lambda a, b: float(b - a) if pd.notna(a) and pd.notna(b) else np.nan))

    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §2  MRI ACQUISITIONS
# =============================================================================
# Source: OASIS3_MR_json.csv
# Columns confirmed: OASIS_ID, session_id, Name, License, filename, ...
# The 'scan category' column identifies MRI type.
# session_id = "OAS3_XXXX_MR_dYYYY" — encodes days_to_session (dYYYY part).
#
# Image ID = session_id (BIDS session identifier).
# Scan type inferred from 'scan category' (BIDS category) or filename.
# =============================================================================

def _classify_mri_type(desc: str) -> str | None:
    if not isinstance(desc, str):
        return None
    d = desc.lower().replace("-", "_").replace(" ", "_")
    # Skip infrastructure
    if any(x in d for x in ["scout", "localizer", "fieldmap", "field_map",
                              "gre_field", "aahscout", "mocoseries"]):
        return None
    # Skip derived DTI maps
    if any(d.endswith(s) for s in ["_fa", "_adc", "_tracew", "_trace"]):
        return None
    for key, canon in MRI_TYPE_MAP.items():
        if key in d:
            return canon
    return None


def _days_from_row(row: pd.Series) -> int | None:
    """Days since entry from any OASIS label: 'OAS30001_MR_d0129',
    'sub-OAS30001_ses-d0129_T1w.json', 'MR ID', 'label' …"""
    for c in ("session_id", "MR ID", "MR_ID", "label", "Label", "filename",
              "session", "Session", "experiment_id"):
        v = row.get(c, None)
        if isinstance(v, str):
            m = re.search(r"(?:_|ses-)d(\d{3,5})(?=\D|$)", v)
            if m:
                return int(m.group(1))
    return None


def _extract_days_from_session(session_id: str) -> int | None:
    """
    Extract days_to_session from OASIS session_id format.
    E.g. "OAS3_0001_MR_d0000" → 0
         "OAS3_0001_MR_d1234" → 1234
    """
    if not isinstance(session_id, str):
        return None
    m = re.search(r"_d(\d+)$", session_id)
    return int(m.group(1)) if m else None


def build_mri(df_raw: pd.DataFrame) -> pd.DataFrame:
    _print("  §2  MRI acquisitions (BIDS JSON metadata)...")
    if df_raw.empty:
        return pd.DataFrame()

    mri = df_raw.copy()
    # Normalise subject ID column (OASIS_ID in MR json)
    col_lower = {c.lower().replace(" ", "_"): c for c in mri.columns}
    for alt in ("oasis_id", "oasisid", "subject_id"):
        if alt in col_lower:
            mri = mri.rename(columns={col_lower[alt]: "PTID"})
            break

    if "PTID" not in mri.columns:
        _print("    [WARN] No subject ID column found in MRI file")
        return pd.DataFrame()

    # Extract scan type from 'scan category' or filename
    type_col = next((c for c in mri.columns
                     if "scan" in c.lower() and "cat" in c.lower()), None)
    if type_col is None:
        type_col = "filename"
    if type_col not in mri.columns:
        _print("    [WARN] No scan type column found in MRI file")
        return pd.DataFrame()

    mri["_canon_type"] = mri[type_col].apply(_classify_mri_type)
    mri = mri.dropna(subset=["_canon_type"])

    # Days since entry: session_id / label / filename (…_d0129 / ses-d0129).
    # Rows without a recoverable day stay NaN (NOT 0, which put every scan
    # on day 0 before).
    mri["_days"] = pd.to_numeric(mri.apply(_days_from_row, axis=1), errors="coerce")
    n_nod = int(mri["_days"].isna().sum())
    if n_nod:
        _print(f"    [WARN] {n_nod} MRI rows without a day code -> dropped")
    mri = mri.dropna(subset=["_days"])
    mri["_days"] = mri["_days"].astype(int)
    mri = mri.sort_values(["PTID", "_days"])

    records = {}
    avail_key = "MRI_type_disponible"

    for ptid, grp in mri.groupby("PTID"):
        row = records.setdefault(ptid, {"ID": ptid, avail_key: set()})
        for _, r in grp.iterrows():
            mtype = r["_canon_type"]
            day   = r["_days"]
            sid   = r.get("session_id", None)
            sid   = sid if isinstance(sid, str) else f"{ptid}_MR_d{day:04d}"
            row[avail_key].add(mtype)
            row.setdefault(f"all_{mtype}_date_list", []).append(day)
            row.setdefault(f"all_{mtype}_id_list",   []).append(str(sid))

    for ptid, row in records.items():
        row[avail_key] = "|".join(sorted(row[avail_key]))
        for mtype in row[avail_key].split("|"):
            dl = row.pop(f"all_{mtype}_date_list", [])
            il = row.pop(f"all_{mtype}_id_list",   [])
            # one ID per DATE (aligned lists); every ID of the date in all_{T}_image_ids
            session_id_columns(row, mtype, dl, il)

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(list(records.values()))
    if "ID" not in result.columns:
        return pd.DataFrame()
    return add_scanner_channels(result.set_index("ID"), mri, "PTID", "_canon_type", "_days")


# =============================================================================
# §3  PET ACQUISITIONS
# =============================================================================
# Source: OASIS3_PET_json.csv + OASIS3_amyloid_centiloid.csv
#
# PET_json: subject_id, session_id, Tracer.Name / TracerName, + parameters
# Centiloid: session-level amyloid centiloid values
#
# Tracer normalisation: same PET_TRACER_MAP as ADNI v2.
# session_id encodes days (OAS3_XXXX_PETYYYYMMDD or OAS3_XXXX_PET_dNNN).
# =============================================================================

# Centiloid column of OASIS3_amyloid_centiloid.csv written per PET session
# (all_{tracer}_centiloid). The *_fSUVR_rsf_* value exists for PIB and AV45;
# the *_fBP_* columns exist only for PIB (dynamic scans).
OASIS_CENTILOID_COL = "Centiloid_fSUVR_rsf_TOT_CORTMEAN"


def _norm_tracer_oasis(raw) -> str | None:
    if pd.isna(raw):
        return None
    t = str(raw).strip().lower()
    if t in ("", "nan", "none"):
        return None
    for key, canon in PET_TRACER_MAP.items():
        if key in t:
            return canon
    return None


def build_pet(df_pet: pd.DataFrame,
              df_centiloid: pd.DataFrame) -> pd.DataFrame:
    _print("  §3  PET acquisitions...")
    records = {}

    if not df_pet.empty:
        pet = df_pet.copy()
        # Normalise subject ID
        col_lower = {c.lower().replace(".", "_"): c for c in pet.columns}
        for alt in ("subject_id", "oasisid", "oasis_id"):
            if alt in col_lower:
                pet = pet.rename(columns={col_lower[alt]: "PTID"})
                break
        if "PTID" not in pet.columns:
            _print("    [WARN] No subject ID column in PET file")
        else:
            # Resolve tracer column (Tracer.Name / TracerName)
            tracer_col = None
            if "Tracer.Name" in pet.columns and "TracerName" in pet.columns:
                pet["_tracer_raw"] = pet["Tracer.Name"].fillna(pet["TracerName"])
                tracer_col = "_tracer_raw"
            elif "Tracer.Name" in pet.columns:
                tracer_col = "Tracer.Name"
            elif "TracerName" in pet.columns:
                tracer_col = "TracerName"
            else:
                for c in pet.columns:
                    if "tracer" in c.lower():
                        tracer_col = c
                        break

            if tracer_col:
                pet["_tracer"] = pet[tracer_col].apply(_norm_tracer_oasis)
                pet = pet.dropna(subset=["_tracer"])

                # Days from session_id
                pet["_days"] = pd.to_numeric(pet.apply(_days_from_row, axis=1),
                                              errors="coerce")
                pet = pet.dropna(subset=["_days"])
                pet["_days"] = pet["_days"].astype(int)
                pet = pet.sort_values(["PTID", "_days"])

                # centiloid of each PET session (amyloid tracers), keyed by session id
                cl_by_session = {}
                if not df_centiloid.empty:
                    _cl = df_centiloid
                    _sid = next((c for c in _cl.columns
                                 if c.lower() in ("oasis_session_id", "session_id")), None)
                    _val = OASIS_CENTILOID_COL if OASIS_CENTILOID_COL in _cl.columns else \
                        next((c for c in _cl.columns if "centiloid" in c.lower()), None)
                    if _sid and _val:
                        cl_by_session = (_cl.dropna(subset=[_val])
                                            .drop_duplicates(_sid).set_index(_sid)[_val].to_dict())

                for ptid, grp in pet.groupby("PTID"):
                    row = records.setdefault(ptid, {
                        "ID": ptid, "PET_tracer_disponible": set()
                    })
                    for tracer, tgrp in grp.groupby("_tracer"):
                        row["PET_tracer_disponible"].add(tracer)
                        days = tgrp["_days"].tolist()
                        ids  = tgrp["session_id"].tolist() if "session_id" in tgrp.columns \
                               else [None] * len(tgrp)
                        session_id_columns(row, tracer, days, ids)
                        days = list(dict.fromkeys(days))
                        row[f"1st_{tracer}_date"]   = days[0]  if days else np.nan
                        row[f"Lst_{tracer}_date"]   = days[-1] if days else np.nan
                        # all_{tracer}_centiloid: one value per session, aligned
                        if cl_by_session and isinstance(row.get(f"{tracer}_IDs"), str):
                            cls = [cl_by_session.get(i) for i in row[f"{tracer}_IDs"].split("|")]
                            if any(v is not None for v in cls):
                                row[f"all_{tracer}_centiloid"] = "|".join(
                                    "NA" if v is None else f"{round(float(v), 4):g}" for v in cls)
                        # Manufacturer if available
                        if "Manufacturer" in tgrp.columns:
                            mfr = tgrp["Manufacturer"].dropna().unique().tolist()
                            row[f"{tracer}_manufacturer"] = list_str(mfr)

    # ── Amyloid centiloid (processed metric) ─────────────────────────────────
    if not df_centiloid.empty:
        cl = df_centiloid.copy()
        col_lower = {c.lower(): c for c in cl.columns}
        for alt in ("oasisid", "oasis_id", "subject_id"):
            if alt in col_lower:
                cl = cl.rename(columns={col_lower[alt]: "PTID"})
                break
        if "PTID" in cl.columns:
            cl_col = next((c for c in cl.columns
                           if "centiloid" in c.lower()), None)
            if cl_col:
                for ptid, grp in cl.groupby("PTID"):
                    row = records.setdefault(ptid, {
                        "ID": ptid, "PET_tracer_disponible": set()
                    })
                    vals = grp[cl_col].dropna().tolist()
                    row["all_PIB_centiloid"] = list_str(vals)
                    row["1st_PIB_centiloid"] = vals[0]  if vals else np.nan
                    row["Lst_PIB_centiloid"] = vals[-1] if vals else np.nan

    for r in records.values():
        r["PET_tracer_disponible"] = (
            "|".join(sorted(r["PET_tracer_disponible"])) or np.nan
        )

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(list(records.values()))
    if "ID" not in result.columns:
        return pd.DataFrame()
    result = result.set_index("ID")
    if not df_pet.empty and "pet" in locals() and "_tracer" in pet.columns:
        result = add_scanner_channels(result, pet, "PTID", "_tracer", "_days")
    return result


# =============================================================================
# §4  NEUROPSYCHOLOGICAL + CDR + GDS + VITALS + MOTOR
# =============================================================================
# All UDS-format: OASISID, OASIS_session_label, days_to_visit, "age at visit"
#
# NPS MAPPING (confirmed from column headers in file tree):
# ┌────────────────────┬──────────────────────────────────────────────────────┐
# │ Source             │ Raw → Canonical                                       │
# ├────────────────────┼──────────────────────────────────────────────────────┤
# │ UDSb4_cdr          │ MMSE → MMSE, CDRTOT → CDR, CDRSUM → CDR_SB          │
# │                    │ memory/orient/judgment/commun/homehobb/perscare       │
# │ UDSb6_gds          │ GDS → GDS                                             │
# │ UDSc1_cognitive    │ ANIMALS → VF_anim, VEG → VF_veg                      │
# │                    │ DIGIF/DIGIB → DS_Fwd/DS_Back                          │
# │                    │ digsym → DigSym                                        │
# │                    │ tma → TMT_A, tmb → TMT_B                               │
# │                    │ bnt → BNT                                              │
# │                    │ LOGIMEM → LM_I, MEMUNITS → LM_II (or lmdelay)        │
# │                    │ mocatots → MoCA                                        │
# │                    │ TRAILA → TMT_A_UDS (alternate)                        │
# │ UDSb1_vitals       │ WEIGHT → Weight, HEIGHT → Height                      │
# │                    │ BPSYS → BP_Sys, BPDIAS → BP_Dias, HRATE → Pulse      │
# │ UDSa1_demo         │ MARISTAT → marital_status                             │
# │                    │ LIVSIT/LIVSITUA → living_situation                    │
# │                    │ INDEPEND → independence_level                          │
# │ UDSb3_motor        │ PDNORMAL → UPDRS_PDNormal (motor exam)                │
# │                    │ BRADYKIN → UPDRS_Brady                                 │
# └────────────────────┴──────────────────────────────────────────────────────┘
# =============================================================================



# UDSa1 non-numeric columns (categorical longitudinal)



def build_all_clinical(dfs: dict) -> pd.DataFrame:
    """Clinical / NPS / vitals / health-history / categorical series.

    Variables come from ``variables.yaml`` (``sources: OASIS:``); every series
    is timed with ``days_to_visit``.
    """
    _print("  §4  Clinical / NPS / Vitals / Functional...")
    if dfs.get("health", pd.DataFrame()).empty:
        _print("    [INFO] OASIS3_UDSa5 (health history) not found -> SEIZURES / DIABETES "
              "skipped. Download it from OASIS-3 XNAT into Tables_OASIS/.")
    parts = []
    for name, src in variable_sources(DATASET, config=_CTX.config):
        df = dfs.get(src.table, pd.DataFrame())
        if df.empty:
            continue
        frame, col = materialise(df, name, src)
        if col is None:
            continue
        part = _nps_from_df(frame, name, col)
        if not part.empty:
            parts.append(part)

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()
    return reduce(lambda a, b: a.join(b, how="outer"), parts)


# =============================================================================
# §5  DEMOGRAPHICS + APOE
# =============================================================================
# Source: OASIS3_demographics.csv (static per subject)
# Confirmed columns: OASISID, Subject_accession, AgeatEntry, AgeatDeath,
#   GENDER (1=M 2=F), EDUC, SES, racecode, race, ETHNIC,
#   AIAN, NHPI, ASIAN, AA, WHITE, daddem, momdem, HAND, APOE
#
# APOE: genotype string (format not specified in dict — likely "e3/e4")
#   APOE_e4_count = count of '4' in string
#
# Family history: daddem (father dementia), momdem (mother dementia)
# =============================================================================

def build_demographics(df_demog: pd.DataFrame) -> pd.DataFrame:
    _print("  §5  Demographics + APOE...")
    if df_demog.empty:
        return pd.DataFrame()

    demog = normalise_id(df_demog.copy())
    demog = demog.drop_duplicates("PTID")

    records = []
    for _, row in demog.iterrows():
        geno = str(row.get("APOE", "") or "").strip()
        records.append({
            "ID"                 : row["PTID"],
            "gender"             : row.get("GENDER",       np.nan),
            "education"          : row.get("EDUC",         np.nan),
            "ses"                : row.get("SES",           np.nan),
            "race"               : row.get("racecode",      np.nan),
            "race_label"         : row.get("race",          np.nan),
            "hispanic"           : row.get("ETHNIC",        np.nan),
            "handedness"         : row.get("HAND",          np.nan),
            "age_at_entry"       : row.get("AgeatEntry",    np.nan),
            "age_at_death"       : row.get("AgeatDeath",    np.nan),
            "famhx_father_dem"   : row.get("daddem",        np.nan),
            "famhx_mother_dem"   : row.get("momdem",        np.nan),
            "APOE_genotype"      : geno if geno and geno.lower() != "nan" else np.nan,
            "APOE_e4_count"      : geno.count("4") if geno else np.nan,
        })

    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §6  AGE AT FIRST DIAGNOSIS  +  BIRTH DATE RECONSTRUCTION
# =============================================================================
# OASIS uses days_to_visit + AgeatEntry to approximate age at 1st diagnosis.
# age_at_1st_diag = AgeatEntry + (1st_diag_date / 365.25)
# (1st_diag_date is days_to_visit at first dx visit)
#
# No birth date available directly. Estimated birth year only:
#   birth_year ≈ entry_year − AgeatEntry   (entry_year not provided)
#   → birth_date stored as NaN (not available in OASIS-3 public)
# =============================================================================

def add_age_at_1st_diag(master: pd.DataFrame) -> pd.DataFrame:
    if "age_at_entry" not in master.columns:
        return master
    age_entry = pd.to_numeric(master["age_at_entry"], errors="coerce")
    days_1st  = pd.to_numeric(master.get("1st_diag_date", 0), errors="coerce").fillna(0)
    master["age_at_1st_diag"] = (age_entry + days_1st / 365.25).round(2)
    bad = (master["age_at_1st_diag"] < 0) | (master["age_at_1st_diag"] > 120)
    master.loc[bad, "age_at_1st_diag"] = np.nan
    master["birth_date"] = np.nan   # not available in OASIS-3 public release
    return master


# =============================================================================
# MAIN
# =============================================================================

def main() -> pd.DataFrame:
    _print("=" * 62)
    _print("OASIS-3 LONGITUDINAL INTEGRATION PIPELINE  v1")
    _print("=" * 62)
    _print("\nSubject ID : OASISID → PTID")
    _print("Time axis  : days_to_visit (days from study entry, day 0 = entry)")
    _print("Diagnosis  : UDS D1 flags (same logic as NACC/ABVIB)")
    _print("No CSF protein table in OASIS-3 public release.")

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load ──────────────────────────────────────────────────────────────────
    _print("\nLoading raw tables...")
    dfs = {k: load(k) for k in PATHS}

    # ── Build domains ─────────────────────────────────────────────────────────
    _print("\nBuilding domain tables...")
    diag  = build_diagnosis(dfs["dx"])
    mri   = build_mri(dfs["mri"])
    pet   = build_pet(dfs["pet"], dfs["centiloid"])
    clin  = build_all_clinical(dfs)
    demog = build_demographics(dfs["demog"])

    # ── Merge ─────────────────────────────────────────────────────────────────
    _print("\nMerging all domains on subject ID...")
    parts = [p for p in [diag, mri, pet, clin, demog]
             if isinstance(p, pd.DataFrame) and not p.empty]

    if not parts:
        _print("ERROR: no domain tables produced.")
        return pd.DataFrame()

    master = reduce(lambda a, b: a.join(b, how="outer"), parts)
    master.index.name = "ID"
    master = master.reset_index()
    gc.collect()

    # ── Post-processing ───────────────────────────────────────────────────────
    master = add_age_at_1st_diag(master)

    if "MRI_type_disponible" not in master.columns:
        master["MRI_type_disponible"] = np.nan
    if "PET_tracer_disponible" not in master.columns:
        master["PET_tracer_disponible"] = np.nan
    master["CSF_available"] = "not_available_in_OASIS3"

    master["Plasma_available"] = "not_available_in_OASIS3"

    LEAD = [
        "ID",
        "gender", "race", "race_label", "hispanic", "education", "ses",
        "handedness", "age_at_entry", "age_at_death",
        "age_at_1st_diag", "birth_date",
        "APOE_genotype", "APOE_e4_count",
        "famhx_father_dem", "famhx_mother_dem",
        "1st_diag", "Lst_diag", "1st_diag_date", "Lst_diag_date",
        "all_diag", "all_diag_date", "nb_of_visits",
        "1st_time_of_change", "1st_time_of_change_days_to",
        "all_time_of_change", "all_time_of_change_days_to",
        "Lst_time_of_change", "nb_of_change", "Diag_total",
        "MRI_type_disponible", "PET_tracer_disponible", "CSF_available",
    ]
    LEAD += ["Plasma_available", "Plasma_analytes_disponible", "Plasma_panels_disponible"]
    exist_lead = [c for c in LEAD if c in master.columns]
    other_cols = [c for c in master.columns if c not in exist_lead]
    master     = master[exist_lead + other_cols]

    # ── Save ──────────────────────────────────────────────────────────────────
    master = attach_generic_variables(master, DATASET, load_generic_tables(DATASET, _CTX),
                                      _CTX.config)
    master.to_csv(OUT_LONG, index=False)

    _print("\n" + "=" * 62)
    _print(f"DONE  {len(master):,} subjects  ×  {len(master.columns):,} columns")
    _print(f"Saved → {OUT_LONG}")
    _print("=" * 62)

    # ── Summary ───────────────────────────────────────────────────────────────
    if "1st_diag" in master.columns:
        _print("\nDiagnosis distribution (1st visit):")
        _print(master["1st_diag"].value_counts().to_string())

    if "nb_of_change" in master.columns:
        _print(f"\nSubjects with ≥1 diagnosis change: "
              f"{master['nb_of_change'].gt(0).sum():,}")

    if "MRI_type_disponible" in master.columns:
        from collections import Counter
        cnt = Counter()
        for v in master["MRI_type_disponible"].dropna():
            cnt.update(str(v).split("|"))
        if cnt:
            _print("\nMRI types:")
            for t, n in cnt.most_common():
                _print(f"  {t:<30s} {n:>5,} subjects")

    if "PET_tracer_disponible" in master.columns:
        from collections import Counter
        cnt = Counter()
        for v in master["PET_tracer_disponible"].dropna():
            cnt.update(str(v).split("|"))
        if cnt:
            _print("\nPET tracers:")
            for t, n in cnt.most_common():
                _print(f"  {t:<30s} {n:>5,} subjects")

    _print("\nTest coverage (subjects with ≥1 valid score):")
    for name in list(dict.fromkeys(n for n, _ in variable_sources(DATASET, config=_CTX.config)))[:12]:
        col = f"nb_of_{name}_visits"
        if col in master.columns:
            n = master[col].gt(0).sum()
            if n > 0:
                _print(f"  {name:<28s} {n:>5,} subjects")

    return master


# =============================================================================
