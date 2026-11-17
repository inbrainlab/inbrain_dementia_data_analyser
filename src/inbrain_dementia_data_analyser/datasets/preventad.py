"""PREVENT-AD LONGITUDINAL INTEGRATION PIPELINE
Applies the same domain-builder logic as ADNI v2, adapted for PREVENT-AD's
multi-CSV structure with CONP_CandID subject IDs and Visit_label time axis.

PREVENT-AD-SPECIFIC ARCHITECTURE
----------------------------------
PREVENT-AD ships one CSV per instrument, all keyed on:
  - CONP_CandID : Primary subject identifier → normalised to PTID
  - CONP_ID     : Alternate ID (Open data repository code)
  - Visit_label : Visit time point (e.g. "PREBL00", "PRENA12", "NAP00")
  - Study_visit_label : Full label with prefix (e.g. "PRE-BL00")
  - Date_taken  : YYYY-MM (no day — de-identified to month precision)
  - Candidate_Age: Age in MONTHS at visit

VISIT LABEL STRUCTURE:
  Prefix PRE = enrolled in PREVENT-AD main study
  Prefix NAP = enrolled in NAP sub-study (Naproxen prevention sub-arm)
  Suffix encodes visit: BL00=baseline, NA12=12mo, NA24=24mo, NA36=36mo,
                        NA48=48mo, FU00=followup, etc.
  Visit ordering uses Candidate_Age (months) when date is partially masked.

DIAGNOSIS CODING:
  Two sources — use together:
  1. MCI_Registered: RC_MCI ("CU" or "MCI") + RC_MCI_subtype (aMCI SD/MD,
     naMCI SD/MD) — research classification per visit.
  2. Clinical_diagnosis_Registered: Clinical_diagnosis (free text from
     external clinic / caregiver report) — not per regular visit schedule.

  Pipeline logic:
    - Primary = MCI RC_MCI column (longitudinal, research-grade)
    - Supplement with Clinical_diagnosis (static, appended when present)
    - CU → CN (cognitively unimpaired = cognitively normal for harmonisation)
    - MCI → MCI
    - Clinical_diagnosis strings map to canonical labels

AVAILABLE FILES (confirmed from PC Databank file tree,
                 all in Tables_PREVENT-AD/):
  AD8_Registered_PREVENTAD_2025-03-21.csv
  APS_Registered_PREVENTAD_2025-03-21.csv
  Auditory_processing_Registered_PREVENTAD_2025-03-21.csv
  BP_Pulse_Weight_Registered_PREVENTAD_2025-03-21.csv
  CDR_FU_Registered_PREVENTAD_2025-03-21.csv
  Clinical_diagnosis_Registered_PREVENTAD_2025-03-21.csv
  CSF_proteins_Registered_PREVENTAD_2025-03-21.csv
  CSF_SomaScan_7K_proteins_Registered_PREVENTAD_2025-03-21.csv
  DKEFS-CWIT_Registered_PREVENTAD_2024-08-01.csv
  Demographics_Registered_PREVENTAD_2025-03-21.csv
  EL_CAIDE_Registered_PREVENTAD_2025-03-21.csv
  EL_CDR_MoCA_Registered_PREVENTAD_2025-03-21.csv
  EL_Medical_history_Registered_PREVENTAD_2025-03-21.csv
  Genetics_Registered_PREVENTAD_2025-03-21.csv
  Lab_Registered_PREVENTAD_2025-03-21.csv
  MCI_Registered_PREVENTAD_2025-03-21.csv
  Med_use_Registered_PREVENTAD_2025-03-21.csv
  PET_info_participants_Registered_PREVENTAD_2024-08-01.csv
  Plasma_4plex_IPMS_Registered_PREVENTAD_2025-03-21.csv
  Plasma_p-tau217_Registered_PREVENTAD_2025-03-21.csv
  RAVLT_Registered_PREVENTAD_2025-03-21.csv
  RBANS_Registered_PREVENTAD_2025-03-21.csv
  Smell_identification_Registered_PREVENTAD_2025-03-21.csv
  TMT_Registered_PREVENTAD_2025-03-21.csv

PET TRACER:
  PET_tracer column: "NAV4694" (18F-NAV4694 = amyloid)
  → canonical: ABETA_NAV4694

CSF PROTEINS:
  tau, ptau, Amyloid_beta_1_42, Amyloid_beta_1_40_lumipulse,
  Amyloid_beta_1_42_lumipulse, ptau217_simoa

PLASMA BIOMARKERS:
  Plasma_4plex: AB_1_40, AB_1_42, AB_ratio, t_tau, ptau181, ptau231,
                NFL, GFAP (IPMS + Simoa platforms)
  Plasma_p-tau217: ptau217_simoa_UGOT

MRI:
  No raw MRI acquisition table in the public release.
  PET_info has Freesurfer_visit_label + Date_MRI (the T1 used for PET
  processing) → used to flag T1w MRI presence.
  Processed MRI metrics in data dictionaries for ASL (CBF), WMH volumes
  — files referenced in data dictionaries but not in the file tree
  (not available in the registered release).

Output:
  <root>/Tables_PREVENT-AD/Integrated_PREVENTAD/PREVENTAD_longitudinal.csv

Author  : VRS_inbrain project
Version : 1  —  2026-06-09

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("PREVENTAD", root_path)
"""
from __future__ import annotations

import gc
import numpy as np
import os
import pandas as pd
from functools import reduce

from ..clinical import AED_TEXT_RX, compute_bmi, text_flag_series
from ..datasets.base import DatasetContext
from ..fluid import plasma_add_summary, plasma_clean_numeric
from ..imaging.sessions import session_id_columns
from ..series import build_scored_columns, d_fmt, list_str, safe_date
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

DATASET = "PREVENTAD"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
PREVAD = OUT_DIR = OUT_LONG = ""


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, PREVAD, OUT_DIR, OUT_LONG
    _CTX = ctx
    PREVAD = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write the per-subject table of PREVENTAD; returns it."""
    configure(ctx)
    return main()



def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns of one measure (see series.build_scored_columns)."""
    return build_scored_columns(scores_by_subj, test_name)




# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================











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
    """
    Rename CONP_CandID → PTID.
    CONP_CandID is the primary per-subject identifier in PREVENT-AD.
    CONP_ID is a secondary (open-data) identifier — kept as a column.
    """
    if "PTID" in df.columns:
        return df
    col_lower = {c.lower(): c for c in df.columns}
    for alt in ("conp_candid", "conp_id", "ptid", "subject_id"):
        if alt in col_lower:
            df = df.rename(columns={col_lower[alt]: "PTID"})
            return df
    return df


def sort_key(df: pd.DataFrame) -> pd.Series:
    """
    Return a numeric sort key. Priority:
      1. Candidate_Age (months) — precise within-subject ordering
      2. Date_taken (YYYY-MM) — partial date
      3. Visit_label lexicographic (BL00 < NA12 < NA24 …)
    """
    if "Candidate_Age" in df.columns:
        return pd.to_numeric(df["Candidate_Age"], errors="coerce").fillna(9999)
    if "Date_taken" in df.columns:
        return pd.to_datetime(df["Date_taken"], errors="coerce")\
                 .astype("int64").div(1e9).fillna(9999)
    if "Visit_label" in df.columns:
        return df["Visit_label"].fillna("ZZZ").astype(str)
    return pd.Series(0, index=df.index)





def _nps_from_df(df: pd.DataFrame, test_name: str,
                 score_col: str,
                 date_col: str = "Date_taken") -> pd.DataFrame:
    """Extract one test longitudinally. Uses Candidate_Age as sort key."""
    if df.empty or score_col not in df.columns:
        _print(f"    [WARN] '{score_col}' not found → {test_name} skipped")
        return pd.DataFrame()
    df = normalise_id(df.copy())
    df["_sort"] = sort_key(df)
    safe_date(df, date_col)
    df = df.sort_values(["PTID", "_sort"])
    date_series = df[date_col] if date_col in df.columns else df.get("Date_taken")
    scores_by_subj = {}
    for subj, grp in df.groupby("PTID"):
        scores_by_subj[subj] = list(zip(
            grp[date_col] if date_col in grp.columns else grp.index,
            grp[score_col]
        ))
    return build_nps_columns(scores_by_subj, test_name)


# =============================================================================
# §1  DIAGNOSIS TRAJECTORY
# =============================================================================
# Source 1: MCI_Registered — RC_MCI per visit
#   CU  → "CN"    (cognitively unimpaired = normal)
#   MCI → "MCI"
#   RC_MCI_subtype: aMCI_SD, aMCI_MD, naMCI_SD, naMCI_MD
#
# Source 2: Clinical_diagnosis_Registered — free-text per event
#   Used to supplement with external clinical labels when RC_MCI is absent
#   or when a formal dementia diagnosis was received externally.
#
# All PREVENT-AD subjects are CU at enrollment (at-risk adults, 60+, FH+).
# Most will remain CU; a subset progresses to MCI or dementia during follow-up.
# =============================================================================

CLINICAL_DX_MAP = {
    # Common strings found in free-text Clinical_diagnosis
    "cu"               : "CN",
    "cognitively unimpaired": "CN",
    "normal"           : "CN",
    "healthy"          : "CN",
    "mci"              : "MCI",
    "mild cognitive"   : "MCI",
    "amci"             : "MCI",
    "naMCI"            : "MCI",
    "dementia"         : "DEMENTIA",
    "alzheimer"        : "AD",
    "ad"               : "AD",
}


def _classify_rc_mci(val) -> str:
    if pd.isna(val):
        return "Unknown"
    try:
        n = int(float(val))
        if n == 0:
            return "CN"
        if n == 1:
            return "MCI"
    except (ValueError, TypeError):
        pass
    t = str(val).strip().upper()
    if t in ("0", "CU", "COGNITIVELY UNIMPAIRED"):
        return "CN"
    if t in ("1",) or t.startswith("MCI") or "MCI" in t:
        return "MCI"
    return "Unknown"


def _classify_clinical_dx(val) -> str | None:
    if pd.isna(val) or str(val).strip() == "":
        return None
    t = str(val).strip().lower()
    for key, label in CLINICAL_DX_MAP.items():
        if key in t:
            return label
    return str(val).strip()[:30]   # keep raw if unmapped


def build_diagnosis(df_mci: pd.DataFrame,
                    df_cdx: pd.DataFrame) -> pd.DataFrame:
    _print("  §1  Diagnosis trajectory (RC_MCI + Clinical_diagnosis)...")
    records = {}

    # ── Primary: MCI research classification ─────────────────────────────────
    if not df_mci.empty:
        mci = normalise_id(df_mci.copy())
        safe_date(mci, "Date")
        mci["_sort"] = sort_key(mci)
        mci["DX_LABEL"] = mci["RC_MCI"].apply(_classify_rc_mci)
        mci = mci.sort_values(["PTID", "_sort"])

        for ptid, grp in mci.groupby("PTID"):
            labels = grp["DX_LABEL"].tolist()
            dates  = grp.get("Date", grp.get("Date_taken",
                             pd.Series([pd.NaT]*len(grp)))).tolist()
            first_date = dates[0]

            changes_dates, changes_days = [], []
            prev = labels[0]
            for lbl, dt in zip(labels[1:], dates[1:]):
                if lbl != prev:
                    changes_dates.append(dt)
                    try:
                        if pd.notna(first_date) and pd.notna(dt):
                            delta = (pd.to_datetime(dt) -
                                     pd.to_datetime(first_date)).days
                        else:
                            delta = np.nan
                    except Exception:
                        delta = np.nan
                    changes_days.append(delta)
                    prev = lbl

            subtypes = (
                grp["RC_MCI_subtype"]
                .replace({0: np.nan, "0": np.nan})
                .dropna()
                .unique()
                .tolist()
            ) if "RC_MCI_subtype" in grp.columns else []

            records[ptid] = {
                "ID"                         : ptid,
                "1st_diag"                   : labels[0],
                "Lst_diag"                   : labels[-1],
                "1st_diag_date"              : dates[0],
                "Lst_diag_date"              : dates[-1],
                "all_diag"                   : list_str(labels),
                "all_diag_date"              : list_str([d_fmt(d)
                                                          for d in dates]),
                "nb_of_visits"               : len(grp),
                "1st_time_of_change"         : changes_dates[0]  if changes_dates else np.nan,
                "1st_time_of_change_days_to" : changes_days[0]   if changes_days  else np.nan,
                "all_time_of_change"         : list_str([d_fmt(d)
                                                          for d in changes_dates]),
                "all_time_of_change_days_to" : list_str(changes_days),
                "Lst_time_of_change"         : changes_dates[-1] if changes_dates else np.nan,
                "nb_of_change"               : len(changes_dates),
                "Diag_total"                 : list_str(sorted(set(labels))),
                "MCI_subtype_all"            : list_str(subtypes),
            }

    # ── Supplement: clinical diagnosis (external) ─────────────────────────────
    if not df_cdx.empty:
        cdx = normalise_id(df_cdx.copy())
        safe_date(cdx, "Date")
        for ptid, grp in cdx.groupby("PTID"):
            rec = records.setdefault(ptid, {"ID": ptid})
            labels_ext = [_classify_clinical_dx(v)
                          for v in grp["Clinical_diagnosis"]
                          if pd.notna(v)]
            labels_ext = [l for l in labels_ext if l]
            if labels_ext:
                rec["clinical_dx_external"] = list_str(labels_ext)
                rec["clinical_dx_source"]   = list_str(
                    grp["Source_clinical_diagnosis"].dropna().tolist()
                    if "Source_clinical_diagnosis" in grp.columns else []
                )

    result = pd.DataFrame(list(records.values()))
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §2  MRI — inferred from PET_info T1 dates
# =============================================================================
# No raw MRI acquisition table in the PREVENT-AD public registered release.
# PET_info_participants has:
#   Freesurfer_visit_label : visit label of T1 used for PET processing
#   Date_MRI               : date of that T1 (YYYY-MM)
# This confirms T1w MRI was acquired at that visit.
# Synthetic image ID: PTID_PETsession
# =============================================================================

def build_mri(df_pet_info: pd.DataFrame) -> pd.DataFrame:
    _print("  §2  MRI (T1w presence from PET_info T1 reference)...")
    if df_pet_info.empty:
        return pd.DataFrame()

    pet = normalise_id(df_pet_info.copy())
    safe_date(pet, "Date_MRI")

    records = {}
    for ptid, grp in pet.sort_values("Date_MRI").groupby("PTID"):
        dates = grp["Date_MRI"].dropna().tolist()
        visits = grp.get("Freesurfer_visit_label",
                          grp.get("PET_session",
                                   pd.Series([np.nan]*len(grp)))).tolist()
        uniq_dates = list(dict.fromkeys(
            str(d.date()) if hasattr(d, "date") and pd.notna(d) else str(d)
            for d in dates if pd.notna(d)
        ))
        if not uniq_dates:
            continue
        rec = {"ID": ptid, "MRI_type_disponible": "MRI_T1w"}
        # one ID per DATE (aligned lists); every ID of the date in all_MRI_T1w_image_ids
        session_id_columns(rec, "MRI_T1w", grp["Date_MRI"].tolist(),
                           [f"{ptid}_{v}" if pd.notna(v) else None for v in visits])
        rec["1st_MRI_T1w_date"] = uniq_dates[0]
        rec["Lst_MRI_T1w_date"] = uniq_dates[-1]
        records[ptid] = rec

    if not records:
        return pd.DataFrame()
    return pd.DataFrame(list(records.values())).set_index("ID")


# =============================================================================
# §3  PET — 18F-NAV4694 (amyloid)
# =============================================================================
# Source: PET_info_participants_Registered_PREVENTAD_2024-08-01.csv
#   CONP_CandID, PET_session, PET_tracer, Date_PET_visit,
#   Age_PET_visit, Weight, Tracer_dose, MMSE_score
#
# PET_tracer value: "NAV4694" → ABETA_NAV4694
# No regional SUVR in PET_info (only acquisition metadata).
# Regional SUVR in data_dictionary but file not in tree (not available).
#
# Columns produced:
#   PET_tracer_disponible, all_ABETA_NAV4694_date, nb_total_ABETA_NAV4694,
#   ABETA_NAV4694_IDs, 1st/Lst_ABETA_NAV4694_date
#   Also: all_ABETA_NAV4694_age (age at scan), all_ABETA_NAV4694_mmse
# =============================================================================

def build_pet(df_pet: pd.DataFrame) -> pd.DataFrame:
    _print("  §3  PET acquisitions (NAV4694 amyloid)...")
    if df_pet.empty:
        return pd.DataFrame()

    pet = normalise_id(df_pet.copy())
    safe_date(pet, "Date_PET_visit")
    pet = pet.sort_values(["PTID", "Date_PET_visit"])

    records = {}
    TRACER = "ABETA_NAV4694"

    for ptid, grp in pet.groupby("PTID"):
        dates = []
        for _, row in grp.iterrows():
            d = row.get("Date_PET_visit", pd.NaT)
            dates.append(str(d.date()) if hasattr(d, "date") and pd.notna(d)
                         else str(d))

        ids = [f"{ptid}_{row.get('PET_session', i)}"
               for i, (_, row) in enumerate(grp.iterrows())]

        records[ptid] = {
            "ID"                        : ptid,
            "PET_tracer_disponible"     : TRACER,
            f"all_{TRACER}_date"        : list_str(dates),
            f"nb_total_{TRACER}"        : len(dates),
            f"{TRACER}_IDs"             : list_str(ids),
            f"1st_{TRACER}_date"        : dates[0]  if dates else np.nan,
            f"Lst_{TRACER}_date"        : dates[-1] if dates else np.nan,
            f"all_{TRACER}_age"         : list_str(
                grp.get("Age_PET_visit", pd.Series()).dropna().tolist()),
            f"all_{TRACER}_mmse"        : list_str(
                grp.get("MMSE_score", pd.Series()).dropna().tolist()),
            f"all_{TRACER}_dose"        : list_str(
                grp.get("Tracer_dose", pd.Series()).dropna().tolist()),
        }

    if not records:
        return pd.DataFrame()
    return pd.DataFrame(list(records.values())).set_index("ID")


# =============================================================================
# §4  CSF PROTEINS
# =============================================================================
# Source: CSF_proteins_Registered_PREVENTAD_2025-03-21.csv
# Columns: tau, ptau, Amyloid_beta_1_42,
#          Amyloid_beta_1_40_lumipulse, Amyloid_beta_1_42_lumipulse,
#          ptau217_simoa, ApoE (protein level), G_CSF, IL_15, IL_8, VEGF, PCSK9
#
# Maps to ADNI v2 CSF canonical names:
#   tau                     → TAU
#   ptau                    → PTAU
#   Amyloid_beta_1_42       → ABETA42
#   Amyloid_beta_1_40_lumipulse → ABETA40_lumipulse
#   Amyloid_beta_1_42_lumipulse → ABETA42_lumipulse
#   ptau217_simoa           → PTAU217
# =============================================================================

CSF_MARKERS = [
    ("tau",                        "TAU"),
    ("ptau",                       "PTAU"),
    ("Amyloid_beta_1_42",          "ABETA42"),
    ("Amyloid_beta_1_40_lumipulse","ABETA40_lumipulse"),
    ("Amyloid_beta_1_42_lumipulse","ABETA42_lumipulse"),
    ("ptau217_simoa",              "PTAU217_CSF"),
    ("ApoE",                       "CSF_ApoE"),
    ("G_CSF",                      "CSF_GCSF"),
    ("IL_15",                      "CSF_IL15"),
    ("IL_8",                       "CSF_IL8"),
    ("VEGF",                       "CSF_VEGF"),
    ("PCSK9",                      "CSF_PCSK9"),
]

# Laboratory / kit per canonical marker — from the 'notes' column of
# Data_Dictionary_Registered_PREVENTAD_20250321.csv
_UGOT = "Univ. Gothenburg Clinical Neurochemistry Lab"
PREVENTAD_ASSAY_META = {
    # CSF
    "TAU":               (None,  "INNOTEST ELISA (Fujirebio)"),
    "PTAU":              (None,  "INNOTEST ELISA (Fujirebio)"),
    "ABETA42":           (None,  "INNOTEST ELISA (Fujirebio)"),
    "ABETA40_lumipulse": (_UGOT, "Fujirebio Lumipulse G (2024)"),
    "ABETA42_lumipulse": (_UGOT, "Fujirebio Lumipulse G (2024)"),
    "PTAU217_CSF":       (_UGOT, "in-house Simoa p-tau217 (2024)"),
    "CSF_ApoE":          (None,  "MILLIPLEX Luminex (Millipore)"),
    "CSF_GCSF":          (None,  "MILLIPLEX Luminex (Millipore)"),
    "CSF_IL15":          (None,  "MILLIPLEX Luminex (Millipore)"),
    "CSF_IL8":           (None,  "MILLIPLEX Luminex (Millipore)"),
    "CSF_VEGF":          (None,  "MILLIPLEX Luminex (Millipore)"),
    "CSF_PCSK9":         (None,  "Quantikine ELISA (R&D Systems)"),
    # Plasma
    "Plasma_AB40_IPMS":      (_UGOT, "IP-MS (2018)"),
    "Plasma_AB42_IPMS":      (_UGOT, "IP-MS (2018)"),
    "Plasma_AB_ratio_IPMS":  (_UGOT, "IP-MS (2018)"),
    "Plasma_tTau_simoa":     (_UGOT, "in-house Simoa (2018)"),
    "Plasma_pTau181_simoa":  (_UGOT, "in-house Simoa (2021)"),
    "Plasma_pTau231_simoa":  (_UGOT, "in-house Simoa (2021)"),
    "Plasma_AB40_4plex":     (_UGOT, "Quanterix Simoa Neuro 4-Plex E, lot 503105, HD-X (2021)"),
    "Plasma_AB42_4plex":     (_UGOT, "Quanterix Simoa Neuro 4-Plex E, lot 503105, HD-X (2021)"),
    "Plasma_NFL":            (_UGOT, "Quanterix Simoa Neuro 4-Plex E, lot 503105, HD-X (2021)"),
    "Plasma_GFAP":           (_UGOT, "Quanterix Simoa Neuro 4-Plex E, lot 503105, HD-X (2021)"),
    "Plasma_AB_ratio_4plex": (_UGOT, "Quanterix Simoa Neuro 4-Plex E, lot 503105, HD-X (2021)"),
    "Plasma_pTau217_simoa":  (_UGOT, "in-house Simoa p-tau217 (2024)"),
}


def _add_assay_meta(df: pd.DataFrame) -> pd.DataFrame:
    """Add {marker}_lab / {marker}_kit for subjects with >= 1 value."""
    for name, (lab, kit) in PREVENTAD_ASSAY_META.items():
        nb = f"nb_of_{name}_visits"
        if nb not in df.columns:
            continue
        has = df[nb].fillna(0).gt(0)
        if lab:
            df.loc[has, f"{name}_lab"] = lab
        if kit:
            df.loc[has, f"{name}_kit"] = kit
    return df


def build_csf(df_csf: pd.DataFrame) -> pd.DataFrame:
    _print("  §4  CSF biomarkers...")
    if df_csf.empty:
        return pd.DataFrame()

    csf = normalise_id(df_csf.copy())
    safe_date(csf, "Date_taken")
    csf["_sort"] = sort_key(csf)
    csf = csf.sort_values(["PTID", "_sort"])

    parts = []
    for raw_col, canonical in CSF_MARKERS:
        if raw_col not in csf.columns:
            continue
        # '>ULoQ' / '<LLoQ' text → numeric (limit) or NaN; count them
        cens = csf[raw_col].astype(str).str.strip().str.match(r"^[<>]")
        csf[raw_col] = plasma_clean_numeric(csf[raw_col])
        scores_by_subj = {}
        for subj, grp in csf.groupby("PTID"):
            scores_by_subj[subj] = list(zip(
                grp["Date_taken"], grp[raw_col]
            ))
        part = build_nps_columns(scores_by_subj, canonical)
        if not part.empty:
            part[f"nb_of_{canonical}_censored"] = (
                cens.groupby(csf["PTID"]).sum().reindex(part.index).fillna(0).astype(int))
            parts.append(part)

    # Summary flag
    if parts:
        all_ids = set()
        for p in parts:
            all_ids.update(p.index.tolist())
        flag = pd.DataFrame(
            {"CSF_available": "yes"}, index=pd.Index(sorted(all_ids), name="ID")
        )
        parts.insert(0, flag)

    if not parts:
        return pd.DataFrame()
    return _add_assay_meta(reduce(lambda a, b: a.join(b, how="outer"), parts))


# =============================================================================
# §5  PLASMA BIOMARKERS
# =============================================================================
# Plasma_4plex: AB_1_40_MS_UGOT, AB_1_42_MS_UGOT, AB_ratio_MS_UGOT,
#               t_tau_simoa_UGOT, ptau181_simoa_UGOT, ptau231_simoa_UGOT,
#               AB40_simoa_4plex, AB42_simoa_4plex, NFL_simoa_4plex,
#               GFAP_simoa_4plex, AB_ratio_simoa_4plex
# Plasma_p-tau217: ptau217_simoa_UGOT
# =============================================================================

PLASMA_MARKERS_4PLEX = [
    ("AB_1_40_MS_UGOT",    "Plasma_AB40_IPMS"),
    ("AB_1_42_MS_UGOT",    "Plasma_AB42_IPMS"),
    ("AB_ratio_MS_UGOT",   "Plasma_AB_ratio_IPMS"),
    ("t_tau_simoa_UGOT",   "Plasma_tTau_simoa"),
    ("ptau181_simoa_UGOT", "Plasma_pTau181_simoa"),
    ("ptau231_simoa_UGOT", "Plasma_pTau231_simoa"),
    ("AB40_simoa_4plex",   "Plasma_AB40_4plex"),
    ("AB42_simoa_4plex",   "Plasma_AB42_4plex"),
    ("NFL_simoa_4plex",    "Plasma_NFL"),
    ("GFAP_simoa_4plex",   "Plasma_GFAP"),
    ("AB_ratio_simoa_4plex","Plasma_AB_ratio_4plex"),
]


def build_plasma(dfs: dict) -> pd.DataFrame:
    _print("  §5  Plasma biomarkers...")
    parts = []

    # 4plex
    p4 = dfs.get("plasma_4plex", pd.DataFrame())
    if not p4.empty:
        p4 = normalise_id(p4.copy())
        safe_date(p4, "Date_taken")
        p4["_sort"] = sort_key(p4)
        p4 = p4.sort_values(["PTID", "_sort"])
        for raw_col, _ in PLASMA_MARKERS_4PLEX:
            if raw_col in p4.columns:
                p4[raw_col] = plasma_clean_numeric(p4[raw_col])
        for raw_col, canonical in PLASMA_MARKERS_4PLEX:
            part = _nps_from_df(p4, canonical, raw_col)
            if not part.empty:
                parts.append(part)

    # p-tau217
    pt = dfs.get("plasma_ptau217", pd.DataFrame())
    if not pt.empty:
        pt = normalise_id(pt.copy())
        safe_date(pt, "Date_taken")
        pt["_sort"] = sort_key(pt)
        if "ptau217_simoa_UGOT" in pt.columns:
            pt["ptau217_simoa_UGOT"] = plasma_clean_numeric(pt["ptau217_simoa_UGOT"])
        part = _nps_from_df(pt, "Plasma_pTau217_simoa", "ptau217_simoa_UGOT")
        if not part.empty:
            parts.append(part)

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()
    return _add_assay_meta(reduce(lambda a, b: a.join(b, how="outer"), parts))


# =============================================================================
# §6  CLINICAL / NPS / FUNCTIONAL / VITALS / LAB
# =============================================================================
# NPS test mapping (raw_col, canonical, table_key):
#
# EL_CDR_MoCA : CDR_total_score → CDR
#               MOCA_total_score → MoCA
#               delayed_recall_score → MoCA_DelayedRecall
#               orientation_score    → MoCA_Orientation
# CDR_FU      : CDR_total_score → CDR   (longitudinal CDR at follow-up)
# RAVLT       : ListA_Sum_Trials_1_5 → RAVLT_I (total learning)
#               ListA_Delayed_Recall_Correct_Answer → RAVLT_D
#               List_A_Recognition_Recognized → RAVLT_Recog
# RBANS       : total_scale_index_score → RBANS_total
#               immediate_memory_index_score → RBANS_ImmMem
#               delayed_memory_index_score → RBANS_DelMem
#               attention_index_score → RBANS_Attn
#               digit_span_score → DS_RBANS
#               semantic_fluency_score → VF_anim_RBANS
#               coding_score → DigSym_RBANS
# TMT         : TMT_A_time → TMT_A
#               TMT_B_time → TMT_B
# DKEFS       : Color_Naming_Time → Stroop_CN
#               Inhibition_Time → Stroop_Inhib
#               Inhibition_Scaled_Score → Stroop_Inhib_scaled
# AD8         : AD8_total_score → AD8
# APS         : APS_score → APS
# Smell       : UPSIT_total_score → UPSIT
# BP/Vitals   : Systolic_blood_pressure → BP_Sys
#               Diastolic_blood_pressure → BP_Dias
#               Pulse → Pulse, Weight → Weight
# Lab         : hba1c_value → HbA1c
#               total_cholesterol_value → Cholesterol
#               HDL_value → HDL, LDL_value → LDL
# EL_CAIDE    : CAIDE_total_score → CAIDE
# =============================================================================




def build_all_clinical(dfs: dict) -> pd.DataFrame:
    """Clinical / NPS / vitals / lab series.

    Variables come from ``variables.yaml`` (``sources: PREVENTAD:``); rows are
    ordered by candidate age and dated with ``Date_taken``.
    """
    _print("  §6  Clinical / NPS / Vitals / Lab...")
    parts = []

    for name, src in variable_sources(DATASET, config=_CTX.config):
        df = dfs.get(src.table, pd.DataFrame())
        if df.empty:
            continue
        frame, col = materialise(df, name, src)
        if col is None:
            continue
        df = normalise_id(frame.copy())
        safe_date(df, "Date_taken")
        df["_sort"] = sort_key(df)
        df = df.sort_values(["PTID", "_sort"])
        part = _nps_from_df(df, name, col, src.date_column or "Date_taken")
        if not part.empty:
            parts.append(part)

    # BMI: Weight (kg, every visit) / Height (cm, eligibility visit)²
    bw  = dfs.get("bp_weight", pd.DataFrame())
    dem = dfs.get("demographics", pd.DataFrame())
    if not bw.empty and not dem.empty and "Weight" in bw.columns and "Height" in dem.columns:
        bw  = normalise_id(bw.copy())
        dem = normalise_id(dem.copy()).drop_duplicates("PTID")
        bw["_height"] = bw["PTID"].map(dem.set_index("PTID")["Height"])
        bw["_BMI"] = compute_bmi(bw["Weight"], bw["_height"], "kg", "cm")
        safe_date(bw, "Date_taken")
        bw["_sort"] = sort_key(bw)
        bw = bw.sort_values(["PTID", "_sort"])
        part = _nps_from_df(bw, "BMI", "_BMI")
        if not part.empty:
            parts.append(part)

    # Antiepileptic text flag from medication lists (added 2026-11-08).
    # PREVENT-AD has no seizure field; this only catches AED names if the
    # SU_/PRN_medication columns hold drug names (not counts) — check AED_text_n.
    mu = dfs.get("med_use", pd.DataFrame())
    if not mu.empty:
        mu = normalise_id(mu.copy())
        safe_date(mu, "Date_taken")
        part = text_flag_series(mu, "PTID", ["SU_medication", "PRN_medication"],
                                AED_TEXT_RX, "AED_text", "Date_taken")
        if not part.empty:
            parts.append(part)

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()
    return reduce(lambda a, b: a.join(b, how="outer"), parts)


# =============================================================================
# §7  DEMOGRAPHICS + GENETICS (APOE)
# =============================================================================
# Demographics_Registered (static per subject — no visit structure):
#   Sex, Ethnicity, Education_years, Education_level,
#   Work_by_category, Retirement_status, Height,
#   Handedness_result, Handedness_interpretation,
#   INTREPAD_Tx_assignment (treatment assignment: Placebo / Active),
#   Family history: father/mother/sibling dx AD dementia flags
#
# Genetics_Registered (static):
#   APOE → genotype string (e.g. "e3/e4")
#   APOE_e4_count = count of "4" or "e4" in APOE string
#   Also: BchE_K_variant, BDNF, HMGCR_Intron_M (secondary genetics)
# =============================================================================

def build_demographics(dfs: dict) -> pd.DataFrame:
    _print("  §7  Demographics + Genetics...")
    records: dict = {}

    # ── Demographics (static — no visit) ──────────────────────────────────────
    dem = dfs.get("demographics", pd.DataFrame())
    if not dem.empty:
        dem = normalise_id(dem.copy())
        # One row per subject (static)
        dem = dem.drop_duplicates("PTID")
        for _, row in dem.iterrows():
            ptid = row["PTID"]
            records[ptid] = {
                "ID"                    : ptid,
                "gender"                : row.get("Sex",              np.nan),
                "ethnicity"             : row.get("Ethnicity",        np.nan),
                "education_years"       : row.get("Education_years",  np.nan),
                "education_level"       : row.get("Education_level",  np.nan),
                "retirement_status"     : row.get("Retirement_status",np.nan),
                "handedness"            : row.get("Handedness_interpretation", np.nan),
                "height"                : row.get("Height",           np.nan),
                "tx_assignment"         : row.get("INTREPAD_Tx_assignment", np.nan),
                "mother_tongue"         : row.get("Mother_tongue",    np.nan),
                "famhx_father_AD"       : row.get("father_dx_ad_dementia", np.nan),
                "famhx_mother_AD"       : row.get("mother_dx_ad_dementia", np.nan),
                "famhx_sibling_AD"      : row.get("sibling_dx_ad_dementia", np.nan),
                "famhx_other_AD"        : row.get("other_family_members_AD", np.nan),
            }

    # ── Genetics ──────────────────────────────────────────────────────────────
    gen = dfs.get("genetics", pd.DataFrame())
    if not gen.empty:
        gen = normalise_id(gen.copy())
        gen = gen.drop_duplicates("PTID")
        for _, row in gen.iterrows():
            ptid = row["PTID"]
            rec  = records.setdefault(ptid, {"ID": ptid})
            geno = str(row.get("APOE", "") or "").strip()
            rec["APOE_genotype"] = geno if geno and geno.lower() != "nan" else np.nan
            # Count e4 alleles — PREVENT-AD may use "e3/e4" or "3/4" format
            rec["APOE_e4_count"] = geno.lower().count("e4") or geno.count("4") \
                                   if geno else np.nan
            rec["BchE_K_variant"]   = row.get("BchE_K_variant",    np.nan)
            rec["BDNF_Val66Met"]    = row.get("BDNF",              np.nan)

    df = pd.DataFrame(list(records.values()))
    return df.set_index("ID") if not df.empty else df


# =============================================================================
# §8  AGE AT FIRST DIAGNOSIS
# =============================================================================
# Candidate_Age is in MONTHS. Convert to years for cross-cohort consistency.
# age_at_1st_diag = Candidate_Age (months) / 12 at first MCI visit.

def add_age_at_1st_diag(master: pd.DataFrame,
                         dfs: dict) -> pd.DataFrame:
    mci = dfs.get("mci", pd.DataFrame())
    if mci.empty:
        return master

    mci = normalise_id(mci.copy())
    mci["_sort"] = sort_key(mci)
    first = mci.sort_values("_sort").drop_duplicates("PTID", keep="first")
    age_map = (first.set_index("PTID")["Candidate_Age"]
                    .dropna() / 12).round(2).to_dict()

    master["age_at_1st_diag"] = master["ID"].map(age_map)
    return master


# =============================================================================
# MAIN
# =============================================================================

def main() -> pd.DataFrame:
    _print("=" * 62)
    _print("PREVENT-AD LONGITUDINAL INTEGRATION PIPELINE  v1")
    _print("=" * 62)
    _print("\nSubject ID : CONP_CandID → PTID")
    _print("Time axis  : Candidate_Age (months) → sorted chronologically")
    _print("Diagnosis  : RC_MCI (research) + Clinical_diagnosis (external)")
    _print("PET tracer : 18F-NAV4694 (amyloid)")
    _print("CSF        : tau, ptau, Aβ42, Aβ40 (Lumipulse), ptau217")
    _print("Plasma     : AB40/42, tTau, pTau181/231/217, NFL, GFAP")

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load ──────────────────────────────────────────────────────────────────
    _print("\nLoading raw tables...")
    dfs = {k: load(k) for k in PATHS}

    # ── Build domains ─────────────────────────────────────────────────────────
    _print("\nBuilding domain tables...")
    diag   = build_diagnosis(dfs["mci"], dfs["clinical_dx"])
    mri    = build_mri(dfs["pet_info"])
    pet    = build_pet(dfs["pet_info"])
    csf    = build_csf(dfs["csf"])
    plasma = plasma_add_summary(build_plasma(dfs))
    clin   = build_all_clinical(dfs)
    demog  = build_demographics(dfs)

    # ── Merge ─────────────────────────────────────────────────────────────────
    _print("\nMerging all domains on subject ID...")
    parts = [p for p in [diag, mri, pet, csf, plasma, clin, demog]
             if isinstance(p, pd.DataFrame) and not p.empty]

    if not parts:
        _print("ERROR: no domain tables produced.")
        return pd.DataFrame()

    master = reduce(lambda a, b: a.join(b, how="outer"), parts)
    master.index.name = "ID"
    master = master.reset_index()
    gc.collect()

    # ── Post-processing ───────────────────────────────────────────────────────
    master = add_age_at_1st_diag(master, dfs)

    if "MRI_type_disponible" not in master.columns:
        master["MRI_type_disponible"] = np.nan
    if "PET_tracer_disponible" not in master.columns:
        master["PET_tracer_disponible"] = np.nan
    if "CSF_available" not in master.columns:
        master["CSF_available"] = np.nan

    master["Plasma_available"] = (
        master["Plasma_panels_disponible"].notna().map({True: "yes", False: "no"})
        if "Plasma_panels_disponible" in master.columns else "no")

    LEAD = [
        "ID",
        "gender", "ethnicity", "education_years", "education_level",
        "retirement_status", "handedness", "height", "mother_tongue",
        "tx_assignment",
        "APOE_genotype", "APOE_e4_count", "BchE_K_variant", "BDNF_Val66Met",
        "age_at_1st_diag",
        "famhx_father_AD", "famhx_mother_AD", "famhx_sibling_AD",
        "famhx_other_AD",
        "1st_diag", "Lst_diag", "1st_diag_date", "Lst_diag_date",
        "all_diag", "all_diag_date", "nb_of_visits",
        "1st_time_of_change", "1st_time_of_change_days_to",
        "all_time_of_change", "all_time_of_change_days_to",
        "Lst_time_of_change", "nb_of_change", "Diag_total",
        "MCI_subtype_all", "clinical_dx_external",
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
              f"{master['nb_of_change'].gt(0).sum()}")

    if "PET_tracer_disponible" in master.columns:
        n_pet = master["PET_tracer_disponible"].notna().sum()
        _print(f"Subjects with NAV4694 PET: {n_pet:,}")

    if "CSF_available" in master.columns:
        n_csf = master["CSF_available"].eq("yes").sum()
        _print(f"Subjects with CSF data: {n_csf:,}")

    _print("\nTest coverage (subjects with ≥1 valid score):")
    for name in list(dict.fromkeys(n for n, _ in variable_sources(DATASET, config=_CTX.config)))[:12]:
        col = f"nb_of_{name}_visits"
        if col in master.columns:
            n = master[col].gt(0).sum()
            if n > 0:
                _print(f"  {name:<28s} {n:>5,} subjects")

    return master


# =============================================================================
