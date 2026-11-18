"""ADNI longitudinal integration pipeline (v2).

One row per subject summarising the full longitudinal record across:
  §1  Diagnosis trajectory   (DXSUM)
  §2  MRI acquisitions       (All_Subjects_Key_MRI)
  §3  PET acquisitions       (All_Subjects_PET_Images)
  §4  Neuropsychological     (NEUROBAT + CDR + GDS + MMSE + MoCA + FAQ + NPI + STAI + ECOG)
  §5  CSF biomarkers         (UPENNBIOMK_ROCHE_ELECSYS)
  §5b Plasma biomarkers      (UPENN Fujirebio/Quanterix, C2N, Lilly, FNIH, UGOT, Blennow ...)
  §5c Routine blood lab      (LABDATA)
  §6  Demographics + APOE    (PTDEMOG + Study_Entry + APOERES)

Plasma, extra CSF, lab and medical-history tables are found by their ADNI
table name anywhere under ``Tables_ADNI/`` (settings ``plasma_specs``,
``csf_extra_specs``, ``lab_specs``, ``medhist_tables`` of
``configs/datasets/adni.yaml``); most of them are keyed on RID only, so RID is
mapped to PTID through DXSUM / PTDEMOG.

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("ADNI", root_path)
"""
from __future__ import annotations

import gc
import numpy as np
import os
import pandas as pd
import re
from collections import Counter
from functools import reduce

from ..clinical import AED_TEXT_RX, SEIZURE_TEXT_RX, text_flag_series
from ..datasets.base import DatasetContext
from ..diagnosis import _get, add_dx_columns, trajectory_record
from ..fluid import (
    plasma_add_summary,
    plasma_clean_numeric,
    plasma_find_tables,
    plasma_from_spec,
    plasma_join,
    plasma_read_csv,
    plasma_series,
    specs_from_config,
)
from ..imaging.scanner import add_scanner_channels
from ..imaging.sessions import session_id_columns
from ..series import build_scored_columns, days_from, list_str, safe_date
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

DATASET = "ADNI"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
ADNI = ADNI_EXPLORE = ADNI_NPS = ADNI_PLASMA_DIR = OUT_DIR = OUT_LONG = ""
ADNI_PLASMA_SPECS: dict = {}
ADNI_LAB_SPEC: dict = {}
ADNI_CSF_EXTRA_SPECS: dict = {}
ADNI_MEDHIST_TABLES: list = []
CSF_MAIN_LAB = CSF_MAIN_KIT = ""


def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns (ADNI keeps the text 'NaT' for a missing first/last date)."""
    return build_scored_columns(scores_by_subj, test_name, missing_date_as_text=True)


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, ADNI, ADNI_EXPLORE, ADNI_NPS, ADNI_PLASMA_DIR, OUT_DIR, OUT_LONG
    global ADNI_PLASMA_SPECS, ADNI_LAB_SPEC, ADNI_CSF_EXTRA_SPECS, ADNI_MEDHIST_TABLES
    global CSF_MAIN_LAB, CSF_MAIN_KIT
    _CTX = ctx
    ADNI = ADNI_EXPLORE = ADNI_NPS = ADNI_PLASMA_DIR = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()
    s = ctx.settings
    ADNI_PLASMA_SPECS = specs_from_config(s.get("plasma_specs", {}))
    ADNI_LAB_SPEC = specs_from_config(s.get("lab_specs", {}))
    ADNI_CSF_EXTRA_SPECS = specs_from_config(s.get("csf_extra_specs", {}))
    ADNI_MEDHIST_TABLES = list(s.get("medhist_tables", []))
    CSF_MAIN_LAB = s.get("csf_main_lab", "")
    CSF_MAIN_KIT = s.get("csf_main_kit", "")


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write ``ADNI_longitudinal.csv``; returns the per-subject table."""
    configure(ctx)
    return main()

# =============================================================================
# ── TYPE / TRACER MAPS ────────────────────────────────────────────────────────
# =============================================================================

MRI_TYPE_MAP = {
    # Standard BIDS-like labels used by ADNI Key_MRI table
    "T1w"        : "T1w",
    "T2w"        : "T2w",
    "FLAIR"      : "T2w_FLAIR",
    "BOLD"       : "fMRI",
    "bold-rest"  : "fMRI",
    "bold_rest"  : "fMRI",
    "dMRI"       : "DWI",
    "DWI"        : "DWI",
    "dwi"        : "DWI",
    "ASL"        : "ASL",
    "asl"        : "ASL",
    "T2starw"    : "T2star",
    "T2star"     : "T2star",
    "SWI"        : "T2star",
    "Swi"        : "T2star",
    "PDT2"       : "pT2",
    "pT2"        : "pT2",
    "MET2starw"  : "other_T2",
}

PET_TYPE_MAP = {
    # ── Tau tracers ───────────────────────────────────────────────────────────
    "18F-AV1451"        : "TAU_Flortaucipir",
    "18F-FLORTAUCIPIR"  : "TAU_Flortaucipir",
    "AV-1451"           : "TAU_Flortaucipir",
    "FLORTAUCIPIR"      : "TAU_Flortaucipir",
    "18F-MK6240"        : "TAU_MK6240",
    "MK-6240"           : "TAU_MK6240",
    "18F-PI2620"        : "TAU_PI2620",
    "PI-2620"           : "TAU_PI2620",
    "PI2620"            : "TAU_PI2620",
    # ── Amyloid tracers ───────────────────────────────────────────────────────
    "18F-AV45"          : "ABETA_AV45",
    "AV-45"             : "ABETA_AV45",
    "18F-FLORBETAPIR"   : "ABETA_AV45",
    "FLORBETAPIR"       : "ABETA_AV45",
    "18F-FBB"           : "ABETA_FBB",
    "18F-FLORBETABEN"   : "ABETA_FBB",
    "FLORBETABEN"       : "ABETA_FBB",
    "11C-PIB"           : "ABETA_PIB",
    "C11-PIB"           : "ABETA_PIB",
    "PIB"               : "ABETA_PIB",
    "18F-NAV4694"       : "ABETA_NAV4694",
    "NAV4694"           : "ABETA_NAV4694",
    "18F-FLUTEMETAMOL"  : "ABETA_Flutemetamol",
    "FLUTEMETAMOL"      : "ABETA_Flutemetamol",
    # ── Metabolism ────────────────────────────────────────────────────────────
    "18F-FDG"           : "Metabolism",
    "FDG"               : "Metabolism",
}

CSF_MARKERS = ["ABETA40", "ABETA42", "TAU", "PTAU"]

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================




def resolve_date_col(df: pd.DataFrame) -> str | None:
    """Return the first usable date column from a priority list."""
    for c in ("EXAMDATE", "VISDATE", "COLDATE", "pet_date", "image_date",
              "APTESTDT", "_date"):
        if c in df.columns:
            return c
    return None


def combine_dates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a unified '_date' column: EXAMDATE if present, else VISDATE.
    Parses both before combining.
    """
    safe_date(df, "EXAMDATE")
    safe_date(df, "VISDATE")
    if "EXAMDATE" in df.columns and "VISDATE" in df.columns:
        df["_date"] = df["EXAMDATE"].combine_first(df["VISDATE"])
    elif "EXAMDATE" in df.columns:
        df["_date"] = df["EXAMDATE"]
    elif "VISDATE" in df.columns:
        df["_date"] = df["VISDATE"]
    return df


def flag_val(series: pd.Series, fn) -> pd.Series:
    """
    Binary flag: 1 if fn(x) is True, 0 if False, -99 if NaN.
    Used by the diagnosis classifier.
    """
    return series.apply(
        lambda x: -99 if pd.isna(x) else (1 if fn(x) else 0)
    )








def load(key: str) -> pd.DataFrame:
    """Load one CSV from PATHS; print a warning and return empty DF if missing."""
    path = PATHS.get(key, "")
    if not path or not os.path.exists(path):
        _print(f"  [SKIP] {key}: not found → {path}")
        return pd.DataFrame()
    df = pd.read_csv(path, low_memory=False)
    _print(f"  [OK]   {key}: {len(df):>7,} rows  {df.shape[1]} cols")
    return df


def resolve_id(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure the subject-ID column is called 'PTID'.
    Case-insensitive: works before or after column lowercasing.
    """
    if "PTID" in df.columns:
        return df
    # Build a case-insensitive lookup of current column names
    col_lower = {c.lower(): c for c in df.columns}
    for alt in ("ptid", "subject_id", "rid", "subid"):
        if alt in col_lower:
            df = df.rename(columns={col_lower[alt]: "PTID"})
            return df
    return df


# =============================================================================
# §1  DIAGNOSIS TRAJECTORY
# =============================================================================
# ADNI stores diagnosis across two parallel coding schemes:
#
#   Modern (ADNI-GO/2/3/4):
#     DIAGNOSIS  integer  1=CN  2=MCI  3=Dementia
#     Dementia sub-typing uses auxiliary flags:
#       DXAD=1   → probable AD (also DXAPP=1 for possible AD)
#       DXODES   → other dementia etiology code
#                  1,4,5 → FTD   ;  9 → VaD
#       DXPCOG=3 → LBD
#
#   Legacy (ADNI-1):
#     DXNORM=1   → CN
#     DXMCI=1    → MCI
#     DXAD=1     → AD    (used without DIAGNOSIS column)
#
#   DX_SUMMARY (compound label):
#     '+' joins every sub-diagnosis that is simultaneously active,
#     e.g. "MCI+AD" when DIAGNOSIS==3 and DXAD==1 while DXMCI was
#     still flagged in the same visit (rare but documented).
#
# =============================================================================

def _classify_visit(row: pd.Series) -> str:
    """
    Map one DXSUM row to a human-readable label.
    Priority:  modern DIAGNOSIS integer  →  legacy flags  →  'Unknown'
    """
    d = row.get("DIAGNOSIS", np.nan)

    # ── Modern coding ────────────────────────────────────────────────────────
    if pd.notna(d):
        d = int(d)
        if d == 1:
            return "CN"
        if d == 2:
            return "MCI"
        if d == 3:
            parts = []
            if row.get("DXAD", 0) == 1 or row.get("DXAPP", 0) == 1:
                parts.append("AD")
            odes = row.get("DXODES", np.nan)
            if pd.notna(odes):
                odes = int(odes)
                if odes in (1, 4, 5):
                    parts.append("FTD")
                elif odes == 9:
                    parts.append("VaD")
            if row.get("DXPCOG", np.nan) == 3:
                parts.append("LBD")
            return "+".join(parts) if parts else "Dementia-NOS"

    # ── Legacy coding (ADNI-1) ───────────────────────────────────────────────
    # Collect ALL active flags → compound label (same logic as AIBL)
    # DXNORM is mutually exclusive with impairment by design; if DXNORM=1
    # and impairment flags are also 1, DXNORM takes sole priority (data error).
    if row.get("DXNORM", 0) == 1:
        return "CN"

    legacy_parts = []
    if row.get("DXMCI", 0) == 1:
        legacy_parts.append("MCI")
    if row.get("DXAD", 0) == 1 or row.get("DXAPP", 0) == 1:
        legacy_parts.append("AD")
    # DXPARK and DXOTHDEM may also appear in legacy ADNI-1 exports
    if row.get("DXPARK", 0) == 1:
        legacy_parts.append("PD")
    if row.get("DXOTHDEM", 0) == 1:
        legacy_parts.append("OtherDem")

    if legacy_parts:
        return "+".join(sorted(legacy_parts))

    return "Unknown"


ADNI_DXODES_ETIOLOGY = {1: "FTD", 4: "FTD", 5: "FTD", 9: "VaD"}
ADNI_DXCURREN_STAGE = {1: "CN", 2: "MCI", 3: "DEM"}


def _stage_adni(row: pd.Series) -> tuple:
    """(stage, etiology, epilepsy) of one DXSUM row (v0.4.0).

    Stage: DIAGNOSIS (1 CN, 2 MCI, 3 dementia); ADNI-1 rows: DXCURREN, then the
    legacy flags. Etiology: for dementia DXAD/DXAPP -> AD, DXODES -> FTD/VaD,
    DXPCOG = 3 -> LBD; for MCI, DXMDUE (1 = MCI due to AD, 2 = other etiology).
    """
    d = _get(row, "DIAGNOSIS")
    if np.isnan(d):
        d = _get(row, "DXCURREN")
    stage = ADNI_DXCURREN_STAGE.get(int(d)) if not np.isnan(d) else None
    if stage is None:
        if _get(row, "DXNORM") == 1:
            stage = "CN"
        elif _get(row, "DXAD") == 1 or _get(row, "DXAPP") == 1 or _get(row, "DXOTHDEM") == 1:
            stage = "DEM"
        elif _get(row, "DXMCI") == 1:
            stage = "MCI"
        elif _get(row, "DXPARK") == 1:
            return "Unknown", "PD", False
        else:
            return "Unknown", "unknown", False
    if stage == "CN":
        return stage, "none", False
    if stage == "MCI":
        due = _get(row, "DXMDUE")
        return stage, {1: "AD", 2: "Other"}.get(int(due), "unknown") if not np.isnan(due) else "unknown", False
    parts = []
    if _get(row, "DXAD") == 1 or _get(row, "DXAPP") == 1:
        parts.append("AD")
    odes = _get(row, "DXODES")
    if not np.isnan(odes) and int(odes) in ADNI_DXODES_ETIOLOGY:
        parts.append(ADNI_DXODES_ETIOLOGY[int(odes)])
    if _get(row, "DXPCOG") == 3:
        parts.append("LBD")
    if _get(row, "DXOTHDEM") == 1 and not parts:
        parts.append("Other")
    return stage, "+".join(sorted(set(parts))) if parts else "unknown", False


def build_diagnosis(df_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Returns one row per subject (indexed by PTID / ID) with columns:
      1st_diag, Lst_diag, 1st_diag_date, Lst_diag_date,
      all_diag, all_diag_date,
      nb_of_visits,
      1st_time_of_change, 1st_time_of_change_days_to,
      all_time_of_change, all_time_of_change_days_to,
      Lst_time_of_change, nb_of_change,
      Diag_total
    """
    _print("  §1  Diagnosis trajectory...")
    dx = df_raw.copy()
    safe_date(dx, "EXAMDATE")
    dx["DX_LABEL"] = dx.apply(_classify_visit, axis=1)       # pre-0.4.0 label -> *_diag_legacy
    dx = add_dx_columns(dx, _stage_adni)

    records = []
    for ptid, grp in dx.sort_values("EXAMDATE").groupby("PTID"):
        grp = grp.dropna(subset=["EXAMDATE"]).sort_values("EXAMDATE")
        if grp.empty:
            continue

        dates = grp["EXAMDATE"].tolist()
        records.append(trajectory_record(ptid, grp, dates, fmt=lambda d: str(d.date())))

    df = pd.DataFrame(records)
    if df.empty:
        return df
    return df.set_index("ID")


# =============================================================================
# §2  MRI ACQUISITIONS
# =============================================================================
# Source table: All_Subjects_Key_MRI_01Feb2026.csv
# Key columns:
#   subject_id  — links to PTID
#   image_date  — acquisition date
#   series_type — raw BIDS-like label (T1w, T2w, FLAIR, BOLD, dMRI …)
#   image_id    — LONI image ID
#
# Mapping: series_type → canonical type via MRI_TYPE_MAP
# One canonical set of columns is produced per unique mapped type:
#   all_{TYPE}_date, nb_total_{TYPE}, {TYPE}_IDs, 1st_{TYPE}_ID, Lst_{TYPE}_ID
# Plus: MRI_type_disponible  (pipe-joined sorted list of types present)
# =============================================================================

def build_acquisition_columns(df: pd.DataFrame,
                               subj_col: str,
                               date_col: str,
                               id_col:   str | list,
                               type_col: str,
                               type_map: dict,
                               prefix:   str) -> pd.DataFrame:
    """
    Generic acquisition builder shared by MRI and PET sections.

    Parameters
    ----------
    id_col : str or list of str
        Column(s) forming the image identifier.  If a list, values are
        joined with '|' to produce a composite ID.
    """
    df = df.copy()
    # Composite image ID
    if isinstance(id_col, (list, tuple)):
        id_cols = [c for c in id_col if c in df.columns]
        df["_img_id"] = (df[id_cols].astype(str).agg("|".join, axis=1)
                         if id_cols else pd.NA)
    else:
        df["_img_id"] = df[id_col] if id_col in df.columns else pd.NA

    # Case-insensitive type mapping: build lowercase lookup of the map keys
    type_map_ci = {k.lower(): v for k, v in type_map.items()}
    df["_type"] = df[type_col].apply(
        lambda x: type_map_ci.get(str(x).strip().lower(), np.nan)
        if pd.notna(x) else np.nan
    )

    n_before = len(df)
    df = df.dropna(subset=["_type"])
    n_after  = len(df)
    if n_after == 0:
        unique_vals = df[type_col].dropna().unique()[:10] if n_before > 0 else []
        _print(f"    [WARN] {prefix}: all rows dropped after type mapping. "
              f"Sample '{type_col}' values: {list(unique_vals)}")
        _print(f"    [WARN] Known map keys: {list(type_map.keys())[:10]}")

    records, avail_key = {}, f"{prefix}_type_disponible"

    if subj_col not in df.columns:
        _print(f"  [WARN] Subject column '{subj_col}' not found in {prefix} table. "
              f"Available: {list(df.columns[:10])}")
        return pd.DataFrame()

    for subj, grp in df.sort_values(date_col).groupby(subj_col):
        row = records.setdefault(subj, {"ID": subj, avail_key: set()})

        for mtype, tgrp in grp.groupby("_type"):
            tgrp = tgrp.sort_values(date_col)
            row[avail_key].add(mtype)

            # one ID per DATE (aligned lists); every ID of the date in all_{T}_image_ids
            session_id_columns(row, mtype, tgrp[date_col].tolist(),
                               tgrp["_img_id"].tolist())

    for r in records.values():
        r[avail_key] = "|".join(sorted(r[avail_key]))

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(list(records.values()))
    if "ID" not in result.columns:
        return pd.DataFrame()
    result = result.set_index("ID")
    # scanner manufacturer / model per acquisition (auto-detected columns,
    # e.g. mri_mfr / mri_mfr_model in Key_MRI)
    mfr_c = next((c for c in df.columns if re.search(r"mfr$|manufact", c, re.I)), None)
    mod_c = next((c for c in df.columns if re.search(r"mfr_?model|model", c, re.I)), None)
    return add_scanner_channels(result, df, subj_col, "_type", date_col, mfr_c, mod_c)


def build_mri(df_raw: pd.DataFrame) -> pd.DataFrame:
    _print("  §2  MRI acquisitions...")
    if df_raw.empty:
        return pd.DataFrame()
    df = df_raw.copy()
    # Lowercase all columns so lookups are case-insensitive
    df.columns = [c.lower() for c in df.columns]
    # After lowercasing, resolve subject ID to a stable name
    df = resolve_id(df)                    # subject_id / ptid → PTID

    # Detect date column (image_date standard; fallback to any *date* column)
    date_col = "image_date"
    if date_col not in df.columns:
        candidates = [c for c in df.columns if "date" in c]
        if candidates:
            date_col = candidates[0]
            _print(f"    [INFO] Using date column: '{date_col}'")
        else:
            _print("    [WARN] No date column found in MRI file — skipping")
            return pd.DataFrame()
    safe_date(df, date_col)

    # Detect series-type column
    type_col = "series_type"
    if type_col not in df.columns:
        for alt in ("scan_category", "scan category", "seriestype", "modality"):
            if alt in df.columns:
                type_col = alt
                _print(f"    [INFO] Using type column: '{type_col}'")
                break

    # Detect image-ID column
    id_col = "image_id"
    if id_col not in df.columns:
        for alt in ("image_id", "imageid", "label", "loni_image"):
            if alt in df.columns:
                id_col = alt
                break

    return build_acquisition_columns(
        df       = df,
        subj_col = "PTID",
        date_col = date_col,
        id_col   = id_col,
        type_col = type_col,
        type_map = MRI_TYPE_MAP,   # case-insensitive matching done inside builder
        prefix   = "MRI",
    )


# =============================================================================
# §3  PET ACQUISITIONS
# =============================================================================
# Source table: All_Subjects_PET_Images_22Mar2026.csv
# Key columns:
#   subject_id     — links to PTID
#   pet_date       — acquisition date
#   pet_radiopharm — radiotracer name (raw, varied formatting)
#   image_id       — LONI image ID
#   pet_mfr        — scanner manufacturer
#   pet_isotope    — isotope (18F, 11C …)
#
# Tracer normalisation: pet_radiopharm → canonical label via PET_TYPE_MAP.
# All raw strings are upper-cased before lookup; substring matching handles
# vendor-specific naming differences.
# =============================================================================

def _norm_tracer(raw) -> str | None:
    """
    Normalise a raw radiotracer string to a canonical label.
    Returns None for null / empty / uninformative values so those rows
    are excluded from the acquisition summary.
    """
    if pd.isna(raw):
        return None
    t = str(raw).strip()
    if t.upper() in ("", "NAN", "NONE", "N/A", "NA", "UNKNOWN"):
        return None
    t_up = t.upper()
    for key, val in PET_TYPE_MAP.items():
        if key.upper() in t_up:
            return val
    # Unknown-but-real tracer: sanitise and keep
    clean = t_up.replace(" ", "_").replace("-", "_").replace("/", "_")
    return clean[:24] if clean else None


def build_pet(df_raw: pd.DataFrame) -> pd.DataFrame:
    _print("  §3  PET acquisitions...")
    if df_raw.empty:
        return pd.DataFrame()
    df = df_raw.copy()
    df.columns = [c.lower() for c in df.columns]
    df = resolve_id(df)                    # subject_id / ptid → PTID
    safe_date(df, "pet_date")

    # Detect tracer column
    tracer_col = "pet_radiopharm"
    if tracer_col not in df.columns:
        for alt in ("pet_radiopharm", "tracer", "tracer.name", "tracername",
                    "radiopharmaceutical"):
            if alt in df.columns:
                tracer_col = alt
                break

    df["_tracer"] = df[tracer_col].apply(_norm_tracer) if tracer_col in df.columns else "Unknown"

    # Drop rows where tracer could not be mapped (None stays as NaN after apply)
    df = df.dropna(subset=["_tracer"])

    records = {}
    for subj, grp in df.sort_values("pet_date").groupby("PTID"):
        row = records.setdefault(subj, {"ID": subj, "PET_tracer_disponible": set()})

        for tracer, tgrp in grp.groupby("_tracer"):
            tgrp = tgrp.sort_values("pet_date")
            row["PET_tracer_disponible"].add(tracer)

            # one ID per DATE (aligned lists); every ID of the date in all_{T}_image_ids
            session_id_columns(
                row, tracer, tgrp["pet_date"].tolist(),
                tgrp["image_id"].tolist() if "image_id" in tgrp.columns else [None] * len(tgrp))

            if "pet_mfr" in tgrp.columns:
                row[f"{tracer}_manufacturer"] = list_str(tgrp["pet_mfr"].dropna().unique())
            if "pet_isotope" in tgrp.columns:
                row[f"{tracer}_isotope"] = list_str(tgrp["pet_isotope"].dropna().unique())

    for r in records.values():
        r["PET_tracer_disponible"] = "|".join(sorted(r["PET_tracer_disponible"]))

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(list(records.values()))
    if "ID" not in result.columns:
        return pd.DataFrame()
    result = result.set_index("ID")
    mod_c = next((c for c in df.columns if re.search(r"mfr_?model|model", c, re.I)), None)
    return add_scanner_channels(result, df, "PTID", "_tracer", "pet_date",
                                "pet_mfr" if "pet_mfr" in df.columns else None, mod_c)


# =============================================================================
# §4  NEUROPSYCHOLOGICAL + FUNCTIONAL + BEHAVIORAL + VITALS
# =============================================================================
# Shared builder: build_nps_columns
#   Input:  {PTID: [(date, value), ...]}  sorted chronologically
#   Output: DataFrame indexed by ID with per-test columns:
#     all_{test}_date, all_{test}_score,
#     1st_{test}_score, 1st_{test}_date,
#     Lst_{test}_score, Lst_{test}_date,
#     nb_of_{test}_visits
#
# NPS_TESTS list drives which column to pull from which source file.
# The NEUROBAT table is the primary neuropsychological battery source for ADNI;
# it contains a superset of test columns collected at every visit.
#
# Functional / behavioral scales (FAQ, NPI-Q, GDS, STAI, ECOG) are separate
# tables with their own date columns (VISDATE or COLDATE).
# =============================================================================




def _nps_from_source(src: pd.DataFrame, ptid_col: str, date_col: str,
                     test_name: str, score_col: str) -> pd.DataFrame:
    """Extract one test from an already-loaded DataFrame."""
    if score_col not in src.columns:
        _print(f"    [WARN] '{score_col}' not found → {test_name} skipped")
        return pd.DataFrame()
    src = src.copy()
    safe_date(src, date_col)
    src = src.sort_values([ptid_col, date_col])
    scores_by_subj = {}
    for subj, grp in src.groupby(ptid_col):
        scores_by_subj[subj] = list(zip(grp[date_col], grp[score_col]))
    return build_nps_columns(scores_by_subj, test_name)


def build_all_clinical(dfs: dict) -> pd.DataFrame:
    """Neuropsychological, functional, behavioural and vital-sign series.

    The variables and their source columns come from ``variables.yaml``
    (``sources: ADNI:``). Visit date: ``EXAMDATE``, else ``VISDATE`` (column
    ``_date``), unless the source sets ``date_column`` (ECog: ``COLDATE``).
    """
    _print("  §4  Clinical / NPS / Functional / Vitals...")
    parts = []
    prepared: dict = {}
    for name, src in variable_sources(DATASET, config=_CTX.config):
        key = src.table
        if key not in prepared:
            df = dfs.get(key, pd.DataFrame())
            if not df.empty:
                if key == "neurobat":
                    df = df.copy()
                    combine_dates(df)
                    safe_date(df, "_date")
                else:
                    combine_dates(df)
            prepared[key] = df
        df = prepared[key]
        if df.empty:
            continue
        frame, col = materialise(df, name, src)
        if col is None:
            continue
        parts.append(_nps_from_source(frame, "PTID", src.date_column or "_date", name, col))

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()

    return reduce(lambda a, b: a.join(b, how="outer"), parts)


# =============================================================================
# §5  CSF BIOMARKERS
# =============================================================================
# Source: All_Subjects_UPENNBIOMK_ROCHE_ELECSYS_22Mar2026.csv
# Columns: PTID, VISCODE2, EXAMDATE, BATCH, ABETA40, ABETA42, TAU, PTAU
#
# BATCH encodes the assay platform / run and is used as the 'method' field.
# All values are numeric concentrations in pg/mL.
#
# Columns produced (per marker):
#   all_{marker}_values  — pipe-separated concentrations (chronological)
#   all_{marker}_method  — corresponding BATCH labels
#   nb_of_{marker}_values
# Plus:
#   all_pct_date_collection  — all lumbar-puncture dates
#   nb_of_date_collection
# =============================================================================

def build_csf(df_raw: pd.DataFrame) -> pd.DataFrame:
    """
    UPENN Roche Elecsys CSF panel.

    Legacy columns (kept for cohort_merge_integration_unified_v1):
      all_pct_date_collection, nb_of_date_collection,
      all_{m}_values, nb_of_{m}_values, all_{m}_method
    New / fixed:
      all_{m}_date     — collection date of EACH value (aligned with
                         all_{m}_values; previously only the shared
                         all_pct_date_collection existed, which is longer
                         than all_{m}_values whenever a marker is missing)
      all_{m}_method   — BATCH of each value (aligned; was a unique list)
      all_{m}_rundate  — assay run date of each value
      nb_of_{m}_censored — values reported as '<x' / '>x' (e.g. ABETA42
                         '>1700', PTAU '<8'); now kept numerically at the limit
                         instead of as text
      {m}_lab, {m}_kit
    Repeated runs of the same sample (same date) are averaged.
    """
    _print("  §5  CSF biomarkers...")
    if df_raw.empty:
        return pd.DataFrame()

    csf = df_raw.copy()
    safe_date(csf, "EXAMDATE")
    if "RUNDATE" in csf.columns:
        csf["RUNDATE"] = pd.to_datetime(csf["RUNDATE"], errors="coerce")
    csf = csf.sort_values("EXAMDATE")

    records = []
    for ptid, grp in csf.groupby("PTID"):
        row = {"ID": ptid}
        all_dates = list(dict.fromkeys(grp["EXAMDATE"].dropna().dt.strftime("%Y-%m-%d")))
        row["all_pct_date_collection"] = "|".join(all_dates) if all_dates else np.nan
        row["nb_of_date_collection"]   = len(all_dates)

        for marker in CSF_MARKERS:
            if marker not in grp.columns:
                continue
            g = grp.copy()
            g["_v"]    = plasma_clean_numeric(g[marker])
            g["_cens"] = g[marker].astype(str).str.strip().str.match(r"^[<>]").astype(int)
            g = g.dropna(subset=["_v"])
            if g.empty:
                continue
            g["_d"] = g["EXAMDATE"].dt.strftime("%Y-%m-%d").fillna("NA")
            agg = {"_v": ("_v", "mean"), "_cens": ("_cens", "max")}
            if "BATCH" in g.columns:
                agg["_b"] = ("BATCH", lambda s: "/".join(sorted(set(s.dropna().astype(str)))) or "NA")
            if "RUNDATE" in g.columns:
                agg["_r"] = ("RUNDATE", lambda s: "/".join(sorted(set(
                    s.dropna().dt.strftime("%Y-%m-%d")))) or "NA")
            a = g.groupby("_d", sort=False).agg(**agg).reset_index()

            row[f"all_{marker}_values"]    = "|".join(map(str, a["_v"].round(4)))
            row[f"all_{marker}_date"]      = "|".join(a["_d"])
            row[f"nb_of_{marker}_values"]  = len(a)
            row[f"nb_of_{marker}_censored"] = int(a["_cens"].sum())
            if "_b" in a.columns:
                row[f"all_{marker}_method"]  = "|".join(a["_b"])
            if "_r" in a.columns:
                row[f"all_{marker}_rundate"] = "|".join(a["_r"])
            row[f"{marker}_lab"] = CSF_MAIN_LAB
            row[f"{marker}_kit"] = CSF_MAIN_KIT

        records.append(row)

    return pd.DataFrame(records).set_index("ID")


def build_csf_extra(dfs: dict) -> pd.DataFrame:
    """Other ADNI CSF tables (AlzBio3, Lumipulse, EUROIMMUN, MSD, NfL, Ng …)."""
    _print("  §5a CSF — additional tables...")
    out = _adni_specs_to_frame(ADNI_CSF_EXTRA_SPECS, dfs, "extra CSF")
    return plasma_add_summary(out, prefix="CSF")


# =============================================================================
# §5b  PLASMA BIOMARKERS  +  §5c  ROUTINE BLOOD LAB
# =============================================================================
# Requires the "Plasma / lab helpers" cell (plasma_* functions).
# Every table listed in ADNI_PLASMA_SPECS / ADNI_LAB_SPEC is searched in
# ADNI_PLASMA_DIR; most plasma tables are keyed on RID only, so RID → PTID is
# mapped from DXSUM / PTDEMOG.
# =============================================================================

def _adni_rid_to_ptid(dfs: dict) -> dict:
    m = {}
    for key in ("dx", "demog"):
        df = dfs.get(key, pd.DataFrame())
        if not df.empty and {"RID", "PTID"} <= set(df.columns):
            sub = df[["RID", "PTID"]].dropna().drop_duplicates("RID")
            m.update(dict(zip(pd.to_numeric(sub["RID"], errors="coerce"), sub["PTID"])))
    return m


def _adni_attach_ptid(df: pd.DataFrame, rid_map: dict) -> pd.DataFrame:
    df = df.copy()
    if "RID" in df.columns:
        mapped = pd.to_numeric(df["RID"], errors="coerce").map(rid_map)
        df["_PTID"] = (df["PTID"].where(df["PTID"].notna(), mapped)
                       if "PTID" in df.columns else mapped)
    elif "PTID" in df.columns:
        df["_PTID"] = df["PTID"]
    else:
        df["_PTID"] = np.nan
    return df


def _adni_specs_to_frame(specs: dict, dfs: dict, label: str) -> pd.DataFrame:
    files = plasma_find_tables(ADNI_PLASMA_DIR, list(specs))
    if not files:
        _print(f"    [INFO] no {label} tables found in {ADNI_PLASMA_DIR}")
        return pd.DataFrame()
    rid_map = _adni_rid_to_ptid(dfs)
    parts = []
    for stem, spec in specs.items():
        path = files.get(stem)
        if not path:
            continue
        raw = plasma_read_csv(path)
        if raw.empty:
            continue
        raw = _adni_attach_ptid(raw, rid_map)
        n_unmapped = raw["_PTID"].isna().sum()
        time_col  = next((c for c in ("EXAMDATE", "DRAWDATE", "DRAWDTE", "DRAW_DATE",
                                      "SPECDATE", "SAMPLEDATE", "Date")
                          if c in raw.columns), None)
        visit_col = next((c for c in ("VISCODE2", "VISCODE") if c in raw.columns), None)
        got = plasma_from_spec(raw, spec, "_PTID", time_col, visit_col)
        parts.extend(got)
        _print(f"    [OK]   {stem}: {len(raw):,} rows → {len(got)} marker(s)"
              + (f"  ({n_unmapped} rows without PTID)" if n_unmapped else ""))
    return plasma_join(parts)


def build_plasma_adni(dfs: dict) -> pd.DataFrame:
    _print("  §5b Plasma biomarkers...")
    out = _adni_specs_to_frame(ADNI_PLASMA_SPECS, dfs, "plasma")
    return plasma_add_summary(out)


def build_lab_adni(dfs: dict) -> pd.DataFrame:
    _print("  §5c Routine blood lab (LABDATA + PETMETA glucose + Nightingale)...")
    out = _adni_specs_to_frame(ADNI_LAB_SPEC, dfs, "LABDATA")
    # Nightingale glucose is mmol/L → convert the whole series to mg/dL
    # (x 18.016) so it is on the same scale as RCT11 and PMBLGLUC.
    def _scale(v, k=18.016):
        if pd.isna(v):
            return v
        return "|".join("" if x in ("", "nan") else str(round(float(x) * k, 1))
                        for x in str(v).split("|"))
    for nm in ("Lab_Glucose_NMR", "Lab_Glucose_NMRv1"):
        if f"all_{nm}_score" in out.columns:
            out[f"all_{nm}_score"] = out[f"all_{nm}_score"].map(_scale)
            for pre in ("1st_", "Lst_"):
                c = f"{pre}{nm}_score"
                out[c] = (pd.to_numeric(out[c], errors="coerce") * 18.016).round(1)
            if f"{nm}_kit" in out.columns:
                out[f"{nm}_kit"] = out[f"{nm}_kit"].astype(str).str.replace(
                    "(mmol/L)", "(converted to mg/dL, x18.016)", regex=False)
    return out


# =============================================================================
# §5d  MEDICAL HISTORY — NEUROLOGIC / SEIZURE EVIDENCE  (added 2026-11-08)
# =============================================================================
# ADNI has no coded seizure variable. Sources (auto-discovered in ADNI_PLASMA_DIR):
#   MEDHIST    MH2NEURL  1/0 "Neurologic (other than AD)"    → Neuro_history (series)
#   RECMHIST   MHDESC    free-text problem description       ┐
#   INITHEALTH IHDESC    free-text condition (ADNI3/4)        ├→ SEIZURE_text_*
#   RECCMEDS   CMREASON  reason a medication was prescribed  ┘
#   RECCMEDS   CMMED/CMMEDO medication name (antiepileptics) → AED_text_*
# =============================================================================

def build_medhist_adni(dfs: dict) -> pd.DataFrame:
    _print("  §5d Medical history / seizure evidence (MEDHIST, RECMHIST, INITHEALTH, RECCMEDS)...")
    files = plasma_find_tables(ADNI_PLASMA_DIR, ADNI_MEDHIST_TABLES)
    if not files:
        _print(f"    [INFO] none of {ADNI_MEDHIST_TABLES} found in {ADNI_PLASMA_DIR}")
        return pd.DataFrame()
    rid_map = _adni_rid_to_ptid(dfs)
    tabs = {}
    for stem, path in files.items():
        raw = plasma_read_csv(path)
        if raw.empty:
            continue
        raw = _adni_attach_ptid(raw, rid_map)
        raw["_date"] = pd.NaT
        for c in ("EXAMDATE", "VISDATE", "USERDATE"):
            if c in raw.columns:
                raw["_date"] = raw["_date"].fillna(pd.to_datetime(raw[c], errors="coerce"))
        tabs[stem] = raw
        _print(f"    [OK]   {stem}: {len(raw):,} rows")

    parts = []
    mh = tabs.get("MEDHIST", pd.DataFrame())
    if not mh.empty and "MH2NEURL" in mh.columns:
        parts.append(plasma_series(mh, "_PTID", "MH2NEURL", "Neuro_history",
                                   "_date", "VISCODE2" if "VISCODE2" in mh.columns else None,
                                   lab="ADNI MEDHIST (1=yes, 0=no)"))

    # free-text seizure evidence across all history tables
    hist = []
    for stem, cols in (("RECMHIST", ["MHDESC"]),
                       ("INITHEALTH", ["IHDESC", "IHCOMM"]),
                       ("RECCMEDS", ["CMREASON", "CMREASONO", "CMCOMM"]),
                       ("MEDHIST", ["MHCOMMEN"])):
        t = tabs.get(stem, pd.DataFrame())
        keep = [c for c in cols if c in t.columns]
        if not t.empty and keep:
            hist.append(t[["_PTID", "_date"] + keep]
                        .assign(_txt=t[keep].fillna("").astype(str)
                                .agg(lambda r: " | ".join(x for x in r if x), axis=1))
                        [["_PTID", "_date", "_txt"]])
    if hist:
        h = pd.concat(hist, ignore_index=True)
        parts.append(text_flag_series(h, "_PTID", ["_txt"], SEIZURE_TEXT_RX,
                                      "SEIZURE_text", "_date"))
    cm = tabs.get("RECCMEDS", pd.DataFrame())
    if not cm.empty:
        parts.append(text_flag_series(cm, "_PTID", ["CMMED", "CMMEDO"], AED_TEXT_RX,
                                      "AED_text", "_date"))
    out = plasma_join(parts)
    for c in ("SEIZURE_text_flag", "AED_text_flag"):
        if c in out.columns:
            _print(f"    {c}: {int(out[c].fillna(0).sum()):,} subjects flagged")
    return out


# =============================================================================
# §6  DEMOGRAPHICS + APOE
# =============================================================================
# PTDEMOG columns used:
#   PTID, PTGENDER (1=Male 2=Female), PTRACCAT (race),
#   PTEDUCAT (years), PTDOBMM + PTDOBYY (birth date),
#   PTMARRY (marital status), PTWORK (employment)
#
# DOB reconstruction:
#   pd.Timestamp(year=PTDOBYY, month=PTDOBMM, day=1)
#   Age at first diagnosis = (1st_diag_date − DOB).days / 365.25
#   Calculated in §7 after the master merge.
#
# APOERES columns used:
#   PTID, GENOTYPE  (e.g. "3/4", "4/4", "3/3" …)
#   APOE_e4_count = number of '4' alleles in GENOTYPE string
# =============================================================================

def build_demographics(df_demog: pd.DataFrame,
                       df_entry: pd.DataFrame,
                       df_apoe: pd.DataFrame) -> pd.DataFrame:
    _print("  §6  Demographics + APOE...")
    records: dict = {}

    # ── PTDEMOG ──────────────────────────────────────────────────────────────
    # DOB strategy (in priority order):
    #   1. PTDOB  — full date string "YYYY-MM-DD" (present in PTDEMOG_24Feb2026)
    #   2. PTDOBYY — year-only integer; stored as "YYYY-01-01"
    #   3. PTDOBMM + PTDOBYY — split month/year (older ADNI exports)
    if not df_demog.empty:
        for _, row in df_demog.drop_duplicates("PTID").iterrows():
            ptid = row["PTID"]

            dob = pd.NaT
            # Priority 1: PTDOB full date
            ptdob = row.get("PTDOB", None)
            if pd.notna(ptdob) and str(ptdob).strip() not in ("", "nan"):
                dob = pd.to_datetime(str(ptdob).strip(), errors="coerce")

            # Priority 2: PTDOBYY as full date string (e.g. "1950-01-01")
            if pd.isna(dob):
                ptdobyy = row.get("PTDOBYY", None)
                if pd.notna(ptdobyy) and str(ptdobyy).strip() not in ("", "nan"):
                    raw = str(ptdobyy).strip()
                    # Could be "1950-01-01" or just "1950"
                    dob = pd.to_datetime(raw, errors="coerce")
                    if pd.isna(dob):
                        try:
                            dob = pd.Timestamp(year=int(float(raw)), month=1, day=1)
                        except Exception:
                            dob = pd.NaT

            # Priority 3: split PTDOBMM + PTDOBYY integers
            if pd.isna(dob):
                try:
                    yy = int(row.get("PTDOBYY", 0) or 0)
                    mm = max(1, int(row.get("PTDOBMM", 1) or 1))
                    if yy > 1800:
                        dob = pd.Timestamp(year=yy, month=mm, day=1)
                except Exception:
                    dob = pd.NaT

            records[ptid] = {
                "ID"             : ptid,
                "gender"         : row.get("PTGENDER", np.nan),   # 1=Male 2=Female
                "race"           : row.get("PTRACCAT", np.nan),
                "education"      : row.get("PTEDUCAT",  np.nan),
                "birth_date"     : dob,
                "marital_status" : row.get("PTMARRY",   np.nan),
                "employment"     : row.get("PTWORK",    np.nan),
            }

    # ── Study entry (fallback age + cohort entry group) ───────────────────────
    if not df_entry.empty:
        for _, row in df_entry.iterrows():
            sid = str(row.get("subject_id", "")).strip()
            if not sid:
                continue
            rec = records.setdefault(sid, {"ID": sid})
            rec["entry_age"]   = row.get("entry_age",            np.nan)
            rec["entry_date"]  = row.get("entry_date",           np.nan)
            rec["entry_group"] = row.get("entry_research_group", np.nan)

    # ── APOE ─────────────────────────────────────────────────────────────────
    if not df_apoe.empty:
        apoe = df_apoe.dropna(subset=["GENOTYPE"]).drop_duplicates("PTID")
        for _, row in apoe.iterrows():
            ptid = row["PTID"]
            rec  = records.setdefault(ptid, {"ID": ptid})
            geno = str(row["GENOTYPE"])
            rec["APOE_genotype"] = geno
            rec["APOE_e4_count"] = geno.count("4")

    df = pd.DataFrame(list(records.values()))
    return df.set_index("ID") if not df.empty else df


# =============================================================================
# §7  POST-MERGE: AGE AT FIRST DIAGNOSIS
# =============================================================================

def add_age_at_1st_diag(master: pd.DataFrame) -> pd.DataFrame:
    """
    Compute age_at_1st_diag = (1st_diag_date − birth_date) in years.
    Uses vectorised datetime arithmetic for speed; falls back row-wise on error.
    """
    if "1st_diag_date" not in master.columns or "birth_date" not in master.columns:
        return master

    bd = pd.to_datetime(master["birth_date"],   errors="coerce")
    d1 = pd.to_datetime(master["1st_diag_date"], errors="coerce")

    # Vectorised path (fast)
    try:
        days = (d1 - bd).dt.days
        master["age_at_1st_diag"] = (days / 365.25).round(2)
        # Sanity check: age must be 0–120
        bad = (master["age_at_1st_diag"] < 0) | (master["age_at_1st_diag"] > 120)
        master.loc[bad, "age_at_1st_diag"] = np.nan
    except Exception:
        # Row-wise fallback
        def _age(row):
            b = pd.to_datetime(row["birth_date"],    errors="coerce")
            d = pd.to_datetime(row["1st_diag_date"], errors="coerce")
            if pd.isna(b) or pd.isna(d):
                return np.nan
            age = round((d - b).days / 365.25, 2)
            return age if 0 <= age <= 120 else np.nan
        master["age_at_1st_diag"] = master.apply(_age, axis=1)

    return master


# =============================================================================
# MAIN
# =============================================================================

def main() -> pd.DataFrame:
    _print("=" * 62)
    _print("ADNI LONGITUDINAL INTEGRATION PIPELINE  v2")
    _print("=" * 62)

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load ──────────────────────────────────────────────────────────────────
    _print("\nLoading raw tables...")
    dfs = {k: resolve_id(load(k)) for k in PATHS}

    # ── Build domain tables ───────────────────────────────────────────────────
    _print("\nBuilding domain tables...")
    diag  = build_diagnosis(dfs["dx"])
    mri   = build_mri(dfs["mri"])
    pet   = build_pet(dfs["pet"])
    clin  = build_all_clinical(dfs)
    csf   = build_csf(dfs["csf"])
    csf_x  = build_csf_extra(dfs)
    plasma = build_plasma_adni(dfs)
    lab    = build_lab_adni(dfs)
    medh   = build_medhist_adni(dfs)
    demog = build_demographics(dfs["demog"], dfs["study_entry"], dfs["apoe"])

    # ── Merge ─────────────────────────────────────────────────────────────────
    _print("\nMerging all domains on subject ID...")
    parts = [p for p in [diag, mri, pet, clin, csf, csf_x, plasma, lab, medh, demog]
             if isinstance(p, pd.DataFrame) and not p.empty]

    if not parts:
        _print("ERROR: no data loaded. Check file paths.")
        return pd.DataFrame()

    master = reduce(lambda a, b: a.join(b, how="outer"), parts)
    master.index.name = "ID"
    master = master.reset_index()
    gc.collect()

    # ── Post-processing ───────────────────────────────────────────────────────
    master = add_age_at_1st_diag(master)

    # Column ordering: identity + demographics first
    master["Plasma_available"] = (
        master["Plasma_panels_disponible"].notna().map({True: "yes", False: "no"})
        if "Plasma_panels_disponible" in master.columns else "no")

    LEAD = [
        "ID",
        "gender", "race", "education", "birth_date",
        "age_at_1st_diag", "entry_age", "entry_group",
        "marital_status", "employment",
        "APOE_genotype", "APOE_e4_count",
        "1st_diag", "Lst_diag", "1st_diag_date", "Lst_diag_date",
        "all_diag", "all_diag_date", "nb_of_visits",
        "1st_time_of_change", "1st_time_of_change_days_to",
        "all_time_of_change", "all_time_of_change_days_to",
        "Lst_time_of_change", "nb_of_change", "Diag_total",
        "MRI_type_disponible", "PET_tracer_disponible",
        "all_pct_date_collection", "nb_of_date_collection",
    ]
    LEAD += ["Plasma_available", "Plasma_analytes_disponible", "Plasma_panels_disponible",
             "CSF_analytes_disponible", "CSF_panels_disponible"]
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

    # ── Summary print ─────────────────────────────────────────────────────────
    if "1st_diag" in master.columns:
        _print("\nDiagnosis distribution (1st visit):")
        _print(master["1st_diag"].value_counts().to_string())

    if "MRI_type_disponible" in master.columns:
        _print("\nMRI types (subjects with >= 1 acquisition):")
        cnt: Counter = Counter()
        for v in master["MRI_type_disponible"].dropna():
            cnt.update(v.split("|"))
        for t, n in cnt.most_common():
            _print(f"  {t:<30s} {n:>5,} subjects")

    if "PET_tracer_disponible" in master.columns:
        _print("\nPET tracers:")
        cnt = Counter()
        for v in master["PET_tracer_disponible"].dropna():
            cnt.update(v.split("|"))
        for t, n in cnt.most_common():
            _print(f"  {t:<30s} {n:>5,} subjects")

    _print("\nPlasma coverage (subjects with >= 1 value):")
    for c in sorted(c for c in master.columns
                    if c.startswith("nb_of_Plasma_") and c.endswith("_visits")):
        n = master[c].fillna(0).gt(0).sum()
        _print(f"  {c[6:-7]:<40s} {n:>5,} subjects")

    _print("\nCSF coverage:")
    for m in CSF_MARKERS:
        col = f"nb_of_{m}_values"
        if col in master.columns:
            n = master[col].gt(0).sum()
            _print(f"  {m:<10s} {n:>5,} subjects with >= 1 value")
        cc = f"nb_of_{m}_censored"
        if cc in master.columns:
            _print(f"  {'':<10s} {int(master[cc].fillna(0).sum()):>5,} censored values (<LoQ / >ULoQ)")
    for c in sorted(c for c in master.columns
                    if c.startswith("nb_of_CSF_") and c.endswith("_visits")):
        n = master[c].fillna(0).gt(0).sum()
        _print(f"  {c[6:-7]:<30s} {n:>5,} subjects")

    return master


# =============================================================================
