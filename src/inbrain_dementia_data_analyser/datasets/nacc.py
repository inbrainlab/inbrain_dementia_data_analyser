"""NACC LONGITUDINAL INTEGRATION PIPELINE
Applies the same domain-builder logic as ADNI v2, adapted for NACC's
large main UDS file with chunk-based reading and separate scan tables.

NACC-SPECIFIC ARCHITECTURE
---------------------------
NACC ships 6 CSV files in Tables_NACC/:

  investigator_ftldlbd_nacc70.csv       — MAIN file: all UDS forms merged
                                          (demographics + diagnosis + NPS +
                                           vitals + medical history)
                                          One row per subject per visit.
  investigator_fcsf_nacc70.csv          — CSF: CSFABETA, CSFPTAU, CSFTTAU
  investigator_scan_mriqc_nacc70.csv    — MRI: SERIESTYPE, SERIESDESCRIPTION,
                                          STUDYDATE, SCANNERMANUFACTURER
  investigator_scan_amyloidpetgaain_nacc71.csv  — Amyloid PET (GAAIN method)
                                          CENTILOIDS, AMYLOID_STATUS, TRACER
  investigator_scan_amyloidpetnpdka_nacc71.csv  — Amyloid PET (NPDKA method)
                                          regional SUVR columns
  investigator_scan_taupetnpdka_nacc71.csv      — Tau PET (NPDKA method)
                                          META_TEMPORAL_SUVR + regional

SUBJECT ID:
  Main file: NACCID → normalised to PTID
  PET/MRI scan files: NACCID
  CSF file: NACCID (NACCADC is the ADC center code, not subject ID)

TIME AXIS:
  VISITYR + VISITMO + VISITDAY → _date (via make_visit_date())
  NACCVNUM: visit number (1=first, 2=second, …)
  CSF: CSFLPYR + CSFLPMO + CSFLPDY → _date

CHUNK READING:
  investigator_ftldlbd_nacc70.csv is very large (~54,000+ rows, 700+ cols).
  It is read in chunks via read_main_chunked() — loaded ONCE per section,
  restricted to the needed columns.
  CHUNK_SIZE = 10,000 rows per chunk.

NACC SENTINEL CODES (replace with NaN):
  All numeric NPS and diagnosis columns: values < 0 → NaN
  Common codes: -4 = not available, -2 = unknown
  888 / 999 handled by < 0 rule not applying, but flagged if seen.

DIAGNOSIS CODING:
  Uses NACC-derived variables (preferred over raw UDS flags):
    NORMCOG = 1          → CN
    NACCTMCI = 1         → MCI
    NACCALZD = 1         → AD
    NACCLBDS = 1         → LBD
    NACCPPA = 1          → PPA (FTD spectrum)
    FTDIF = 1            → FTD
    PSPIF = 1            → PSP
    MSAIF = 1            → MSA
    CORTIF = 1           → CBS (Corticobasal syndrome)
    CVDIF / VASCIF = 1   → VaD
    DEMENTED=1 (no match) → Dementia-NOS

  Falls back to full UDS D1 boolean flags (PROBAD, POSSAD, DLB, VASC, FTD,
  PSP, MCIAMEM, etc.) when NACC-derived variables are absent.
  Compound labels via sorted '+' join (same as all other pipelines).

NPS VARIABLES (from main file):
  NACC uses derived variables (NACCMMSE, NACCGDS, NACCMOCA) preferred over
  raw form scores; raw scores (BOSTON, ANIMALS, LOGIMEM, MEMUNITS, DIGIF,
  DIGIB, WAIS) used where NACC-derived are absent.
  TMT: TRAILA/OTRAILA and TRAILB/OTRAILB (v2/v3 variants) → combined.
  RAVLT: REY1REC-REY5REC (trials 1-5), REYDREC (delayed recall), REYTCOR (recognition).
  All NPS: negative values → NaN.

MRI: investigator_scan_mriqc_nacc70.csv
  SERIESTYPE → canonical via MRI_TYPE_MAP
  STUDYDATE → date
  SERIESTYPE values: T1w, T2w, FLAIR, DWI, fMRI, ASL, etc.

PET: Three separate scan files
  Amyloid GAAIN: TRACER, SCANDATE, CENTILOIDS, AMYLOID_STATUS
  Amyloid NPDKA: TRACER, SCANDATE, NPDKA_SUMMARY_SUVR + regional
  Tau NPDKA:     TRACER, SCANDATE, META_TEMPORAL_SUVR + regional
  Tracers from NACC: PIB, AV45, Flortaucipir, FBB, FDG

CSF: investigator_fcsf_nacc70.csv
  CSFABETA, CSFPTAU, CSFTTAU
  Date from CSFLPYR/CSFLPMO/CSFLPDY

Demographics (from main file, first visit per subject):
  BIRTHMO, BIRTHYR → birth_date
  SEX (1=M 2=F), EDUC, RACE, HISPANIC, MARISTAT, HANDED,
  NACCAGE (age at visit), NACCAPOE (APOE genotype code)
  NACCAPOE coding: 1=e3/e3, 2=e3/e4, 3=e3/e2, 4=e4/e4, 5=e4/e2, 6=e2/e2, 9=missing

Output:
  <root>/Tables_NACC/Integrated_NACC/NACC_longitudinal.csv

Author  : VRS_inbrain project
Version : 1  —  2026-06-09

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("NACC", root_path)
"""
from __future__ import annotations

import gc
import numpy as np
import os
import pandas as pd
from functools import reduce

from ..datasets.base import DatasetContext
from ..diagnosis import (UDS_ETIOLOGIES, UDS_MCI_COLS, add_dx_columns,
                         nacc_stage_etiology, trajectory_record)
from ..fluid import plasma_series
from ..imaging.scanner import add_scanner_channels, norm_vendor
from ..imaging.sessions import session_id_columns
from ..series import build_scored_columns, d_fmt, list_str
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

DATASET = "NACC"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
NACC = OUT_DIR = OUT_LONG = ""


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, NACC, OUT_DIR, OUT_LONG, CHUNK_SIZE
    _CTX = ctx
    NACC = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()
    CHUNK_SIZE = int(ctx.settings.get("chunk_size", 10_000))


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write the per-subject table of NACC; returns it."""
    configure(ctx)
    return main()



def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns of one measure (see series.build_scored_columns)."""
    return build_scored_columns(scores_by_subj, test_name)




# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================











CHUNK_SIZE = 10_000   # rows per chunk when reading main file

# =============================================================================
# ── COLUMN DEFINITIONS ────────────────────────────────────────────────────────
# =============================================================================

# Identifiers + date
ID_DATE_COLS = ["NACCID", "NACCVNUM", "VISITYR", "VISITMO", "VISITDAY"]

# Diagnosis columns — NACC-derived + full UDS D1 flags
DX_DERIVED_COLS = [
    "NORMCOG",
    "NACCTMCI",   # MCI (NACC derived type)
    "NACCALZD",   # Alzheimer dementia (NACC derived)
    "NACCLBDS",   # LBD spectrum (NACC derived)
    "NACCPPA",    # Primary Progressive Aphasia (NACC derived)
    "FTDIF",      # FTD impairment flag
    "PSPIF",      # PSP impairment flag
    "MSAIF",      # MSA impairment flag
    "CORTIF",     # Corticobasal syndrome
    "CVDIF",      # Cerebrovascular disease
    "VASCIF",     # Vascular dementia
    'EPILEPIF',
    "NACCUDSD",   # v0.4.0: cognitive stage (1 CN, 2 impaired-not-MCI, 3 MCI, 4 dementia)
    "NACCETPR",   # v0.4.0: primary etiologic diagnosis (one per visit)
]

# Full UDS D1 boolean flags (fallback when NACC-derived absent)
DX_UDS_COLS = [
    "DEMENTED",
    "PROBAD", "POSSAD",
    "NACCLBDE",
    "VASC", "VASCPS",
    "FTD", "PPAPH", "PNAPH", "SEMDEMAN", "SEMDEMAG", "PPAOTHR",
    "PSP", "CORT", "HUNT", "PARK",
    "MCIAMEM", "MCIAPLUS", "MCIAPLAN", "MCIAPATT", "MCIAPEX", "MCIAPVIS",
    "MCINON1", "MCIN1LAN", "MCIN1ATT", "MCIN1EX", "MCIN1VIS",
    "MCINON2", "MCIN2LAN", "MCIN2ATT", "MCIN2EX", "MCIN2VIS", "IMPNOMCI", 'EPILEP',
]
# v0.4.0: the UDS D1 fallback also needs the primary / contributing (…IF) columns
DX_FLAG_COLS = list(dict.fromkeys(
    DX_DERIVED_COLS + DX_UDS_COLS
    + [c for _e, flag_, if_ in UDS_ETIOLOGIES for c in (flag_, if_)] + UDS_MCI_COLS))

# NPS columns (all numeric; negative values → NaN)
NPS_NUMERIC_COLS = [
    "NACCMMSE", "NACCMOCA", "NACCGDS",
    "CDRGLOB", "CDRSUM",
    "BOSTON", "MINTTOTS",
    "ANIMALS",
    "TRAILA", "OTRAILA",
    "TRAILB", "OTRAILB",
    "LOGIMEM", "MEMUNITS", "MEMTIME",
    "DIGIF", "DIGIB", "WAIS",
    "REY1REC", "REY2REC", "REY3REC", "REY4REC", "REY5REC",
    "REYDREC", "REYTCOR",
    "CRAFTVRS", "CRAFTDVR",
    "CRAFTURS", "CRAFTDRE",
]

# Demographics columns from main file
DEMOG_COLS = [
    "NACCAGE", "SEX", "EDUC", "RACE", "HISPANIC",
    "MARISTAT", "HANDED", "NACCAPOE",
    "BIRTHMO", "BIRTHYR",
    "NACCLIVS",   # living situation (NACC derived)
    "INDEPEND",   # independence level
    "RESIDENC",   # residence type
    "HEIGHT", "WEIGHT", "BPSYS", "BPDIAS", "HRATE",
    # Family history
    "NACCFAM",    # any first-degree relative with dementia
    "NACCMOM",    # mother with dementia
    "NACCDAD",    # father with dementia
    "NACCFADM",   # familial AD mutation
    "NACCFFTD",   # familial FTLD mutation
]

# NACCAPOE genotype codes (NACC RDD):
#   1 = e3,e3   2 = e3,e4   3 = e3,e2   4 = e4,e4   5 = e4,e2   6 = e2,e2
#   9 = missing / unknown / not assessed
# (until 2026-11-16 this map read 1=e2/e2 … 6=e4/e4, which labelled the most
#  common genotype, e3/e3, as e2/e2 and under-counted the e4 carriers)
APOE_CODE_MAP = {
    1: "e3/e3", 2: "e3/e4", 3: "e2/e3",
    4: "e4/e4", 5: "e2/e4", 6: "e2/e2",
}
APOE_E4_COUNT = {
    1: 0, 2: 1, 3: 0, 4: 2, 5: 1, 6: 0,
}

# MRI series type → canonical (from NACC SERIESTYPE column)
MRI_TYPE_MAP = {
    "t1w"     : "MRI_T1w",
    "t2w"     : "MRI_T2w",
    "flair"   : "MRI_T2w_FLAIR",
    "t2*w"    : "MRI_T2star",
    "swi"     : "MRI_T2star",
    "dwi"     : "MRI_DTI",
    "dti"     : "MRI_DTI",
    "bold"    : "MRI_fMRI",
    "rsfmri"  : "MRI_fMRI",
    "asl"     : "MRI_ASL",
    "perf"    : "MRI_ASL",
}

# PET tracer → canonical
# NACC SCAN PET RDD — TRACER codes (amyloid GAAIN/NPDKA, tau, FDG files)
# petqc SCANNER_MODEL codes (SCAN PET RDD)
SCAN_PET_SCANNER = {
    1: "Biograph Vision", 2: "Biograph Horizon", 3: "Biograph TruPoint",
    4: "Biograph mCT", 5: "Biograph mMR", 6: "ECAT HR+", 7: "Discovery MI",
    8: "Discovery MI DR", 9: "Discovery 610", 10: "Discovery 710",
    11: "Discovery 690XT", 12: "Discovery STE", 13: "Discovery ST",
    14: "SIGNA PET/MR", 15: "Ingenuity TF", 16: "Vereos",
    17: "United Imaging uMI 550",
}

SCAN_TRACER_CODES = {
    1: "PET_FDG",           2: "ABETA_PIB",      3: "ABETA_AV45",        # FBP
    4: "ABETA_FBB",         5: "ABETA_NAV4694",  6: "TAU_Flortaucipir",
    7: "TAU_MK6240",        8: "TAU_PI2620",     9: "TAU_GTP1",
    10: "ABETA_Flutemetamol",
}

PET_TRACER_MAP = {
    "pib"          : "ABETA_PIB",
    "11c-pib"      : "ABETA_PIB",
    "c11-pib"      : "ABETA_PIB",
    "fdg"          : "Metabolism",
    "18f-fdg"      : "Metabolism",
    "av45"         : "ABETA_AV45",
    "florbetapir"  : "ABETA_AV45",
    "18f-av45"     : "ABETA_AV45",
    "fbb"          : "ABETA_FBB",
    "florbetaben"  : "ABETA_FBB",
    "18f-fbb"      : "ABETA_FBB",
    "flortaucipir" : "TAU_Flortaucipir",
    "av1451"       : "TAU_Flortaucipir",
    "18f-av1451"   : "TAU_Flortaucipir",
    "mk6240"       : "TAU_MK6240",
    "18f-mk6240"   : "TAU_MK6240",
    "nav4694"      : "ABETA_NAV4694",
    "18f-nav4694"  : "ABETA_NAV4694",
}

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================







def make_visit_date(df: pd.DataFrame) -> pd.DataFrame:
    """Build _date from VISITYR/VISITMO/VISITDAY (same as ABVIB pattern)."""
    if all(c in df.columns for c in ["VISITYR", "VISITMO", "VISITDAY"]):
        yr = pd.to_numeric(df["VISITYR"],  errors="coerce").fillna(0).astype(int)
        mo = pd.to_numeric(df["VISITMO"],  errors="coerce").fillna(1).astype(int).clip(1, 12)
        dy = pd.to_numeric(df["VISITDAY"], errors="coerce").fillna(1).astype(int).clip(1, 28)
        df["_date"] = pd.to_datetime(dict(year=yr, month=mo, day=dy),
                                      errors="coerce")
    else:
        df["_date"] = pd.NaT
    return df


def make_csf_date(df: pd.DataFrame) -> pd.DataFrame:
    """Build _date from CSFLPYR/CSFLPMO/CSFLPDY (CSF lumbar puncture date)."""
    if all(c in df.columns for c in ["CSFLPYR", "CSFLPMO", "CSFLPDY"]):
        yr = pd.to_numeric(df["CSFLPYR"], errors="coerce").fillna(0).astype(int)
        mo = pd.to_numeric(df["CSFLPMO"], errors="coerce").fillna(1).astype(int).clip(1, 12)
        dy = pd.to_numeric(df["CSFLPDY"], errors="coerce").fillna(1).astype(int).clip(1, 28)
        df["_date"] = pd.to_datetime(dict(year=yr, month=mo, day=dy),
                                      errors="coerce")
    else:
        df["_date"] = pd.NaT
    return df


def clean_negatives(df: pd.DataFrame, cols: list) -> pd.DataFrame:
    """Replace NACC sentinel negative codes with NaN in specified columns."""
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
            df.loc[df[col] < 0, col] = np.nan
    return df


def load_small(key: str) -> pd.DataFrame:
    """Load a small supplementary CSV (MRI, PET, CSF) in one shot."""
    path = PATHS.get(key, "")
    if not path or not os.path.exists(path):
        _print(f"  [SKIP] {key}: not found")
        return pd.DataFrame()
    for i, kwargs in enumerate([
        dict(engine="c",      on_bad_lines="error", encoding="utf-8"),
        dict(engine="c",      on_bad_lines="skip",  encoding="utf-8"),
        dict(engine="python", on_bad_lines="skip",  encoding="utf-8"),
        dict(engine="python", on_bad_lines="skip",  encoding="latin-1"),
    ]):
        try:
            df = pd.read_csv(path, low_memory=False, **kwargs)
            if i > 0:
                _print(f"  [WARN] {key}: loaded via fallback strategy {i+1}")
            _print(f"  [OK]   {key}: {len(df):>7,} rows  {df.shape[1]} cols")
            return df
        except Exception:
            continue
    _print(f"  [ERROR] {key}: all strategies failed")
    return pd.DataFrame()


# =============================================================================
# ── CHUNK READER (main file only) ─────────────────────────────────────────────
# =============================================================================

def read_main_chunked(keep_cols: list | None = None) -> pd.DataFrame:
    """
    Read investigator_ftldlbd_nacc70.csv in CHUNK_SIZE-row chunks.

    Steps per chunk:
      1. Restrict to keep_cols (where available)
      2. Replace NACC negative sentinel codes with NaN in NPS + DX cols
      3. Build _date from VISITYR/VISITMO/VISITDAY

    Returns a single concatenated DataFrame.
    """
    path = PATHS["main"]
    if not os.path.exists(path):
        _print(f"  [ERROR] Main NACC file not found → {path}")
        return pd.DataFrame()

    # Read header to discover available columns
    header = pd.read_csv(path, nrows=0)
    available = list(header.columns)

    if keep_cols is None:
        use_cols = available
    else:
        use_cols = [c for c in keep_cols if c in available]
        missing  = [c for c in keep_cols if c not in available]
        if missing:
            _print(f"    [WARN] Columns not in main file: {missing}")

    nps_to_clean = [c for c in NPS_NUMERIC_COLS if c in use_cols]
    dx_to_clean  = [c for c in DX_FLAG_COLS      if c in use_cols]

    chunks = []
    n_rows = 0
    for chunk in pd.read_csv(path, usecols=use_cols,
                              chunksize=CHUNK_SIZE, low_memory=False):
        chunk = clean_negatives(chunk, nps_to_clean + dx_to_clean)
        chunk = make_visit_date(chunk)
        chunks.append(chunk)
        n_rows += len(chunk)

    if not chunks:
        return pd.DataFrame()

    df = pd.concat(chunks, ignore_index=True)
    gc.collect()
    _print(f"  [OK]   main: {n_rows:>7,} rows  {df.shape[1]} cols (chunked)")
    return df





# =============================================================================
# §1  DIAGNOSIS TRAJECTORY
# =============================================================================

MCI_COLS_UDS = [
    "MCIAMEM", "MCIAPLUS", "MCIAPLAN", "MCIAPATT", "MCIAPEX", "MCIAPVIS",
    "MCINON1", "MCIN1LAN", "MCIN1ATT", "MCIN1EX",  "MCIN1VIS",
    "MCINON2", "MCIN2LAN", "MCIN2ATT", "MCIN2EX",  "MCIN2VIS", "IMPNOMCI",
]


def _classify_nacc(row: pd.Series) -> str:
    """
    Compound-label classifier for NACC.
    Priority 1: NORMCOG=1 → CN (sole priority, no compound)
    Priority 2: NACC-derived variables (NACCTMCI, NACCALZD, NACCLBDS,
                NACCPPA, FTDIF, PSPIF, MSAIF, CORTIF, CVDIF)
    Priority 3: Full UDS D1 boolean flags (fallback)
    """
    def flag(col):
        v = row.get(col, None)
        return pd.notna(v) and float(v or 0) == 1

    # NORMCOG: sole priority — but epilepsy is recorded in D1 even when
    # cognition is normal, so keep it: "CN+EPY"
    if flag("NORMCOG"):
        return "CN+EPY" if (flag("EPILEP") or flag("EPILEPIF")) else "CN"

    parts = []

    # NACC-derived variables (more reliable, use first)
    if flag("NACCTMCI"):  parts.append("MCI")
    if flag("NACCALZD"):  parts.append("AD")
    if flag("NACCLBDS"):  parts.append("LBD")
    if flag("NACCPPA"):   parts.append("FTD")   # PPA is FTD spectrum
    if flag("FTDIF"):     parts.append("FTD")
    if flag("PSPIF"):     parts.append("PSP")
    if flag("MSAIF"):     parts.append("MSA")
    if flag("CORTIF"):    parts.append("CBS")
    if flag("EPILEP"):    parts.append("EPY")
    if flag("CVDIF") or flag("VASCIF"): parts.append("VaD")

    # If NACC-derived gave nothing, fall back to full UDS D1
    if not parts:
        if any(flag(c) for c in MCI_COLS_UDS):
            parts.append("MCI")
        if flag("PROBAD") or flag("POSSAD"):
            parts.append("AD")
        if flag("NACCLBDE"):
            parts.append("LBD")
        if flag("VASC") or flag("VASCPS"):
            parts.append("VaD")
        if flag("FTD") or flag("PSP") or flag("CORT") or \
           flag("PPAPH") or flag("PNAPH") or flag("SEMDEMAN"):
            parts.append("FTD")
        if flag("HUNT"):
            parts.append("HD")
        if flag("PARK"):
            parts.append("PD")
        if flag("EPILEPIF"):
            parts.append("EPY")

    if not parts:
        if flag("DEMENTED"):
            return "Dementia-NOS"
        return "Unknown"

    return "+".join(sorted(set(parts)))


def build_diagnosis(df: pd.DataFrame) -> pd.DataFrame:
    _print("  §1  Diagnosis trajectory...")
    df = df.sort_values(["NACCID", "_date"])
    df["DX_LABEL"] = df.apply(_classify_nacc, axis=1)    # pre-0.4.0 label -> *_diag_legacy
    df = add_dx_columns(df, nacc_stage_etiology)

    records = []
    for ptid, grp in df.groupby("NACCID"):
        grp    = grp.sort_values("_date")
        dates  = grp["_date"].tolist()
        records.append(trajectory_record(ptid, grp, dates))

    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §2  MRI ACQUISITIONS
# =============================================================================
# Source: investigator_scan_mriqc_nacc70.csv
# Key columns: NACCID, STUDYDATE, SERIESTYPE, SERIESDESCRIPTION,
#              SCANNERMANUFACTURER, SCANNERMODEL, SERIESCHOSEN (=1 chosen)
#
# SERIESTYPE → canonical via MRI_TYPE_MAP (lowercase substring match)
# Date from STUDYDATE (format: YYYYMMDD integer or string)
# =============================================================================

def _classify_mri_series(series_type: str,
                          series_desc: str = "") -> str | None:
    """Map SERIESTYPE (and optionally SERIESDESCRIPTION) to canonical type."""
    for s in [series_type, series_desc]:
        if not isinstance(s, str) or not s.strip():
            continue
        d = s.lower().replace("-", "_").replace(" ", "_")
        # Skip scouts and field maps
        if any(x in d for x in ["scout", "localizer", "field_map",
                                  "gre_field", "aahscout"]):
            return None
        # Skip derived DTI maps
        if any(d.endswith(x) for x in ["_fa", "_adc", "_tracew"]):
            return None
        for key, canon in MRI_TYPE_MAP.items():
            if key in d:
                return canon
    return None


def build_mri(df_mri: pd.DataFrame) -> pd.DataFrame:
    _print("  §2  MRI acquisitions (MRIQC scan table)...")
    if df_mri.empty:
        return pd.DataFrame()

    mri = df_mri.copy()
    mri.columns = [c.strip() for c in mri.columns]

    # Parse date (STUDYDATE is YYYYMMDD integer)
    if "STUDYDATE" in mri.columns:
        mri["_date"] = pd.to_datetime(mri["STUDYDATE"], errors="coerce")

    # Classify series type
    st_col   = "SERIESTYPE"      if "SERIESTYPE"      in mri.columns else None
    desc_col = "SERIESDESCRIPTION" if "SERIESDESCRIPTION" in mri.columns else None

    mri["_canon"] = mri.apply(
        lambda r: _classify_mri_series(
            str(r.get(st_col, "")) if st_col else "",
            str(r.get(desc_col, "")) if desc_col else ""
        ),
        axis=1
    )
    mri = mri.dropna(subset=["_canon"])
    if mri.empty:
        return pd.DataFrame()

    # Filter to chosen series where flag available
    if "SERIESCHOSEN" in mri.columns:
        chosen = mri[mri["SERIESCHOSEN"] == 1]
        if not chosen.empty:
            mri = chosen

    mri = mri.sort_values(["NACCID", "_date"])

    records = {}
    avail_key = "MRI_type_disponible"

    for ptid, grp in mri.groupby("NACCID"):
        row = records.setdefault(ptid, {"ID": ptid, avail_key: set()})
        for _, r in grp.iterrows():
            mtype = r["_canon"]
            date  = r["_date"]
            loni  = next((r.get(c) for c in ("LONIIMAGE", "LONIUID", "LONISERIES")
                          if pd.notna(r.get(c, np.nan))), None)
            img_id = str(loni) if pd.notna(loni) else f"{ptid}_{d_fmt(date)}"
            row[avail_key].add(mtype)
            row.setdefault(f"all_{mtype}_date_list", []).append(d_fmt(date))
            row.setdefault(f"all_{mtype}_id_list",   []).append(img_id)

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
    return add_scanner_channels(result.set_index("ID"), mri, "NACCID", "_canon", "_date",
                                "SCANNERMANUFACTURER", "SCANNERMODEL")


# =============================================================================
# §3  PET ACQUISITIONS
# =============================================================================
# Three PET files, all keyed on NACCID + SCANDATE + TRACER
#
# GAAIN:  CENTILOIDS, AMYLOID_STATUS, GAAIN_SUMMARY_SUVR
# NPDKA:  NPDKA_SUMMARY_SUVR + regional (amyloid)
# TAU:    META_TEMPORAL_SUVR + regional
# =============================================================================

def _norm_tracer(raw) -> str | None:
    """NACC SCAN TRACER is a numeric code (SCAN PET RDD); names also accepted."""
    if pd.isna(raw):
        return None
    try:
        return SCAN_TRACER_CODES.get(int(float(raw)))
    except (TypeError, ValueError):
        pass
    t = str(raw).strip().lower()
    for key, canon in PET_TRACER_MAP.items():
        if key in t:
            return canon
    return None


def _pet_date(df: pd.DataFrame) -> pd.DataFrame:
    """SCANDATE is 'YYYY-MM-DD' in SCAN files (older exports: YYYYMMDD)."""
    if "SCANDATE" in df.columns:
        s = df["SCANDATE"].astype(str).str.strip()
        d = pd.to_datetime(s, format="%Y-%m-%d", errors="coerce")
        d = d.fillna(pd.to_datetime(s.str.replace(r"\.0$", "", regex=True).str.zfill(8),
                                    format="%Y%m%d", errors="coerce"))
        df["_date"] = d
    else:
        df["_date"] = pd.NaT
    return df


def _aligned(series) -> str | float:
    """Pipe list that KEEPS positions (missing -> 'NA'), aligned with dates/IDs."""
    vals = ["NA" if pd.isna(v) else str(v) for v in series]
    return "|".join(vals) if vals else np.nan


def build_pet(dfs: dict) -> pd.DataFrame:
    """
    SCAN PET — one block per tracer (ABETA_*, TAU_*, PET_FDG):
      all_{T}_date, nb_total_{T}, {T}_IDs (LONIUID), 1st/Lst_{T}_ID, 1st/Lst_{T}_date
      all_{T}_suvr       amyloid: GAAIN summary SUVR | tau: META_TEMPORAL_SUVR
                         | FDG: FDG_METAROI_SUVR
      all_{T}_centiloid  amyloid only (GAAIN)
      all_{T}_status     amyloid positivity (GAAIN AMYLOID_STATUS)
      all_{T}_suvr_npdka amyloid NPDKA summary SUVR (joined by LONIUID)
      all_{T}_entorhinal_suvr  tau only
      all_{T}_qc         PASS_FAIL from petqc (1 pass / 0 fail), when available
    Every all_* list is aligned scan by scan with all_{T}_date / {T}_IDs.
    """
    _print("  §3  PET acquisitions (SCAN amyloid + tau + FDG)...")

    qc = dfs.get("pet_qc", pd.DataFrame())
    # petqc uses the RAW image LONIUID (≠ processed-scan LONIUID) -> join on
    # NACCID + scan date + tracer code
    qc_map = {}
    if not qc.empty and {"NACCID", "SCAN_DATE", "RADIOTRACER", "PASS_FAIL"} <= set(qc.columns):
        sm = (pd.to_numeric(qc["SCANNER_MODEL"], errors="coerce") if "SCANNER_MODEL" in qc.columns
              else pd.Series(np.nan, index=qc.index))
        sm = [SCAN_PET_SCANNER.get(int(x), f"code {int(x)}") if pd.notna(x) else np.nan for x in sm]
        qd = _pet_date(qc.rename(columns={"SCAN_DATE": "SCANDATE"}))["_date"]
        qc_map = {(n, d_fmt(d), str(int(t))): (v, m) for n, d, t, v, m in zip(
            qc["NACCID"], qd, pd.to_numeric(qc["RADIOTRACER"], errors="coerce").fillna(-1),
            pd.to_numeric(qc["PASS_FAIL"], errors="coerce"), sm)}

    npdka = dfs.get("pet_npdka", pd.DataFrame())
    npdka_map = {}
    if not npdka.empty and "NPDKA_SUMMARY_SUVR" in npdka.columns:
        npdka_map = dict(zip(npdka["LONIUID"].astype(str), npdka["NPDKA_SUMMARY_SUVR"]))

    frames = []
    gaain = dfs.get("pet_gaain", pd.DataFrame())
    if not gaain.empty:
        g = gaain.copy()
        g["_suvr"]  = g.get("GAAIN_SUMMARY_SUVR")
        g["_cl"]    = g.get("CENTILOIDS")
        g["_stat"]  = g.get("AMYLOID_STATUS")
        g["_npdka"] = g["LONIUID"].astype(str).map(npdka_map)
        frames.append(g)
    elif not npdka.empty:                     # NPDKA only (no GAAIN file)
        g = npdka.copy()
        g["_suvr"] = g.get("NPDKA_SUMMARY_SUVR")
        frames.append(g)
    tau = dfs.get("pet_tau", pd.DataFrame())
    if not tau.empty:
        t = tau.copy()
        t["_suvr"] = t.get("META_TEMPORAL_SUVR")
        t["_ento"] = t.get("CTX_ENTORHINAL_SUVR")
        frames.append(t)
    fdg = dfs.get("pet_fdg", pd.DataFrame())
    if not fdg.empty:
        f = fdg.copy()
        f["_suvr"] = f.get("FDG_METAROI_SUVR")
        frames.append(f)

    if not frames:
        return pd.DataFrame()
    pet = pd.concat(frames, ignore_index=True, sort=False)
    for c in ("_suvr", "_cl", "_stat", "_npdka", "_ento"):   # absent when GAAIN is missing
        if c not in pet.columns:
            pet[c] = np.nan
    pet = _pet_date(pet)
    pet["_tracer"] = pet["TRACER"].apply(_norm_tracer) if "TRACER" in pet.columns else None
    n_unk = pet["_tracer"].isna().sum()
    if n_unk:
        _print(f"    [WARN] {n_unk} PET rows with unknown TRACER code: "
              f"{pet.loc[pet['_tracer'].isna(), 'TRACER'].value_counts().to_dict()}")
    pet = pet.dropna(subset=["_tracer"])
    pet["_id"] = pet["LONIUID"].astype(str)
    if qc_map:
        tcode = pd.to_numeric(pet["TRACER"], errors="coerce").fillna(-1).astype(int).astype(str)
        keys = list(zip(pet["NACCID"], pet["_date"].map(d_fmt), tcode))
        hits = [qc_map.get(k, (np.nan, np.nan)) for k in keys]
        pet["_qc"]    = [h[0] for h in hits]
        pet["_model"] = [h[1] for h in hits]
        n_hit = sum(k in qc_map for k in keys)
        _print(f"    petqc matched {n_hit}/{len(keys)} scans (NACCID + date + tracer)")
    else:
        pet["_qc"] = pet["_model"] = np.nan
    pet["_mfr"] = pet["_model"].map(norm_vendor)
    pet = pet.drop_duplicates(["NACCID", "_tracer", "_id"]).sort_values(["NACCID", "_date"])
    _print("    scans per tracer: " + str(pet["_tracer"].value_counts().to_dict()))

    records = {}
    for (ptid, tracer), tg in pet.groupby(["NACCID", "_tracer"], sort=False):
        r = records.setdefault(ptid, {"ID": ptid, "PET_tracer_disponible": set()})
        r["PET_tracer_disponible"].add(tracer)
        dates = ["NA" if pd.isna(d) else str(d.date()) for d in tg["_date"]]
        ids   = tg["_id"].tolist()
        r[f"all_{tracer}_date"] = "|".join(dates)
        r[f"nb_total_{tracer}"] = len(dates)
        r[f"{tracer}_IDs"]      = "|".join(ids)
        r[f"1st_{tracer}_ID"],   r[f"Lst_{tracer}_ID"]   = ids[0],   ids[-1]
        r[f"1st_{tracer}_date"], r[f"Lst_{tracer}_date"] = dates[0], dates[-1]
        r[f"all_{tracer}_suvr"] = _aligned(tg["_suvr"])
        if tracer.startswith("ABETA_"):
            r[f"all_{tracer}_centiloid"]  = _aligned(tg["_cl"])
            r[f"all_{tracer}_status"]     = _aligned(tg["_stat"])
            r[f"all_{tracer}_suvr_npdka"] = _aligned(tg["_npdka"])
            st = tg["_stat"].dropna()
            r[f"Lst_{tracer}_status"] = st.iloc[-1] if len(st) else np.nan
        if tracer.startswith("TAU_") and "_ento" in tg.columns:
            r[f"all_{tracer}_entorhinal_suvr"] = _aligned(tg["_ento"])
        if qc_map:
            r[f"all_{tracer}_qc"] = _aligned(tg["_qc"])

    for r in records.values():
        r["PET_tracer_disponible"] = "|".join(sorted(r["PET_tracer_disponible"])) or np.nan
    result = pd.DataFrame(list(records.values()))
    if "ID" not in result.columns:
        return pd.DataFrame()
    # scanner model/vendor from petqc (scan-level, aligned with all_{T}_date)
    return add_scanner_channels(result.set_index("ID"), pet, "NACCID", "_tracer", "_date",
                                "_mfr", "_model")


# =============================================================================
# §4  CSF BIOMARKERS
# =============================================================================
# Source: investigator_fcsf_nacc70.csv
# NACCID, CSFABETA (Aβ42), CSFPTAU (p-tau), CSFTTAU (t-tau)
# Date: CSFLPYR/CSFLPMO/CSFLPDY (lumbar puncture date)
# Measurement date per analyte: CSFABYR/CSFABMO, CSFPTYR/CSFPTMO, etc.
# =============================================================================

CSF_MARKERS = [
    ("CSFABETA", "ABETA42"),
    ("CSFPTAU",  "PTAU"),
    ("CSFTTAU",  "TAU"),
]
# Per-analyte metadata columns in investigator_fcsf (prefix CSFAB / CSFPT / CSFTT):
#   {p}MD   assay method code, {p}MDX  method free text (when MD = 8 'Other')
#   {p}MO / {p}DY / {p}YR   assay date
# Method codes per the NACC FCSF RDD (verify against your RDD version):
NACC_CSF_METHOD = {1: "ELISA", 2: "Luminex xMAP", 8: "Other"}
NACC_CSF_PREFIX = {"CSFABETA": "CSFAB", "CSFPTAU": "CSFPT", "CSFTTAU": "CSFTT"}


def _nacc_ymd(df, p):
    cols = [f"{p}YR", f"{p}MO", f"{p}DY"]
    if not all(c in df.columns for c in cols):
        return pd.Series(pd.NaT, index=df.index)
    y, m, d = (pd.to_numeric(df[c], errors="coerce") for c in cols)
    ok = y.between(1900, 2100)
    return pd.to_datetime(dict(year=y.where(ok), month=m.where(m.between(1, 12), 1),
                               day=d.where(d.between(1, 31), 1)), errors="coerce")


def build_csf(df_csf: pd.DataFrame) -> pd.DataFrame:
    """
    Per analyte (aligned lists): all_{m}_date (lumbar puncture date),
    all_{m}_score, all_{m}_rundate (assay date), {m}_kit (assay method),
    {m}_lab (ADC centre that sent the result), nb_of_{m}_censored.
    """
    _print("  §4  CSF biomarkers...")
    if df_csf.empty:
        return pd.DataFrame()

    csf = make_csf_date(df_csf.copy())
    csf = csf.sort_values(["NACCID", "_date"])

    parts = []
    for raw_col, canonical in CSF_MARKERS:
        if raw_col not in csf.columns:
            continue
        p = NACC_CSF_PREFIX[raw_col]
        t = csf[["NACCID", "_date", raw_col]].copy()
        md = pd.to_numeric(csf.get(f"{p}MD"), errors="coerce") if f"{p}MD" in csf.columns else None
        if md is not None:
            lab_ = md.map(NACC_CSF_METHOD)
            if f"{p}MDX" in csf.columns:
                txt = csf[f"{p}MDX"].astype(str).str.strip()
                txt = txt.where(~txt.isin(["", "nan", "-4", "-4.0"]))
                lab_ = lab_.where(md != 8, "Other: " + txt.fillna("unspecified"))
            t["_method"] = lab_
        t["RUNDATE"] = _nacc_ymd(csf, p)
        if "NACCADC" in csf.columns:
            t["LAB"] = "NACC ADC " + csf["NACCADC"].astype(str)
        part = plasma_series(t, "NACCID", raw_col, canonical, "_date", None,
                             method_col="_method" if md is not None else None)
        if not part.empty:
            parts.append(part)

    if not parts:
        return pd.DataFrame()

    result = reduce(lambda a, b: a.join(b, how="outer"), parts)
    result["CSF_available"] = "yes"
    return result


# =============================================================================
# §5  NEUROPSYCHOLOGICAL + VITALS  (from main file — chunked)
# =============================================================================



# UDS A5 health history + NACC-derived medication flags (0/1/2 codes;
# 9 / -4 / 8 -> NaN). SEIZURES: 0 absent, 1 recent/active, 2 remote/inactive.


# Valid codes per health variable (everything else → NaN); default (0, 1, 2)





def build_all_clinical() -> pd.DataFrame:
    """§5  NPS, vitals and health history from the main UDS file (read in chunks).

    Only the columns used by the variables configured for NACC in
    ``variables.yaml`` are read. Negative NACC codes are set to missing for
    the columns in ``NPS_NUMERIC_COLS``; other cleaning rules come from the
    variable entries (``valid_codes``, ``missing_values``, ``valid_range``).
    """
    _print("  §5  Clinical / NPS / Vitals (chunked main file)...")
    sources = variable_sources(DATASET, config=_CTX.config)
    raw_cols = []
    for _, src in sources:
        raw_cols += src.raw_columns()
    keep = ID_DATE_COLS + list(dict.fromkeys(raw_cols))
    df   = read_main_chunked(keep_cols=keep)
    if df.empty:
        return pd.DataFrame()

    df = df.sort_values(["NACCID", "_date"])

    parts = []
    for name, src in sources:
        frame, col = materialise(df, name, src)
        if col is None:
            continue
        scores_by_subj = {}
        for subj, grp in frame.groupby("NACCID"):
            scores_by_subj[subj] = list(zip(grp["_date"], grp[col]))
        part = build_nps_columns(scores_by_subj, name)
        if not part.empty:
            parts.append(part)

    del df
    gc.collect()

    if not parts:
        return pd.DataFrame()
    return reduce(lambda a, b: a.join(b, how="outer"), parts)


# =============================================================================
# §6  DEMOGRAPHICS + APOE  (from main file — chunked)
# =============================================================================

def build_demographics() -> pd.DataFrame:
    """
    §6  Demographics + APOE.
    Reads main file chunked — restricted to demog columns.
    Uses first visit per subject for static fields.
    BIRTHMO + BIRTHYR → birth_date.
    NACCAPOE code → genotype string + e4 count.
    Born year estimated: VISITYR − NACCAGE (when DOB missing).
    """
    _print("  §6  Demographics + APOE (chunked main file)...")

    keep = ID_DATE_COLS + DEMOG_COLS + ["BIRTHMO", "BIRTHYR"]
    df   = read_main_chunked(keep_cols=keep)
    if df.empty:
        return pd.DataFrame()

    df = df.sort_values(["NACCID", "_date"])
    first = df.drop_duplicates("NACCID", keep="first")

    records = []
    for _, row in first.iterrows():
        ptid = row["NACCID"]

        # Birth date from BIRTHMO + BIRTHYR
        dob = pd.NaT
        try:
            yr = int(row.get("BIRTHYR", 0) or 0)
            mo = max(1, int(row.get("BIRTHMO", 1) or 1))
            if yr > 1800:
                dob = pd.Timestamp(year=yr, month=mo, day=1)
        except Exception:
            dob = pd.NaT

        # APOE genotype from numeric code
        apoe_code = row.get("NACCAPOE", np.nan)
        try:
            ac = int(float(apoe_code))
        except Exception:
            ac = None
        apoe_geno   = APOE_CODE_MAP.get(ac, np.nan)
        apoe_e4     = APOE_E4_COUNT.get(ac, np.nan)

        # Estimated birth year (fallback when DOB not complete)
        born_yr = np.nan
        if pd.isna(dob):
            try:
                born_yr = int(float(row.get("VISITYR", 0)) -
                              float(row.get("NACCAGE", 0)))
            except Exception:
                pass

        records.append({
            "ID"              : ptid,
            "gender"          : row.get("SEX",      np.nan),
            "race"            : row.get("RACE",     np.nan),
            "hispanic"        : row.get("HISPANIC", np.nan),
            "education"       : row.get("EDUC",     np.nan),
            "marital_status"  : row.get("MARISTAT", np.nan),
            "handedness"      : row.get("HANDED",   np.nan),
            "living_situation": row.get("NACCLIVS", np.nan),
            "independence"    : row.get("INDEPEND", np.nan),
            "residence_type"  : row.get("RESIDENC", np.nan),
            "birth_date"      : dob,
            "birth_year_est"  : born_yr,
            "age_at_1st_visit": row.get("NACCAGE",  np.nan),
            "APOE_genotype"   : apoe_geno,
            "APOE_code"       : ac if ac in APOE_CODE_MAP else np.nan,   # raw NACCAPOE
            "APOE_e4_count"   : apoe_e4,
            "famhx_any_dem"   : row.get("NACCFAM",  np.nan),
            "famhx_mother_dem": row.get("NACCMOM",  np.nan),
            "famhx_father_dem": row.get("NACCDAD",  np.nan),
            "famhx_AD_mut"    : row.get("NACCFADM", np.nan),
            "famhx_FTD_mut"   : row.get("NACCFFTD", np.nan),
        })

    del df, first
    gc.collect()

    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §7  AGE AT FIRST DIAGNOSIS + BIRTH DATE (post-merge)
# =============================================================================

def add_age_at_1st_diag(master: pd.DataFrame) -> pd.DataFrame:
    """
    age_at_1st_diag = birth_date + (1st_diag_date / 365.25)
    If birth_date unavailable, estimate from birth_year_est.
    """
    if "birth_date" in master.columns and "1st_diag_date" in master.columns:
        bd = pd.to_datetime(master["birth_date"], errors="coerce")
        d1 = pd.to_datetime(master["1st_diag_date"], errors="coerce")
        try:
            days = (d1 - bd).dt.days
            master["age_at_1st_diag"] = (days / 365.25).round(2)
            bad = (master["age_at_1st_diag"] < 0) | (master["age_at_1st_diag"] > 120)
            master.loc[bad, "age_at_1st_diag"] = np.nan
        except Exception:
            master["age_at_1st_diag"] = np.nan

    # Fallback: use age_at_1st_visit
    if "age_at_1st_diag" not in master.columns or \
       master["age_at_1st_diag"].isna().all():
        master["age_at_1st_diag"] = master.get("age_at_1st_visit", np.nan)

    return master


# =============================================================================
# MAIN
# =============================================================================

def main() -> pd.DataFrame:
    _print("=" * 62)
    _print("NACC LONGITUDINAL INTEGRATION PIPELINE  v1")
    _print("=" * 62)
    _print("\nSubject ID  : NACCID → PTID")
    _print("Time axis   : VISITYR/VISITMO/VISITDAY → _date")
    _print("Main file   : read in chunks of {:,} rows".format(CHUNK_SIZE))
    _print("Diagnosis   : NACC-derived vars → UDS D1 fallback")

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load small files ──────────────────────────────────────────────────────
    _print("\nLoading scan / CSF tables...")
    df_mri      = load_small("mri")
    df_csf      = load_small("csf")
    df_pet_gaain= load_small("pet_gaain")
    df_pet_npdka= load_small("pet_npdka")
    df_pet_tau  = load_small("pet_tau")
    df_pet_fdg  = load_small("pet_fdg")
    df_pet_qc   = load_small("pet_qc")

    # ── Build domains — main file read separately per section ─────────────────
    _print("\nBuilding domain tables...")

    # §1 Diagnosis — read main file with DX cols only
    _print("  Loading main file for diagnosis...")
    dx_keep = ID_DATE_COLS + DX_FLAG_COLS
    df_dx = read_main_chunked(keep_cols=dx_keep)
    diag  = build_diagnosis(df_dx) if not df_dx.empty else pd.DataFrame()
    del df_dx
    gc.collect()

    # §2 MRI
    mri = build_mri(df_mri)
    del df_mri

    # §3 PET
    pet = build_pet({
        "pet_gaain": df_pet_gaain,
        "pet_npdka": df_pet_npdka,
        "pet_tau"  : df_pet_tau,
        "pet_fdg"  : df_pet_fdg,
        "pet_qc"   : df_pet_qc,
    })
    del df_pet_gaain, df_pet_npdka, df_pet_tau, df_pet_fdg, df_pet_qc
    gc.collect()

    # §4 CSF
    csf = build_csf(df_csf)
    del df_csf
    gc.collect()

    # §5 NPS + Vitals (reads main file internally)
    clin = build_all_clinical()

    # §6 Demographics (reads main file internally)
    demog = build_demographics()

    # ── Merge ─────────────────────────────────────────────────────────────────
    _print("\nMerging all domains on subject ID...")
    parts = [p for p in [diag, mri, pet, csf, clin, demog]
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
    if "CSF_available" not in master.columns:
        master["CSF_available"] = np.nan

    master["Plasma_available"] = "not_available_in_NACC"

    LEAD = [
        "ID",
        "gender", "race", "hispanic", "education",
        "marital_status", "handedness", "living_situation",
        "independence", "residence_type",
        "birth_date", "birth_year_est",
        "age_at_1st_visit", "age_at_1st_diag",
        "APOE_genotype", "APOE_e4_count",
        "famhx_any_dem", "famhx_mother_dem", "famhx_father_dem",
        "famhx_AD_mut", "famhx_FTD_mut",
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
        n_conv = master["nb_of_change"].gt(0).sum()
        _print(f"\nSubjects with ≥1 diagnosis change: {n_conv:,}")

    if "MRI_type_disponible" in master.columns:
        from collections import Counter
        cnt = Counter()
        for v in master["MRI_type_disponible"].dropna():
            cnt.update(str(v).split("|"))
        if cnt:
            _print("\nMRI types:")
            for t, n in cnt.most_common():
                _print(f"  {t:<30s} {n:>6,} subjects")

    if "PET_tracer_disponible" in master.columns:
        from collections import Counter
        cnt = Counter()
        for v in master["PET_tracer_disponible"].dropna():
            cnt.update(str(v).split("|"))
        if cnt:
            _print("\nPET tracers:")
            for t, n in cnt.most_common():
                _print(f"  {t:<30s} {n:>6,} subjects")

    if "CSF_available" in master.columns:
        n_csf = master["CSF_available"].eq("yes").sum()
        _print(f"\nSubjects with CSF data: {n_csf:,}")

    _print("\nTest coverage (subjects with ≥1 valid score):")
    for name in list(dict.fromkeys(n for n, _ in variable_sources(DATASET, config=_CTX.config)))[:10]:
        col = f"nb_of_{name}_visits"
        if col in master.columns:
            n = master[col].gt(0).sum()
            if n > 0:
                _print(f"  {name:<25s} {n:>6,} subjects")

    return master


# =============================================================================
