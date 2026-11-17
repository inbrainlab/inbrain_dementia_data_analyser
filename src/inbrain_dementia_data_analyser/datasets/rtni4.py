"""4RTNI LONGITUDINAL INTEGRATION PIPELINE
Applies the same domain-builder logic as ADNI v2, adapted for 4RTNI's
unique single-wide-table with repeated timepoint column suffixes.

4RTNI-SPECIFIC ARCHITECTURE
----------------------------
4RTNI ships ONE wide CSV:  4RTNI_DATA_22May2026.csv
  - ONE ROW PER SUBJECT (not one row per visit)
  - All timepoints encoded as column name suffixes:
      _0     → Baseline (TP0)
      _6MO   → 6-month follow-up
      _12MO  → 12-month follow-up
  - SUBID as subject key
  - No EXAMDATE column — each instrument has its own {INSTRUMENT}_{TP}_DATE

DIAGNOSIS CODING (DX column — static, one label per subject):
  PSP   → Progressive Supranuclear Palsy
  CBD/S → Corticobasal Degeneration / Syndrome  (coded as "CBD" in column)
  Oth   → Other atypical parkinsonian disorder

Because DX is STATIC (one value per subject, not longitudinal), the
diagnosis trajectory approach is simplified:
  - 1st_diag = Lst_diag = DX (no changes expected)
  - nb_of_visits = number of timepoints with data (estimated from SCANID)
  - no change columns (all zeros / NaN)

TIMEPOINTS:
  TP0 (baseline), 6MO (6 months), 12MO (12 months)
  Visit ordering is fixed: 0 → 6MO → 12MO
  Visit dates: each instrument records its own date column per timepoint

MRI:
  SCANID_{TP} — MRI scan date at each timepoint (non-null = MRI performed)
  LONI IDA search export (single file since 2026-11-16, found by ida_load):
  sequence type, Study Date, Image ID and scanner (vendor / model / field)
  of every series. The old per-visit files live in old_ida_search/.

PET: Not available in the 4RTNI public release

CSF:
  BIOSPECIMENS_{TP}_CSF = 1 means CSF collected
  No protein concentration values (Aβ, tau, p-tau) in the public release.
  Serum, plasma, urine collection flags are also stored.

Demographics (static, per subject):
  SUBID, SEX, AGE AT TP0 (age at baseline), EDUCATION, RACE, LATINO
  DDURATION (disease duration at enrollment), AUTOPSYDATE, AUTOPSYPX

NPS / Motor / Disease-specific (per timepoint):
  PSPRS   — PSP Rating Scale (total + 6 subscores)
            NOTE: subscore columns use DOT separator: PSPRS_0_SUBSCORE.HISTORY
  UPDRS   — Unified Parkinson's Disease Rating Scale
  SEADL   — Schwab & England Activities of Daily Living
  CGI_S   — Clinician's Global Impression - Severity (col: CGI_S_{TP})
  CGI_C   — Clinician's Global Impression - Change (some timepoints)
  MOCA    — Montreal Cognitive Assessment (with education adjustment)
  MMSE    — Mini-Mental State Examination
  CVLT    — California Verbal Learning Test
  BENSON  — Modified Rey-Osterrieth Complex Figure
            (MODREY=copy score, REY10M=10min recall, REYRECG=recognition binary)
  DSYM    — Digit Symbol Substitution
  MTRAILS — Modified Trails (time + correct)
  DSPAN   — Digit Span (forward + backward)
  UDSTRAILS — UDS Trails A and B
  VERBFLU_DWORDS — Verbal Fluency D-words
  VERBFLU_AN — Verbal Fluency Animals
  BNT     — Boston Naming Test (total + corrected)
  STROOP  — Stroop test (color-word + color)
  VOSP    — Visual Object and Space Perception
  WRAT4   — Wide Range Achievement Test 4 (reading)
  GDS     — Geriatric Depression Scale (total + 15-item)
  SNQ22   — Social Norms Questionnaire
  NPI_Q   — Neuropsychiatric Inventory-Q (12 domains, severity + distress)
  CDR     — Clinical Dementia Rating (total + box score)
  FAQ     — Functional Activities Questionnaire
  IRI     — Interpersonal Reactivity Index
  RSMS    — Revised Self-Monitoring Scale
  BIS     — Behavioral Inhibition Scale

Output:
  <root>/Tables_4TRNI/Integrated_4RTNI/4RTNI_longitudinal.csv

Author  : VRS_inbrain project
Version : 2  —  2026-11-16  (single-file IDA search: dates, Image IDs, scanner)

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("4RTNI", root_path)
"""
from __future__ import annotations

import gc
import numpy as np
import os
import pandas as pd
import re
from functools import reduce

from ..datasets.base import DatasetContext
from ..imaging.ida import (
    ida_build_acquisitions,
    ida_classify,
    ida_combine_blocks,
    ida_load,
    ida_report_unclassified,
)
from ..series import build_scored_columns, d_fmt, list_str
from ..variables import (attach_generic_variables, clean_values, load_generic_tables,
                         variable_sources)

import logging as _logging
from ..logging_utils import make_print as _make_print

log = _logging.getLogger(__name__)
_print = _make_print(log)

DATASET = "4RTNI"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
FRTNI = OUT_DIR = OUT_LONG = ""


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, FRTNI, OUT_DIR, OUT_LONG
    _CTX = ctx
    FRTNI = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write the per-subject table of 4RTNI; returns it."""
    configure(ctx)
    return main()



def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns of one measure (see series.build_scored_columns)."""
    return build_scored_columns(scores_by_subj, test_name)


def _tp_item_sum(row, instrument: str, tp: str, items) -> float:
    """Sum of ``<instrument>_<tp>_<item>`` over ``items`` (NaN when none is filled)."""
    total = 0.0
    found = False
    for suf in items:
        val = row.get(f"{instrument}_{tp}_{suf}", np.nan)
        if pd.notna(val):
            try:
                total += float(val)
                found = True
            except (ValueError, TypeError):
                pass
    return total if found else np.nan


def extract_tp_sum_series(df: pd.DataFrame, instrument: str, items, date_suffix: str,
                          test_name: str) -> pd.DataFrame:
    """Per-timepoint item sum (e.g. NPI-Q total of 12 severity items) as a scored series."""
    scores_by_subj = {}
    for ptid, row in df.set_index("PTID").iterrows():
        visits = []
        for tp in TIMEPOINTS:
            d_col = f"{instrument}_{tp}_{date_suffix}"
            date_val = pd.NaT
            if d_col in df.columns and pd.notna(row.get(d_col)):
                date_val = pd.to_datetime(row.get(d_col), errors="coerce")
            if pd.isna(date_val):
                date_val = TP_MONTHS[tp]
            visits.append((date_val, _tp_item_sum(row, instrument, tp, items)))
        scores_by_subj[ptid] = visits
    return build_nps_columns(scores_by_subj, test_name)





# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================









# IDA search: the idaSearch_*.csv of Tables_4TRNI/ is discovered by ida_load()
# (shared-helpers cell). IDA Visit label -> timepoint suffix of the wide table
IDA_VISIT_TP = {"baseline": "0", "month 6": "6MO", "month 12": "12MO"}


def _4rtni_visit_tp(v):
    """'Baseline' -> '0', 'Month 6' -> '6MO', 'Month 12' -> '12MO'."""
    return IDA_VISIT_TP.get(str(v).strip().lower(), np.nan)


def _4rtni_legacy_visit(fname: str) -> str:
    """Old per-visit file name -> Visit label of the new export."""
    m = re.search(r"_(\d+)_months", fname, re.I)
    if not m:
        return ""
    return "Baseline" if int(m.group(1)) == 0 else f"Month {int(m.group(1))}"

# Fixed timepoint labels in chronological order
TIMEPOINTS = ["0", "6MO", "12MO"]
TP_MONTHS  = {"0": 0, "6MO": 6, "12MO": 12}

# IDA Description classification (used by classify_4rtni_description)
# Infrastructure / derived series — skipped (same rules as the NIFD pipeline)
SKIP_PATTERNS_4RTNI = [
    "scout", "localizer", "field_map", "field map", "fieldmap",
    "gre_field", "asset cal", "mocoseries", "aahscout", "aascout",
    "field mapping", "3db0",
]
SKIP_SUFFIXES_4RTNI = ["_fa", "_adc", "_tracew"]

# Canonical MRI types returned by classify_4rtni_description (output order)
MRI_TYPES_4RTNI = ["MRI_T1w", "MRI_T2w", "MRI_T2w_FLAIR",
                   "MRI_DTI", "MRI_fMRI", "MRI_ASL"]

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================










def load_data() -> pd.DataFrame:
    """Multi-strategy CSV loader."""
    path = PATHS["data"]
    if not os.path.exists(path):
        _print(f"  [ERROR] Data file not found → {path}")
        return pd.DataFrame()
    STRATEGIES = [
        dict(engine="c",      on_bad_lines="error", encoding="utf-8"),
        dict(engine="c",      on_bad_lines="skip",  encoding="utf-8"),
        dict(engine="python", on_bad_lines="skip",  encoding="utf-8"),
        dict(engine="python", on_bad_lines="skip",  encoding="latin-1"),
        dict(engine="python", on_bad_lines="skip",  encoding="cp1252"),
    ]
    last_err = None
    for i, kwargs in enumerate(STRATEGIES):
        try:
            df = pd.read_csv(path, low_memory=False, **kwargs)
            if i > 0:
                _print(f"  [WARN] data: loaded with fallback (strategy {i+1})")
            _print(f"  [OK]   data: {len(df):>5,} rows  {df.shape[1]} cols")
            return df
        except Exception as e:
            last_err = e
    _print(f"  [ERROR] data: all strategies failed — {last_err}")
    return pd.DataFrame()


def normalise_id(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise subject ID column to 'PTID'. 4RTNI uses SUBID."""
    if "PTID" in df.columns:
        return df
    col_lower = {c.lower(): c for c in df.columns}
    for alt in ("subid", "ptid", "subject_id"):
        if alt in col_lower:
            df = df.rename(columns={col_lower[alt]: "PTID"})
            return df
    return df


def get_tp_col(df: pd.DataFrame, instrument: str, tp: str,
               suffix: str = "") -> str | None:
    """
    Resolve the column name for a given instrument, timepoint, and optional
    suffix. Returns None if the column does not exist.
    Example: get_tp_col(df, "MMSE", "6MO", "MMSETOT")
             → "MMSE_6MO_MMSETOT" if it exists
    """
    col = f"{instrument}_{tp}_{suffix}" if suffix else f"{instrument}_{tp}"
    # Try exact case first, then case-insensitive
    if col in df.columns:
        return col
    col_lower = col.lower()
    for c in df.columns:
        if c.lower() == col_lower:
            return c
    return None





def _make_col(instrument: str, tp: str, suffix: str) -> str:
    """
    Build the 4RTNI column name for a given instrument, timepoint, suffix.
    When suffix is empty → {instrument}_{tp}   (no trailing underscore)
    When suffix present  → {instrument}_{tp}_{suffix}
    Examples:
      SEADL, 0,   ""      → "SEADL_0"
      MMSE,  6MO, "MMSETOT" → "MMSE_6MO_MMSETOT"
    """
    return f"{instrument}_{tp}_{suffix}" if suffix else f"{instrument}_{tp}"


def _find_col(df: pd.DataFrame, instrument: str, tp: str,
              suffix: str) -> str | None:
    """Case-insensitive column lookup via _make_col()."""
    target = _make_col(instrument, tp, suffix).lower()
    for c in df.columns:
        if c.lower() == target:
            return c
    return None


def extract_tp_series(df: pd.DataFrame, instrument: str,
                      score_suffix: str, date_suffix: str,
                      test_name: str) -> pd.DataFrame:
    """
    Extract a longitudinal series from the 4RTNI wide-table timepoint columns.

    Column name rules (no trailing underscores when suffix is empty):
      score → _make_col(instrument, tp, score_suffix)
      date  → _make_col(instrument, tp, date_suffix)

    Examples with empty score_suffix:
      SEADL, "", "DATE"  → score=SEADL_0,      date=SEADL_0_DATE
      CGIS,  "", "DATE"  → score=CGIS_0,        date=CGIS_0_DATE
      DSYM,  "", "DATE"  → score=DSYM_0,        date=DSYM_0_DATE

    Examples with non-empty score_suffix:
      MMSE, "MMSETOT", "DATE" → score=MMSE_0_MMSETOT, date=MMSE_0_DATE
    """
    df = normalise_id(df)
    any_found = False
    scores_by_subj = {}

    for ptid, row in df.set_index("PTID").iterrows():
        visits = []
        for tp in TIMEPOINTS:
            d_col = _find_col(df, instrument, tp, date_suffix)
            s_col = _find_col(df, instrument, tp, score_suffix)

            if d_col is None and s_col is None:
                visits.append((TP_MONTHS[tp], np.nan))
                continue
            any_found = True

            # Parse date
            date_val = pd.NaT
            if d_col:
                raw_d = row.get(d_col, None)
                if pd.notna(raw_d):
                    date_val = pd.to_datetime(raw_d, errors="coerce")
            if pd.isna(date_val):
                date_val = TP_MONTHS[tp]   # integer month offset fallback

            # Parse score
            score_val = np.nan
            if s_col:
                raw_s = row.get(s_col, None)
                if pd.notna(raw_s):
                    try:
                        score_val = float(raw_s)
                    except (ValueError, TypeError):
                        score_val = raw_s

            visits.append((date_val, score_val))

        scores_by_subj[ptid] = visits

    if not any_found:
        s_ex = _make_col(instrument, "0", score_suffix)
        d_ex = _make_col(instrument, "0", date_suffix)
        _print(f"    [WARN] Columns not found → {test_name} skipped"
              f"  (looked for: '{s_ex}', '{d_ex}')")
        return pd.DataFrame()

    return build_nps_columns(scores_by_subj, test_name)


def classify_4rtni_description(desc: str) -> str | None:
    """Return canonical MRI type for an IDA Description string, or None to skip."""
    if not isinstance(desc, str) or not desc.strip():
        return None
    d_lo = desc.lower()
    d_up = desc.upper()

    if any(p in d_lo for p in SKIP_PATTERNS_4RTNI):
        return None
    if any(d_up.endswith(s.upper()) for s in SKIP_SUFFIXES_4RTNI):
        return None
    if any(x in d_up for x in ["FRACTIONAL ANISO", "ISOTROPIC IMAGE",
                                  "T2-WEIGHTED TRACE", "AVERAGE DC"]):
        return None

    # DTI primary acquisitions
    if any(x in d_up for x in ["DTI", "DWI", "DIFFUSION", "ADVDIFF",
                                  "B2000", "B1000", "B=2000", "B=1000"]):
        return "MRI_DTI"
    # T1w
    if any(x in d_lo for x in ["t1_mprage", "t1mprage", "mprage",
                                  "mp-rage", "fl3d_tra",
                                  "spgr"]):          # GE: Sag IR-SPGR
        return "MRI_T1w"
    # T2w FLAIR
    if any(x in d_lo for x in ["t2_flair", "t2flair", "flair"]):
        return "MRI_T2w_FLAIR"
    # T2w structural
    if any(x in d_lo for x in ["t2_spc", "t2_tse", "ax fse t2",
                                  "axial t2", "t2star", "t2 var",
                                  "t2 fse", "t2 frfse"]):   # GE: Ax T2 FSE / FRFSE
        return "MRI_T2w"
    # fMRI / resting-state
    if any(x in d_lo for x in ["resting", "rsbold", "fmri", "bold",
                                  "epi_", "ep2d_pace", "ep2d_te",
                                  "rswhole", "rsbrain", "ge_func"]):
        if "field" not in d_lo:
            return "MRI_fMRI"
    # ASL
    if any(x in d_lo for x in ["asl", "pasl", "perfusion",
                                  "relcbf", "cbf", "perfusion_weighted",
                                  "cerebral blood flow"]):
        return "MRI_ASL"

    return None   # unknown — skip


def is_infra_4rtni(desc) -> bool:
    """True for scouts / field maps / derived DTI maps (skipped on purpose)."""
    if not isinstance(desc, str) or not desc.strip():
        return True
    d_lo, d_up = desc.lower(), desc.upper()
    return (any(p in d_lo for p in SKIP_PATTERNS_4RTNI)
            or any(d_up.endswith(x.upper()) for x in SKIP_SUFFIXES_4RTNI)
            or any(x in d_up for x in ["FRACTIONAL ANISO", "ISOTROPIC IMAGE",
                                       "T2-WEIGHTED TRACE", "AVERAGE DC"]))


# =============================================================================
# §1  DIAGNOSIS
# =============================================================================
# Source: DX column (static — one label per subject, no longitudinal change)
#
# DX values:
#   PSP   → PSP (Progressive Supranuclear Palsy)
#   CBD   → CBD (Corticobasal Degeneration/Syndrome)
#   Oth   → OtherDx
#
# Because diagnosis is STATIC in 4RTNI:
#   - All trajectory fields are simplified (no changes)
#   - nb_of_visits estimated from non-null SCANID columns
#   - DDURATION (disease duration at enrollment) added as extra field
# =============================================================================

DX_MAP_4RTNI = {
    "PSP" : "PSP",
    "CBD" : "CBD",
    "OTH" : "OtherDx",
    "OTHER": "OtherDx",
}


def _classify_4rtni(raw) -> str:
    if pd.isna(raw):
        return "Unknown"
    t = str(raw).strip().upper()
    return DX_MAP_4RTNI.get(t, t)


def build_diagnosis(df: pd.DataFrame) -> pd.DataFrame:
    _print("  §1  Diagnosis (static, one label per subject)...")
    df = normalise_id(df)

    # Count timepoints with imaging data as proxy for nb_of_visits
    scanid_cols = [f"SCANID_{tp}" for tp in TIMEPOINTS if f"SCANID_{tp}" in df.columns]

    records = []
    for ptid, row in df.set_index("PTID").iterrows():
        dx_label = _classify_4rtni(row.get("DX", np.nan))

        # Estimate visit count from non-null SCANIDs
        n_visits = sum(
            1 for c in scanid_cols
            if pd.notna(row.get(c)) and str(row.get(c)).strip() not in ("", "nan")
        ) if scanid_cols else 1

        # First date: use MOCA_0_DATE or MMSE_0_DATE or PSPRS_0_DATE as baseline proxy
        baseline_date = pd.NaT
        for date_probe in ["MOCA_0_DATE", "MMSE_0_DATE", "PSPRS_0_DATE",
                           "CDR_0_DATE", "UPDRS_0_DATE"]:
            if date_probe in df.columns and pd.notna(row.get(date_probe)):
                baseline_date = pd.to_datetime(row.get(date_probe), errors="coerce")
                if pd.notna(baseline_date):
                    break

        records.append({
            "ID"                         : ptid,
            "1st_diag"                   : dx_label,
            "Lst_diag"                   : dx_label,   # static — no change
            "1st_diag_date"              : baseline_date,
            "Lst_diag_date"              : baseline_date,
            "all_diag"                   : dx_label,
            "all_diag_date"              : d_fmt(baseline_date),
            "nb_of_visits"               : n_visits,
            "1st_time_of_change"         : np.nan,
            "1st_time_of_change_days_to" : np.nan,
            "all_time_of_change"         : np.nan,
            "all_time_of_change_days_to" : np.nan,
            "Lst_time_of_change"         : np.nan,
            "nb_of_change"               : 0,
            "Diag_total"                 : dx_label,
            "disease_duration_at_enroll" : row.get("DDURATION", np.nan),
            "autopsy_date"               : row.get("AUTOPSYDATE", np.nan),
            "autopsy_pathology"          : row.get("AUTOPSYPX", np.nan),
        })

    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §2  MRI — IDA search export (primary) + SCANID_{TP} of the wide table
# =============================================================================
# IDA search (ida_load / ida_build_acquisitions, shared-helpers cell):
#   one row per image series with Description, Study Date, Image ID and
#   scanner. A SESSION = subject x type x Study Date.
# SCANID_{TP}: scan date recorded in the wide table (non-null = MRI performed).
#   Kept as raw columns, and used as T1w evidence for subjects with a SCANID
#   but no row at all in the IDA export.
#
# Columns produced per canonical type (every list has nb_total_{TYPE} tokens):
#   all_{TYPE}_date        study date of each session (YYYY-MM-DD)
#   all_{TYPE}_visit       timepoint of each session ("0|6MO|12MO")
#   {TYPE}_IDs             LONI Image ID of each session
#   all_{TYPE}_image_ids   every Image ID of the session ('/'-joined)
#   all_{TYPE}_manufacturer / _model / _field_strength, {TYPE}_manufacturer
#   nb_total_{TYPE}, nb_images_{TYPE}, 1st_/Lst_{TYPE}_visit, 1st_/Lst_{TYPE}_ID
# Plus: SCANID_0, SCANID_6MO, SCANID_12MO and MRI_type_disponible
# =============================================================================

def load_ida() -> pd.DataFrame:
    """IDA search of 4RTNI, classified (canon) with timepoint tokens."""
    ida = ida_load(FRTNI, legacy_visit=_4rtni_legacy_visit, label="4RTNI")
    if ida.empty:
        return ida
    ida["canon"]     = ida_classify(ida, classify_4rtni_description)
    ida["visit_tok"] = ida["visit"].map(_4rtni_visit_tp)
    ida["visit_ord"] = ida["visit_tok"].map(TP_MONTHS)
    return ida


def build_mri(df_wide: pd.DataFrame) -> pd.DataFrame:
    """
    §2  MRI acquisitions — IDA search export + SCANID_{TP}.
    IDA gives the sequence-level sessions (dates, Image IDs, scanner);
    SCANID_{TP} is carried as raw columns and backs up subjects that are
    absent from the IDA export (T1w session at the SCANID date).
    """
    _print("  §2  MRI acquisitions (SCANID + IDA)...")
    df_wide = normalise_id(df_wide.copy())

    # ── Part A: SCANID values per subject (wide -> long, no row loop) ─────────
    scan_cols = [f"SCANID_{tp}" for tp in TIMEPOINTS if f"SCANID_{tp}" in df_wide.columns]
    scan = df_wide.drop_duplicates("PTID").set_index("PTID")[scan_cols].astype(object)
    scan = scan.apply(lambda c: c.where(c.notna() & ~c.astype(str).str.strip()
                                        .isin(["", "nan"])).map(
        lambda v: str(v).strip() if pd.notna(v) else np.nan))
    scan.index.name = "ID"

    # ── Part B: IDA search export ────────────────────────────────────────────
    ida = load_ida()
    if not ida.empty:
        ida_report_unclassified(ida, "canon", skip_fn=is_infra_4rtni)
    img = ida_build_acquisitions(ida, "canon") if not ida.empty else pd.DataFrame()

    # ── Part C: subjects with a SCANID but no IDA row -> T1w at the SCANID date
    long = scan.stack().rename("scanid").reset_index()
    long.columns = ["subject", "col", "scanid"]
    long = long[~long["subject"].isin(img.index)] if not img.empty else long
    if not long.empty:
        long["visit_tok"] = long["col"].str.replace("SCANID_", "", regex=False)
        long["visit_ord"] = long["visit_tok"].map(TP_MONTHS)
        long["date"]      = pd.to_datetime(long["scanid"], errors="coerce")
        long["canon"]     = "MRI_T1w"
        fb = ida_build_acquisitions(long, "canon")
        _print(f"    [INFO] {fb.shape[0]} subjects with SCANID but not in the IDA "
              f"export -> MRI_T1w from SCANID dates")
        img = ida_combine_blocks(img, fb)

    result = scan.join(img, how="outer") if not img.empty else scan
    result = result.dropna(how="all")
    result.index.name = "ID"
    if "MRI_type_disponible" in result.columns:      # same column order as before
        lead = ["MRI_type_disponible"] + scan_cols
        result = result[lead + [c for c in result.columns if c not in lead]]
    return result


# =============================================================================
# §3  BIOSPECIMEN COLLECTION FLAGS
# =============================================================================
# BIOSPECIMENS_{TP}_{TYPE} = 1 means collected, 0 = not collected
# Types: SERUM, PLASMA, URINE, CSF
#
# No protein concentrations are in the 4RTNI public release.
# Columns produced (per specimen type):
#   all_{type}_collected_visits  → pipe-separated TPs where collected
#   nb_{type}_collections        → count of TPs where collected
# =============================================================================

BIO_TYPES = ["SERUM", "PLASMA", "URINE", "CSF"]


def build_biospecimens(df: pd.DataFrame) -> pd.DataFrame:
    _print("  §3  Biospecimen collection flags...")
    df = normalise_id(df)

    records = {}
    for ptid, row in df.set_index("PTID").iterrows():
        rec = {"ID": ptid}
        for btype in BIO_TYPES:
            tps_collected = []
            for tp in TIMEPOINTS:
                col = f"BIOSPECIMENS_{tp}_{btype}"
                if col in df.columns:
                    val = row.get(col, np.nan)
                    try:
                        if pd.notna(val) and int(float(val)) == 1:
                            tps_collected.append(tp)
                    except (ValueError, TypeError):
                        pass
            rec[f"all_{btype.lower()}_collected_visits"] = list_str(tps_collected)
            rec[f"nb_{btype.lower()}_collections"]       = len(tps_collected)
        records[ptid] = rec

    if not records:
        return pd.DataFrame()
    return pd.DataFrame(list(records.values())).set_index("ID")


# =============================================================================
# §4  CLINICAL / NPS / MOTOR SCALES
# =============================================================================
# All instruments follow the pattern: {INSTRUMENT}_{TP}_{SUFFIX}
# Each entry in NPS_TESTS:
#   (instrument, score_suffix, date_suffix, canonical_name)
#
# PSPRS subscores are also extracted individually.
# NPI_Q total = sum of severity scores across 12 domains.
# =============================================================================



# NPI-Q severity columns (12 domains) — summed to NPI total
# NPI-Q severity columns confirmed from dictionary (12 domains):






def build_all_clinical(df: pd.DataFrame) -> pd.DataFrame:
    """Clinical / NPS / motor series from the wide 4RTNI table.

    Each variable in ``variables.yaml`` (``sources: 4RTNI:``) names an
    ``instrument``, a score ``suffix`` and a ``date_suffix``; the columns read
    are ``<instrument>_<timepoint>_<suffix>`` for the timepoints 0, 6MO, 12MO.
    ``transform: sum_timepoint_items`` sums ``<instrument>_<tp>_<item>`` over
    ``items`` at each timepoint (NPI-Q total).
    """
    _print("  §4  Clinical / NPS / Motor scales...")
    df = normalise_id(df)
    parts = []
    for name, src in variable_sources(DATASET, config=_CTX.config):
        if src.transform == "sum_timepoint_items":
            part = extract_tp_sum_series(df, src.instrument, src.items, src.date_suffix or "DATE", name)
        else:
            frame = df
            if src.missing_values or src.valid_range or src.valid_codes is not None:
                # cleaning rules of the variable (e.g. TMT 995-998 codes)
                cols = [c for c in (_find_col(df, src.instrument, tp, src.suffix or "")
                                    for tp in TIMEPOINTS) if c]
                if cols:
                    frame = df.copy()
                    for c in cols:
                        frame[c] = clean_values(frame[c], src)
            part = extract_tp_series(frame, src.instrument, src.suffix or "", src.date_suffix or "DATE", name)
        if not part.empty:
            parts.append(part)

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()
    return reduce(lambda a, b: a.join(b, how="outer"), parts)


# =============================================================================
# §5  DEMOGRAPHICS
# =============================================================================
# All static (one value per subject):
#   SUBID    → PTID → ID
#   SEX      → gender
#   AGE AT TP0 → age_at_baseline  (NOTE: column name has spaces)
#   EDUCATION → education
#   RACE     → race
#   LATINO   → hispanic
# =============================================================================

def build_demographics(df: pd.DataFrame) -> pd.DataFrame:
    _print("  §5  Demographics...")
    df = normalise_id(df)

    # Handle "AGE AT TP0" column (spaces in name)
    # AGE column: file uses "AGE_AT_TP0" (underscores)
    # dict shows "AGE AT TP0" (spaces) — actual header confirmed underscores
    age_col = next(
        (c for c in df.columns
         if c.upper().replace(" ", "_") == "AGE_AT_TP0"),
        None
    )

    records = []
    for ptid, row in df.set_index("PTID").iterrows():
        records.append({
            "ID"              : ptid,
            "gender"          : row.get("SEX", np.nan),
            "race"            : row.get("RACE", np.nan),
            "hispanic"        : row.get("LATINO", np.nan),
            "education"       : row.get("EDUCATION", np.nan),
            "age_at_baseline" : row.get(age_col, np.nan) if age_col else np.nan,
        })

    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §6  AGE AT FIRST DIAGNOSIS
# =============================================================================
# 4RTNI provides AGE_AT_TP0 (age at baseline enrollment).
# This IS the age at first clinical assessment — used directly.
# =============================================================================

def add_age_at_1st_diag(master: pd.DataFrame) -> pd.DataFrame:
    if "age_at_baseline" in master.columns:
        master["age_at_1st_diag"] = master["age_at_baseline"]
    return master


# =============================================================================
# MAIN
# =============================================================================

def main() -> pd.DataFrame:
    _print("=" * 62)
    _print("4RTNI LONGITUDINAL INTEGRATION PIPELINE  v2")
    _print("=" * 62)
    _print("\nNOTE: 4RTNI is a STATIC-DIAGNOSIS cohort (PSP/CBD).")
    _print("      One row per subject; timepoints encoded as column suffixes.")
    _print("      No PET or CSF protein data in public release.")
    _print("      IDA search export: sequence-level MRI sessions (date, Image ID, scanner).")

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load ──────────────────────────────────────────────────────────────────
    _print("\nLoading raw table...")
    df = load_data()
    if df.empty:
        _print("ERROR: data file could not be loaded.")
        return pd.DataFrame()

    df.columns = [str(c).strip() for c in df.columns]
    df = normalise_id(df)

    _print(f"\n  Subjects       : {len(df):,}")
    _print(f"  Columns        : {df.shape[1]:,}")
    if "DX" in df.columns:
        _print(f"  DX distribution:")
        _print(df["DX"].value_counts().to_string())

    # ── Build domains ─────────────────────────────────────────────────────────
    _print("\nBuilding domain tables...")
    diag = build_diagnosis(df)
    mri  = build_mri(df)
    bio  = build_biospecimens(df)
    clin = build_all_clinical(df)
    demog = build_demographics(df)

    # ── Merge ─────────────────────────────────────────────────────────────────
    _print("\nMerging all domains on subject ID...")
    parts = [p for p in [diag, mri, bio, clin, demog]
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

    # Fill unavailable flags
    if "MRI_type_disponible" not in master.columns:
        master["MRI_type_disponible"] = np.nan
    master["PET_tracer_disponible"] = "not_available_in_4RTNI"
    master["CSF_proteins_available"] = "not_available_in_4RTNI_public"

    master["Plasma_available"] = "not_available_in_4RTNI_public (collection flags only)"

    LEAD = [
        "ID",
        "gender", "race", "hispanic", "education",
        "age_at_baseline", "age_at_1st_diag",
        "1st_diag", "Lst_diag", "1st_diag_date",
        "all_diag", "nb_of_visits",
        "nb_of_change", "Diag_total",
        "disease_duration_at_enroll", "autopsy_date", "autopsy_pathology",
        "MRI_type_disponible",
        "SCANID_0", "SCANID_6MO", "SCANID_12MO",
        "PET_tracer_disponible", "CSF_proteins_available",
        # Biospecimen collection flags
        "all_csf_collected_visits", "nb_csf_collections",
        "all_plasma_collected_visits", "nb_plasma_collections",
        "all_serum_collected_visits", "nb_serum_collections",
        "all_urine_collected_visits", "nb_urine_collections",
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
        _print("\nDiagnosis distribution:")
        _print(master["1st_diag"].value_counts().to_string())

    if "nb_total_T1w" in master.columns:
        n_mri = master["nb_total_T1w"].gt(0).sum()
        _print(f"\nSubjects with ≥1 MRI: {n_mri}")

    if "nb_csf_collections" in master.columns:
        n_csf = master["nb_csf_collections"].gt(0).sum()
        _print(f"Subjects with ≥1 CSF collection: {n_csf}")

    _print("\nTest coverage (subjects with ≥1 valid score):")
    for name in list(dict.fromkeys(n for n, _ in variable_sources(DATASET, config=_CTX.config)))[:12]:
        col = f"nb_of_{name}_visits"
        if col in master.columns:
            n = master[col].gt(0).sum()
            if n > 0:
                _print(f"  {name:<25s} {n:>5,} subjects")

    return master


# =============================================================================
