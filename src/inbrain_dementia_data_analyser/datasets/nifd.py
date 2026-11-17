"""NIFD LONGITUDINAL INTEGRATION PIPELINE  v2
MRI / PET acquisitions come from the LONI IDA search export. Since the
2026-11-16 download this is ONE file (any idaSearch_*.csv in Tables_NIFD/ that
has 'Study Date' + 'Image ID'; found automatically by ida_load) with columns:
  'Subject ID', 'Sex', 'Research Group', 'Visit', 'Study Date', 'Age',
  'Modality', 'Description', 'Type', 'Imaging Protocol', 'Image ID'
The four old per-visit files (idaSearch_6_08_2026_{0,6,12,18}_months.csv, only
'Subject ID', 'Sex', 'Age', 'Description') now live in old_ida_search/ and are
read ONLY to fill rows the new export does not have (legacy fill).

The 'Description' column encodes the sequence type; the PET tracer is read
from the Radiopharmaceutical field of 'Imaging Protocol' when present.
Each acquisition is dated with its real 'Study Date' (same calendar as
CLINICAL_LINKDATE) — the IDA 'Visit' label is kept as a token only, because
in NIFD one label can cover scans many months apart.

Canonical type mapping (derived from full description inventory):

  MRI_T1w        : T1_mprage*, t1_mprage*, MP-RAGE, fl3d_tra, Sagittal MP-Rage
  MRI_T2w_FLAIR  : T2_flair*, t2_flair*, AXIAL FLAIR, T2_FLAIR
  MRI_T2w        : T2_spc*, t2_spc*, T2 FSE, Ax FSE T2, AXIAL T2, t2star*
  MRI_DTI        : DTI*, dti_*, ep2d-advdiff*, DIFFUSION* (excludes _FA/_ADC/_TRACEW)
  MRI_fMRI       : Resting*, RSWhole*, RSBrain*, Bold*, fmri*, EPI_*, ep2d_pace*
  MRI_ASL        : ASL*, pASL*, pasl*, perfusion*, relCBF, 3DASL, Perfusion_Weighted
  MRI_T2star     : Gradient echo, fl3d_tra_p2_swi*, t2star*
  PET_PIB        : *PIB*            (or Radiopharmaceutical = PIB)
  PET_FDG        : *FDG*            (or Radiopharmaceutical = FDG)
  ABETA_AV45     : Radiopharmaceutical = 18F-AV45 (florbetapir)

SKIP (infrastructure): scouts, field maps, MoCoSeries, localizers
SKIP (derived):  _FA, _ADC, _TRACEW, _TRACEW maps from DTI

Author  : VRS_inbrain project
Version : 3  —  2026-11-16  (single-file IDA search: dates, Image IDs, scanner)

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("NIFD", root_path)
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
    ida_load,
    ida_report_unclassified,
)
from ..series import build_scored_columns, d_fmt, days_from, list_str, safe_date
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

DATASET = "NIFD"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
NIFD = OUT_DIR = OUT_LONG = ""


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, NIFD, OUT_DIR, OUT_LONG
    _CTX = ctx
    NIFD = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write the per-subject table of NIFD; returns it."""
    configure(ctx)
    return main()



def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns of one measure (see series.build_scored_columns)."""
    return build_scored_columns(scores_by_subj, test_name)




# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================









# IDA search: the idaSearch_*.csv of Tables_NIFD/ is discovered by ida_load()
# (shared-helpers cell) — no download date hard-coded here.
# PET tracer (Radiopharmaceutical field) -> column family used by this cohort
NIFD_TRACER_MAP = {"ABETA_PIB": "PET_PIB", "PET_FDG": "PET_FDG",
                   "ABETA_AV45": "ABETA_AV45"}


def _nifd_visit_month(v):
    """IDA Visit label 'Month 12' -> 12 (visit token, months as in v2)."""
    m = re.search(r"(\d+)", str(v))
    return int(m.group(1)) if m else np.nan


def _nifd_legacy_visit(fname: str) -> str:
    """Old per-visit file name -> Visit label of the new export."""
    m = re.search(r"_(\d+)_months", fname, re.I)
    return f"Month {int(m.group(1))}" if m else ""

# =============================================================================
# ── CANONICAL TYPE RULES ──────────────────────────────────────────────────────
# Rules are evaluated in ORDER — first match wins.
# Each rule: (canonical_label, list_of_match_functions)
# A description matches if ANY function returns True.
# =============================================================================

def _up(d): return d.upper()
def _lo(d): return d.lower()

CANONICAL_RULES = [
    # ── PET (check before MRI to avoid cross-contamination) ──────────────────
    ("PET_PIB",       [lambda d: "PIB" in _up(d)]),
    ("PET_FDG",       [lambda d: "FDG" in _up(d)]),

    # ── MRI: DTI (primary acquisition only — skip derived maps) ──────────────
    ("MRI_DTI", [
        lambda d: (
            any(x in _up(d) for x in
                ["DTI", "DWI", "DIFFUSION", "ADVDIFF", "B2000", "B1000",
                 "B=2000", "B=1000", "BP+_TE86_B2000"])
            and not any(x in _up(d) for x in
                        ["_FA", "_ADC", "_TRACEW", "FRACTIONAL", "ANISO",
                         "ISOTROPIC", "AVERAGE DC", "T2-WEIGHTED TRACE",
                         "AVERAGE DCM"])
        )
    ]),

    # ── MRI: T1w ─────────────────────────────────────────────────────────────
    ("MRI_T1w", [
        lambda d: any(x in _lo(d) for x in
                      ["t1_mprage", "t1mprage", "mprage", "mp-rage",
                       "fl3d_tra_iso", "sagittal mp"]),
    ]),

    # ── MRI: T2w FLAIR ───────────────────────────────────────────────────────
    ("MRI_T2w_FLAIR", [
        lambda d: any(x in _lo(d) for x in
                      ["t2_flair", "t2flair", "flair"]),
        lambda d: _up(d) in ["T2_FLAIR", "AXIAL FLAIR"],
    ]),

    # ── MRI: T2w (structural, non-FLAIR) ─────────────────────────────────────
    ("MRI_T2w", [
        lambda d: any(x in _lo(d) for x in
                      ["t2_spc", "t2_tse", "t2 var", "ax fse t2",
                       "axial t2", "t2star", "ax t2"]),
        lambda d: re.search(r"ax(ial)?\s+t2", _lo(d)) is not None,
    ]),

    # ── MRI: fMRI / resting-state ─────────────────────────────────────────────
    ("MRI_fMRI", [
        lambda d: any(x in _lo(d) for x in
                      ["resting", "rsbold", "rswhole", "rsbrainstem",
                       "fmri_run", "ax rsfmri", "bold"]),
        lambda d: any(x in _lo(d) for x in
                      ["epi_syncomp", "epi_repetition", "epi_standard",
                       "epi_pseudoread", "ep2d_pace", "ep2d_te"]),
        lambda d: re.search(r"resting.state", _lo(d)) is not None,
    ]),

    # ── MRI: ASL / perfusion ─────────────────────────────────────────────────
    ("MRI_ASL", [
        lambda d: any(x in _lo(d) for x in
                      ["asl", "pasl", "perfusion", "relcbf", "cbf",
                       "3dasl", "perfusion_weighted"]),
    ]),

    # ── MRI: T2* / SWI / GRE ────────────────────────────────────────────────
    ("MRI_T2star", [
        lambda d: any(x in _lo(d) for x in
                      ["swi", "fl3d_tra_p2_swi", "gradient echo",
                       "t2star_fl2d"]),
        lambda d: _lo(d).startswith("gradient echo"),
    ]),
]

# Infrastructure / derived — SKIP entirely
SKIP_PATTERNS = [
    "scout", "localizer", "field_map", "field map", "fieldmap",
    "gre_field", "asset cal", "mocoseries", "aahscout", "aascout",
    "field mapping", "3db0",
]
SKIP_SUFFIXES = ["_fa", "_adc", "_tracew"]


def classify_description(desc: str) -> str | None:
    """
    Return canonical label for a Description string, or None to skip.
    """
    if not isinstance(desc, str) or not desc.strip():
        return None

    d_lo = desc.lower()
    d_up = desc.upper()

    # Infrastructure skip
    if any(p in d_lo for p in SKIP_PATTERNS):
        return None

    # Derived-map skip (DTI post-processed outputs)
    if any(d_up.endswith(s.upper()) for s in SKIP_SUFFIXES):
        return None
    if any(x in d_up for x in ["FRACTIONAL ANISO", "FRACTIONAL ANSIO",
                                  "AVERAGE DC", "ISOTROPIC IMAGE",
                                  "T2-WEIGHTED TRACE"]):
        return None

    # Apply canonical rules
    for label, fns in CANONICAL_RULES:
        if any(fn(desc) for fn in fns):
            return label

    return None   # unknown — skip silently


def is_infra_description(desc) -> bool:
    """True for scouts / field maps / derived DTI maps (skipped on purpose)."""
    if not isinstance(desc, str) or not desc.strip():
        return True
    d_lo, d_up = desc.lower(), desc.upper()
    return (any(p in d_lo for p in SKIP_PATTERNS)
            or any(d_up.endswith(x.upper()) for x in SKIP_SUFFIXES)
            or any(x in d_up for x in ["FRACTIONAL ANISO", "FRACTIONAL ANSIO",
                                       "AVERAGE DC", "ISOTROPIC IMAGE",
                                       "T2-WEIGHTED TRACE"]))


# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================













def load_clinical() -> pd.DataFrame:
    path = PATHS["clinical"]
    if not os.path.exists(path):
        _print(f"  [ERROR] Clinical file not found → {path}")
        return pd.DataFrame()
    for sheet in (0, "Sheet1", "Sheet3", 1):
        try:
            df = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
            if "LONI_ID" in df.columns and len(df) > 1:
                _print(f"  [OK]   clinical (sheet={sheet!r}): "
                      f"{len(df):>5,} rows  {df.shape[1]} cols")
                return df
        except Exception:
            continue
    _print(f"  [ERROR] Cannot read clinical data → {path}")
    return pd.DataFrame()


def load_ida() -> pd.DataFrame:
    """
    IDA search of NIFD, one row per image series, classified:
      canon      canonical type (Description rules; PET tracer from the
                 Radiopharmaceutical field when the export has it)
      visit_tok  month of the IDA Visit label (0 / 6 / 12 / 18)
    Requires the shared-helpers cell (ida_load, ida_classify …).
    """
    ida = ida_load(NIFD, legacy_visit=_nifd_legacy_visit, label="NIFD")
    if ida.empty:
        return ida
    ida["canon"] = ida_classify(ida, classify_description)
    tracer = ida["tracer"].map(NIFD_TRACER_MAP)
    use_tr = ida["modality"].eq("PET") & tracer.notna()
    n_new  = int((use_tr & ida["canon"].isna()).sum())
    ida["canon"] = ida["canon"].where(~use_tr, tracer)
    if n_new:
        _print(f"    [INFO] {n_new:,} PET rows typed by Radiopharmaceutical "
              f"(description alone did not say the tracer)")
    ida["visit_tok"] = ida["visit"].map(_nifd_visit_month)
    ida["visit_ord"] = ida["visit_tok"]
    return ida





def _nps_from_wide(df, test_name, score_col,
                   date_col="CLINICAL_LINKDATE") -> pd.DataFrame:
    if score_col not in df.columns:
        _print(f"    [WARN] '{score_col}' not found → {test_name} skipped")
        return pd.DataFrame()
    actual_date = date_col if date_col in df.columns else "CLINICAL_LINKDATE"
    sub = df[["LONI_ID", actual_date, score_col]].copy()
    safe_date(sub, actual_date)
    sub = sub.sort_values(["LONI_ID", actual_date])
    scores_by_subj = {}
    for subj, grp in sub.groupby("LONI_ID"):
        scores_by_subj[subj] = list(zip(grp[actual_date], grp[score_col]))
    return build_nps_columns(scores_by_subj, test_name)


# =============================================================================
# §1  DIAGNOSIS TRAJECTORY
# =============================================================================

DX_MAP = {
    "CON"  : "CON",
    "BV"   : "bvFTD",
    "PNFA" : "nfvPPA",
    "SV"   : "svPPA",
    "OTH"  : "OtherDem",
    "OTHER": "OtherDem",
}


def _classify_nifd(raw_dx) -> str:
    if pd.isna(raw_dx):
        return "Unknown"
    parts = [p.strip().upper() for p in str(raw_dx).split("+")]
    labels = sorted({DX_MAP.get(p, p) for p in parts if p})
    return "+".join(labels) if labels else "Unknown"


def build_diagnosis(df: pd.DataFrame) -> pd.DataFrame:
    _print("  §1  Diagnosis trajectory...")
    safe_date(df, "CLINICAL_LINKDATE")
    dx = df[["LONI_ID", "CLINICAL_LINKDATE", "VISIT_NUMBER", "DX"]].copy()
    dx["DX_LABEL"] = dx["DX"].apply(_classify_nifd)
    dx = dx.sort_values(["LONI_ID", "CLINICAL_LINKDATE", "VISIT_NUMBER"])
    records = []
    for ptid, grp in dx.groupby("LONI_ID"):
        grp    = grp.sort_values(["CLINICAL_LINKDATE", "VISIT_NUMBER"])
        labels = grp["DX_LABEL"].tolist()
        dates  = grp["CLINICAL_LINKDATE"].tolist()
        first_date = dates[0]
        changes_dates, changes_days = [], []
        prev = labels[0]
        for lbl, dt in zip(labels[1:], dates[1:]):
            if lbl != prev:
                changes_dates.append(dt)
                changes_days.append(days_from(first_date, dt))
                prev = lbl
        records.append({
            "ID"                         : ptid,
            "1st_diag"                   : labels[0],
            "Lst_diag"                   : labels[-1],
            "1st_diag_date"              : dates[0],
            "Lst_diag_date"              : dates[-1],
            "all_diag"                   : list_str(labels),
            "all_diag_date"              : list_str([d_fmt(d) for d in dates]),
            "nb_of_visits"               : len(grp),
            "1st_time_of_change"         : changes_dates[0]  if changes_dates else np.nan,
            "1st_time_of_change_days_to" : changes_days[0]   if changes_days  else np.nan,
            "all_time_of_change"         : list_str([d_fmt(d) for d in changes_dates]),
            "all_time_of_change_days_to" : list_str(changes_days),
            "Lst_time_of_change"         : changes_dates[-1] if changes_dates else np.nan,
            "nb_of_change"               : len(changes_dates),
            "Diag_total"                 : list_str(sorted(set(labels))),
        })
    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §2  MRI & PET — from the IDA search export
# =============================================================================
#
# Logic (ida_build_acquisitions, shared-helpers cell):
#   1. Each row's Description (+ PET tracer) -> canonical label (or None)
#   2. A SESSION = subject x type x Study Date; one token per session
#
# Per-modality columns produced (every list has nb_total_{TYPE} tokens):
#   all_{TYPE}_date          — study date of each session (YYYY-MM-DD)
#   all_{TYPE}_visit         — IDA visit month of each session (e.g. "0|6|18")
#   {TYPE}_IDs               — LONI Image ID of each session
#   all_{TYPE}_image_ids     — every Image ID of the session ('/'-joined)
#   all_{TYPE}_manufacturer / _model / _field_strength — scanner per session
#   {TYPE}_manufacturer      — vendor(s) across sessions
#   nb_total_{TYPE}          — number of sessions
#   nb_images_{TYPE}         — number of image series
#   1st_/Lst_{TYPE}_visit, 1st_/Lst_{TYPE}_ID
#
# Summary columns:
#   MRI_type_disponible    — sorted pipe-joined set of MRI types present
#   PET_tracer_disponible  — sorted pipe-joined set of PET types present
#   IDA_legacy_types       — types with sessions taken from the old per-visit
#                            files (date inferred from the visit, no Image ID)
# =============================================================================

MRI_TYPES = ["MRI_T1w", "MRI_T2w_FLAIR", "MRI_T2w",
             "MRI_DTI", "MRI_fMRI", "MRI_ASL", "MRI_T2star"]
PET_TYPES = ["PET_PIB", "PET_FDG", "ABETA_AV45"]
ALL_TYPES = MRI_TYPES + PET_TYPES


def build_imaging(ida: pd.DataFrame) -> pd.DataFrame:
    """
    §2  MRI + PET acquisitions from the IDA search export.
    Returns one row per subject indexed by LONI_ID (normalised to 'ID').
    """
    _print("  §2  MRI & PET acquisitions (IDA)...")

    if ida is None or ida.empty:
        _print("    [WARN] No IDA file loaded — imaging section skipped.")
        return pd.DataFrame()

    ida_report_unclassified(ida, "canon", skip_fn=is_infra_description)
    result = ida_build_acquisitions(ida, "canon")
    if result.empty:
        _print("    [WARN] No usable IDA rows after classification — imaging skipped.")
    return result


# =============================================================================
# §3  NEUROPSYCHOLOGICAL + FUNCTIONAL + BEHAVIORAL + MOTOR
# =============================================================================






def build_all_clinical(df: pd.DataFrame) -> pd.DataFrame:
    """Clinical / NPS / functional / motor series from the clinical spreadsheet.

    Variables, source columns and their date columns come from
    ``variables.yaml`` (``sources: NIFD:``); default date ``CLINICAL_LINKDATE``.
    """
    _print("  §3  Clinical / NPS / Functional / Motor...")
    parts = []
    for name, src in variable_sources(DATASET, config=_CTX.config):
        frame, col = materialise(df, name, src)
        if col is None:
            continue
        part = _nps_from_wide(frame, name, col, src.date_column or "CLINICAL_LINKDATE")
        if not part.empty:
            parts.append(part)

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()
    return reduce(lambda a, b: a.join(b, how="outer"), parts)


# =============================================================================
# §4  DEMOGRAPHICS
# =============================================================================

def build_demographics(df: pd.DataFrame) -> pd.DataFrame:
    _print("  §4  Demographics...")
    safe_date(df, "CLINICAL_LINKDATE")
    safe_date(df, "DOB")
    df_sorted = df.sort_values(["LONI_ID", "CLINICAL_LINKDATE", "VISIT_NUMBER"])
    first = df_sorted.drop_duplicates("LONI_ID", keep="first")
    records = []
    for _, row in first.iterrows():
        records.append({
            "ID"        : row["LONI_ID"],
            "gender"    : row.get("GENDER",    np.nan),
            "race"      : row.get("RACE",      np.nan),
            "education" : row.get("EDUCATION", np.nan),
            "birth_date": row.get("DOB",       pd.NaT),
        })
    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §5  AGE AT FIRST DIAGNOSIS
# =============================================================================

def add_age_at_1st_diag(master: pd.DataFrame) -> pd.DataFrame:
    if "1st_diag_date" not in master.columns or "birth_date" not in master.columns:
        return master
    bd = pd.to_datetime(master["birth_date"],    errors="coerce")
    d1 = pd.to_datetime(master["1st_diag_date"], errors="coerce")
    try:
        days = (d1 - bd).dt.days
        master["age_at_1st_diag"] = (days / 365.25).round(2)
        bad = (master["age_at_1st_diag"] < 0) | (master["age_at_1st_diag"] > 120)
        master.loc[bad, "age_at_1st_diag"] = np.nan
    except Exception:
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
    _print("NIFD LONGITUDINAL INTEGRATION PIPELINE  v3")
    _print("=" * 62)

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load ──────────────────────────────────────────────────────────────────
    _print("\nLoading raw tables...")
    df       = load_clinical()
    ida      = load_ida()

    if df.empty:
        _print("ERROR: clinical file could not be loaded.")
        return pd.DataFrame()

    df.columns = [str(c).strip() for c in df.columns]
    safe_date(df, "CLINICAL_LINKDATE")
    safe_date(df, "DOB")

    _print(f"\n  Unique subjects (clinical): {df['LONI_ID'].nunique():,}")
    _print(f"  Visit range: "
          f"{df['VISIT_NUMBER'].min()} – {df['VISIT_NUMBER'].max()}")
    _print(f"  DX distribution (unique subjects):")
    _print(df.drop_duplicates('LONI_ID')['DX']
            .value_counts()
            .to_string(index=True))

    # ── Build domains ─────────────────────────────────────────────────────────
    _print("\nBuilding domain tables...")
    diag    = build_diagnosis(df)
    imaging = build_imaging(ida)
    clin    = build_all_clinical(df)
    demog   = build_demographics(df)

    # ── Merge ─────────────────────────────────────────────────────────────────
    _print("\nMerging all domains on subject ID...")
    parts = [p for p in [diag, imaging, clin, demog]
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

    # Fill imaging flags for subjects with no IDA data
    if "MRI_type_disponible" not in master.columns:
        master["MRI_type_disponible"]  = np.nan
    if "PET_tracer_disponible" not in master.columns:
        master["PET_tracer_disponible"] = np.nan

    # No CSF in NIFD public release
    master["CSF_available"] = "not_available_in_NIFD_public"

    master["Plasma_available"] = "not_available_in_NIFD"

    LEAD = [
        "ID",
        "gender", "race", "education", "birth_date", "age_at_1st_diag",
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

    if "MRI_type_disponible" in master.columns:
        from collections import Counter
        _print("\nMRI types available:")
        cnt = Counter()
        for v in master["MRI_type_disponible"].dropna():
            cnt.update(str(v).split("|"))
        for t, n in cnt.most_common():
            _print(f"  {t:<30s} {n:>5,} subjects")

    if "PET_tracer_disponible" in master.columns:
        from collections import Counter
        _print("\nPET tracers available:")
        cnt = Counter()
        for v in master["PET_tracer_disponible"].dropna():
            cnt.update(str(v).split("|"))
        for t, n in cnt.most_common():
            _print(f"  {t:<30s} {n:>5,} subjects")

    if "nb_of_visits" in master.columns:
        _print(f"\nVisit counts: mean={master['nb_of_visits'].mean():.1f}  "
              f"max={master['nb_of_visits'].max()}  "
              f"median={master['nb_of_visits'].median():.0f}")

    return master


# =============================================================================
