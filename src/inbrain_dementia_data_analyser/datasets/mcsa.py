"""MCSA LONGITUDINAL INTEGRATION PIPELINE  v2
MRI and PET acquisitions come from the LONI IDA search export.

IDA SEARCH (single file since the 2026-11-16 download — any idaSearch_*.csv in
Tables_MCSA/ that has 'Study Date' + 'Image ID'; found by ida_load):
  Columns: 'Subject ID', 'Sex', 'Weight', 'Visit', 'Study Date', 'Age',
           'Modality', 'Description', 'Type', 'Imaging Protocol', 'Image ID'
  The 17 old files (idaSearch_6_09_2026 (N).csv, N = visit number, only
  'Subject ID', 'Sex', 'Age', 'Description', 'Type') now live in
  old_ida_search/ and are read only to fill rows missing from the new export.

DESCRIPTION → CANONICAL TYPE (full inventory analysis):

  MRI_T1w       : Sag 3D MP-RAGE, Accel Sagittal 3D MPRAGE*, Sagittal MPRAGE*
  MRI_T2w_FLAIR : AXIAL FLAIR*, Sagittal 3D FLAIR, 2nd AXIAL FLAIR
  MRI_DTI       : AX DTI*, Axial MB DTI, AX Flair DTI 128
  PET_PIB       : *PIB* (TR:BRAIN 3D:PIB, AX:BRAIN:PIB, TRANSAXIAL:BRAIN 3D:PIB)

SKIP: none — every description maps to a canonical type or is unknown.

SESSION RULE: one session per subject x type x Study Date.
  The IDA 'Visit' (= visit_num of the main CSV) is kept as a token.

CRITICAL MCSA NOTES:
  - Subject ID in the IDA export is MCSA_ID (matching mcsa_id in main CSV)
  - Dates are de-identified per subject (first visit = 1950-01-01), but the
    IDA 'Study Date' is in the SAME shifted calendar as VisitDate/ImagingDate
    of the main CSV (it equals ImagingDate for ~94% of the sessions), so
    intervals are real and the acquisitions are dated with it.
  - In the IDA export Sex is 'X' and Age is 0 for everybody (de-identified)
  - PIB in IDA = same tracer as PIB_SUVR in main CSV → confirms PET_PIB
  - No FDG, no tau PET, no AV45 in MCSA IDA

ARCHITECTURE CHANGE v1→v2:
  §2 MRI now REPLACES the MRIMANUFACTURER/HPVOL presence flag approach
     with IDA-derived sequence-level data (T1w, T2w_FLAIR, DTI)
  §3 PET now SUPPLEMENTS PIB_SUVR presence with IDA confirmation
  Both sections merged: IDA primary, flag-based fallback for unmatched subjects.

Output:
  <root>/Tables_MCSA/Integrated_MCSA/MCSA_longitudinal.csv

Author  : VRS_inbrain project
Version : 3  —  2026-11-16  (single-file IDA search: dates, Image IDs, scanner)

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("MCSA", root_path)
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
from ..imaging.scanner import norm_vendor
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

DATASET = "MCSA"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
MCSA = OUT_DIR = OUT_LONG = ""


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, MCSA, OUT_DIR, OUT_LONG
    _CTX = ctx
    MCSA = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write the per-subject table of MCSA; returns it."""
    configure(ctx)
    return main()



def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns of one measure (see series.build_scored_columns)."""
    return build_scored_columns(scores_by_subj, test_name)




# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================









# IDA search: the idaSearch_*.csv of Tables_MCSA/ is discovered by ida_load()
# (shared-helpers cell) — no download date / file numbers hard-coded here.
def _mcsa_legacy_visit(fname: str) -> str:
    """Old file 'idaSearch_6_09_2026 (N).csv' -> Visit label 'N' of the new export."""
    m = re.search(r"\((\d+)\)", fname)
    return m.group(1) if m else ""

MRI_TYPES_MCSA = ["MRI_T1w", "MRI_T2w_FLAIR", "MRI_DTI"]
PET_TYPES_MCSA = ["PET_PIB"]

# =============================================================================
# ── CANONICAL TYPE RULES ──────────────────────────────────────────────────────
# =============================================================================

def classify_mcsa_description(desc: str) -> str | None:
    """
    Map a Description string to a canonical type, or None to skip.
    Based on full MCSA IDA description inventory.
    """
    if not isinstance(desc, str) or not desc.strip():
        return None
    dl = desc.lower().strip()

    # PET — PIB tracer (check before MRI to avoid any cross-match)
    if "pib" in dl:
        return "PET_PIB"

    # MRI T1w
    if any(x in dl for x in ["mprage", "mp-rage", "mp rage"]):
        return "MRI_T1w"

    # MRI FLAIR (exclude DTI variants)
    if "flair" in dl and "dti" not in dl:
        return "MRI_T2w_FLAIR"

    # MRI DTI
    if "dti" in dl:
        return "MRI_DTI"

    return None   # unknown — skip silently


# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================













def load_data() -> pd.DataFrame:
    path = PATHS["data"]
    if not os.path.exists(path):
        _print(f"  [ERROR] Data file not found → {path}")
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
                _print(f"  [WARN] data: loaded via fallback strategy {i+1}")
            _print(f"  [OK]   data: {len(df):>6,} rows  {df.shape[1]} cols")
            return df
        except Exception:
            continue
    _print("  [ERROR] data: all strategies failed")
    return pd.DataFrame()


def load_ida() -> pd.DataFrame:
    """IDA search of MCSA, classified (canon) with visit-number tokens.
    Requires the shared-helpers cell (ida_load, ida_classify …)."""
    ida = ida_load(MCSA, legacy_visit=_mcsa_legacy_visit, label="MCSA")
    if ida.empty:
        return ida
    ida["canon"]     = ida_classify(ida, classify_mcsa_description)
    ida["visit_tok"] = pd.to_numeric(ida["visit"], errors="coerce").astype("Int64")
    ida["visit_ord"] = ida["visit_tok"]
    return ida


def normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise MCSA column names: strip whitespace, MCSA_ID → mcsa_id."""
    df.columns = [c.strip() for c in df.columns]
    col_lower = {c.lower(): c for c in df.columns}
    if "mcsa_id" in col_lower and col_lower["mcsa_id"] != "mcsa_id":
        df = df.rename(columns={col_lower["mcsa_id"]: "mcsa_id"})
    if "visit_num" not in df.columns:
        for alt in ("VISIT_NUM", "visitnum", "VisitNum", "Visit_Num"):
            if alt in df.columns:
                df = df.rename(columns={alt: "visit_num"})
                break
    return df





def _nps_from_wide(df: pd.DataFrame, test_name: str,
                   score_col: str, date_col: str = "VisitDate") -> pd.DataFrame:
    df = normalise_columns(df)
    if score_col not in df.columns:
        _print(f"    [WARN] '{score_col}' not found → {test_name} skipped")
        return pd.DataFrame()
    actual_date = date_col if date_col in df.columns else "VisitDate"
    sub = df[["mcsa_id", actual_date, score_col]].copy()
    safe_date(sub, actual_date)
    sub = sub.sort_values(["mcsa_id", actual_date])
    scores_by_subj = {}
    for subj, grp in sub.groupby("mcsa_id"):
        scores_by_subj[subj] = list(zip(grp[actual_date], grp[score_col]))
    return build_nps_columns(scores_by_subj, test_name)


# =============================================================================
# §1  DIAGNOSIS TRAJECTORY (unchanged from v1)
# =============================================================================

FCOGIMPR_MAP = {0: "CN", 1: "MCI", 3: "DEMENTIA"}


def _classify_mcsa(val) -> str:
    if pd.isna(val):
        return "Unknown"
    try:
        return FCOGIMPR_MAP.get(int(val), "Unknown")
    except (ValueError, TypeError):
        return "Unknown"


def build_diagnosis(df: pd.DataFrame) -> pd.DataFrame:
    _print("  §1  Diagnosis trajectory...")
    df = normalise_columns(df)
    safe_date(df, "VisitDate")
    dx = df[["mcsa_id", "visit_num", "VisitDate", "FCOGIMPR"]].copy()
    dx["DX_LABEL"] = dx["FCOGIMPR"].apply(_classify_mcsa)
    dx = dx.sort_values(["mcsa_id", "visit_num"])
    records = []
    for subj, grp in dx.groupby("mcsa_id"):
        grp       = grp.sort_values("visit_num")
        labels    = grp["DX_LABEL"].tolist()
        dates     = grp["VisitDate"].tolist()
        first_date = dates[0]
        changes_dates, changes_days = [], []
        prev = labels[0]
        for lbl, dt in zip(labels[1:], dates[1:]):
            if lbl != prev:
                changes_dates.append(dt)
                changes_days.append(days_from(first_date, dt))
                prev = lbl
        records.append({
            "ID"                         : subj,
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
# §2 + §3  MRI & PET — IDA search export (primary) + main CSV flags (fallback)
# =============================================================================
# PRIMARY — IDA search export:
#   Subject ID: 'Subject ID' column (= mcsa_id)
#   A SESSION = subject x type x Study Date (ida_build_acquisitions)
#
# FALLBACK — main CSV:
#   MRI: MRIMANUFACTURER not-null OR HPVOL not-null → MRI_T1w
#   PET: PIB_SUVR not-null → PET_PIB
#   Applied only for subjects NOT found in the IDA export; sessions are dated
#   with ImagingDate (same de-identified calendar as the IDA Study Date).
#
# Columns produced per canonical type (every list has nb_total_{TYPE} tokens):
#   MRI_type_disponible, PET_tracer_disponible
#   all_{TYPE}_date        study date of each session (de-identified calendar)
#   all_{TYPE}_visit       visit number of each session
#   {TYPE}_IDs             LONI Image ID of each session        (IDA only)
#   all_{TYPE}_image_ids   every Image ID of the session        (IDA only)
#   all_{TYPE}_manufacturer / _model / _field_strength, {TYPE}_manufacturer
#   nb_total_{TYPE}, nb_images_{TYPE}, 1st_/Lst_{TYPE}_visit, 1st_/Lst_{TYPE}_ID
# =============================================================================

def build_imaging_ida(ida: pd.DataFrame, df_main: pd.DataFrame = None) -> pd.DataFrame:
    """MRI + PET sessions from the IDA search export. When the main CSV is
    given, PIB_SUVR / PIB_CENTILOID of the same subject + visit are written per
    PET session (all_PET_PIB_suvr / all_PET_PIB_centiloid, aligned with the dates)."""
    if ida is None or ida.empty:
        return pd.DataFrame()
    ida_report_unclassified(ida, "canon")
    if df_main is not None and {"PIB_SUVR", "mcsa_id", "visit_num"}.issubset(df_main.columns):
        q = df_main[["mcsa_id", "visit_num"] + [c for c in ("PIB_SUVR", "PIB_CENTILOID")
                                                 if c in df_main.columns]].copy()
        q["mcsa_id"] = q["mcsa_id"].astype(str)
        q = q.dropna(subset=["PIB_SUVR"]).drop_duplicates(["mcsa_id", "visit_num"]) \
             .set_index(["mcsa_id", "visit_num"])
        key = pd.MultiIndex.from_arrays([ida["subject"],
                                         pd.to_numeric(ida["visit"], errors="coerce")])
        is_pib = ida["canon"].eq("PET_PIB").values
        ida["suvr"] = np.where(is_pib, q["PIB_SUVR"].reindex(key).values, np.nan)
        if "PIB_CENTILOID" in q.columns:
            ida["centiloid"] = np.where(is_pib, q["PIB_CENTILOID"].reindex(key).values, np.nan)
    return ida_build_acquisitions(ida, "canon")


def build_imaging_flags(df: pd.DataFrame,
                         ida_subjects: set) -> pd.DataFrame:
    """
    Fallback: MRI + PET presence from main CSV flags, for subjects
    NOT already in the IDA export. Same column schema as the IDA builder.
    """
    df = normalise_columns(df)
    safe_date(df, "ImagingDate")

    # Restrict to subjects not covered by IDA
    df = df[~df["mcsa_id"].astype(str).isin(ida_subjects)]
    if df.empty:
        return pd.DataFrame()

    base = pd.DataFrame({
        "subject"  : df["mcsa_id"].astype(str),
        "date"     : df["ImagingDate"] if "ImagingDate" in df.columns else pd.NaT,
        "visit_tok": pd.to_numeric(df["visit_num"], errors="coerce").astype("Int64"),
    })
    base["visit_ord"] = base["visit_tok"]
    parts = []

    # MRI flag
    has_mfr = "MRIMANUFACTURER" in df.columns
    mri_mask = df["MRIMANUFACTURER"].notna() if has_mfr else pd.Series(False, index=df.index)
    if "HPVOL" in df.columns:
        mri_mask = mri_mask | df["HPVOL"].notna()
    if mri_mask.any():
        m = base[mri_mask].assign(canon="MRI_T1w")
        if has_mfr:
            m["vendor"] = df.loc[mri_mask, "MRIMANUFACTURER"].map(norm_vendor)
        parts.append(m)

    # PET flag
    pet_mask = df["PIB_SUVR"].notna() if "PIB_SUVR" in df.columns \
               else pd.Series(False, index=df.index)
    if pet_mask.any():
        p = base[pet_mask].assign(canon="PET_PIB")
        p["suvr"] = pd.to_numeric(df.loc[pet_mask, "PIB_SUVR"], errors="coerce")
        if "PIB_CENTILOID" in df.columns:
            p["centiloid"] = pd.to_numeric(df.loc[pet_mask, "PIB_CENTILOID"], errors="coerce")
        parts.append(p)

    if not parts:
        return pd.DataFrame()
    result = ida_build_acquisitions(pd.concat(parts, ignore_index=True), "canon")

    # extras kept from v2 (subject-level lists from the main CSV)
    order = df.sort_values(["mcsa_id", "visit_num"])
    ids = order["mcsa_id"].astype(str)
    if has_mfr:
        result["MRI_manufacturer"] = (order[mri_mask.reindex(order.index)]
            .groupby(ids)["MRIMANUFACTURER"]
            .agg(lambda x: list_str(x.dropna().unique().tolist())))
    if "PIB_CENTILOID" in df.columns:
        result["all_PIB_CENTILOID"] = (order[pet_mask.reindex(order.index)]
            .groupby(ids)["PIB_CENTILOID"]
            .agg(lambda x: list_str(x.dropna().tolist())))
    return result


# =============================================================================
# §4  CLINICAL SCALES (unchanged from v1)
# =============================================================================




def build_all_clinical(df: pd.DataFrame) -> pd.DataFrame:
    """Clinical / NPS / functional / vitals / imaging-metric series.

    Variables come from ``variables.yaml`` (``sources: MCSA:``); default date
    ``VisitDate`` (imaging metrics use ``ImagingDate``).
    """
    _print("  §4  Clinical / NPS / Functional / Vitals / Imaging metrics...")
    df = normalise_columns(df)
    safe_date(df, "VisitDate")
    safe_date(df, "ImagingDate")
    parts = []
    for name, src in variable_sources(DATASET, config=_CTX.config):
        frame, col = materialise(df, name, src)
        if col is None:
            continue
        part = _nps_from_wide(frame, name, col, src.date_column or "VisitDate")
        if not part.empty:
            parts.append(part)

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()
    return reduce(lambda a, b: a.join(b, how="outer"), parts)


# =============================================================================
# §5  DEMOGRAPHICS (unchanged from v1)
# =============================================================================

def build_demographics(df: pd.DataFrame) -> pd.DataFrame:
    _print("  §5  Demographics...")
    df = normalise_columns(df)
    safe_date(df, "VisitDate")
    first = df.sort_values(["mcsa_id", "visit_num"]).drop_duplicates(
        "mcsa_id", keep="first")
    records = []
    for _, row in first.iterrows():
        records.append({
            "ID"               : row["mcsa_id"],
            "gender"           : row.get("Male",        np.nan),
            "race"             : row.get("race_bin",    np.nan),
            "hispanic"         : row.get("ETHNICTY",    np.nan),
            "education"        : row.get("EDUC",        np.nan),
            "APOE_e4_any"      : row.get("Any_E4",      np.nan),
            "age_at_1st_visit" : row.get("calc_age_vis",np.nan),
        })
    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


def add_age_at_1st_diag(master: pd.DataFrame,
                         df_raw: pd.DataFrame) -> pd.DataFrame:
    if "age_at_1st_visit" in master.columns:
        master["age_at_1st_diag"] = master["age_at_1st_visit"]
        return master
    df_raw = normalise_columns(df_raw)
    age_map = (
        df_raw.sort_values(["mcsa_id", "visit_num"])
        .drop_duplicates("mcsa_id", keep="first")
        .set_index("mcsa_id")["calc_age_vis"].to_dict()
    )
    master["age_at_1st_diag"] = master["ID"].map(age_map)
    return master


# =============================================================================
# MAIN
# =============================================================================

def main() -> pd.DataFrame:
    _print("=" * 62)
    _print("MCSA LONGITUDINAL INTEGRATION PIPELINE  v3")
    _print("=" * 62)
    _print("\nMRI+PET: IDA search export (primary) + main CSV flags (fallback).")
    _print("Dates de-identified per subject (same shifted calendar in IDA and main CSV).")

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load main data ────────────────────────────────────────────────────────
    _print("\nLoading main data table...")
    df = load_data()
    if df.empty:
        _print("ERROR: main data file not found.")
        return pd.DataFrame()
    df = normalise_columns(df)
    _print(f"  ID column confirmed: {[c for c in df.columns if c.lower()=='mcsa_id']}")
    safe_date(df, "VisitDate")
    safe_date(df, "ImagingDate")

    # ── Load IDA search export ───────────────────────────────────────────────
    _print("\nLoading IDA search export...")
    ida = load_ida()

    # ── Build domains ─────────────────────────────────────────────────────────
    _print("\nBuilding domain tables...")
    diag = build_diagnosis(df)

    # Imaging: IDA primary, main CSV flags fallback
    imaging_ida = build_imaging_ida(ida, df)
    ida_subjects = set(imaging_ida.index.astype(str)) if not imaging_ida.empty else set()
    imaging_flags = build_imaging_flags(df, ida_subjects)
    if not imaging_flags.empty:
        _print(f"    [INFO] {len(imaging_flags):,} subjects not in the IDA export -> "
              f"MRI/PET presence from main CSV flags (dated with ImagingDate)")

    # Merge IDA + flags (whole subject x type blocks, IDA takes priority)
    imaging = ida_combine_blocks(imaging_ida, imaging_flags)

    clin  = build_all_clinical(df)
    demog = build_demographics(df)

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
    master = add_age_at_1st_diag(master, df)

    if "MRI_type_disponible" not in master.columns:
        master["MRI_type_disponible"] = np.nan
    if "PET_tracer_disponible" not in master.columns:
        master["PET_tracer_disponible"] = np.nan
    master["CSF_available"] = "not_available_in_MCSA"

    master["Plasma_available"] = "not_available_in_MCSA"

    LEAD = [
        "ID",
        "gender", "race", "hispanic", "education",
        "APOE_e4_any", "age_at_1st_visit", "age_at_1st_diag",
        "1st_diag", "Lst_diag", "1st_diag_date", "Lst_diag_date",
        "all_diag", "all_diag_date", "nb_of_visits",
        "1st_time_of_change", "1st_time_of_change_days_to",
        "all_time_of_change", "all_time_of_change_days_to",
        "Lst_time_of_change", "nb_of_change", "Diag_total",
        "MRI_type_disponible", "MRI_manufacturer",
        "PET_tracer_disponible", "CSF_available",
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

    if "MRI_type_disponible" in master.columns:
        from collections import Counter
        cnt = Counter()
        for v in master["MRI_type_disponible"].dropna():
            cnt.update(str(v).split("|"))
        if cnt:
            _print("\nMRI types (IDA + flags):")
            for t, n in cnt.most_common():
                _print(f"  {t:<30s} {n:>5,} subjects")

    if "PET_tracer_disponible" in master.columns:
        from collections import Counter
        cnt = Counter()
        for v in master["PET_tracer_disponible"].dropna():
            cnt.update(str(v).split("|"))
        if cnt:
            _print("\nPET tracers (IDA + flags):")
            for t, n in cnt.most_common():
                _print(f"  {t:<30s} {n:>5,} subjects")

    _print("\nTest coverage (subjects with ≥1 valid score):")
    for canonical in list(dict.fromkeys(n for n, _ in variable_sources(DATASET, config=_CTX.config)))[:12]:
        col = f"nb_of_{canonical}_visits"
        if col in master.columns:
            n = master[col].gt(0).sum()
            if n > 0:
                _print(f"  {canonical:<25s} {n:>5,} subjects")

    return master


# =============================================================================
