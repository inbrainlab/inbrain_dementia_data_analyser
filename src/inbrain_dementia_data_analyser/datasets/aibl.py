"""AIBL LONGITUDINAL INTEGRATION PIPELINE
Applies the same domain-builder logic as the ADNI v2 pipeline, adapted for
AIBL's specific table structure (RID-based IDs, VISCODE visit ordering,
separate PET tables per tracer, split APOE allele columns, etc.).

One row per subject summarising the full longitudinal record across:
  §1  Diagnosis trajectory   (aibl_pdxconv)
  §2  MRI acquisitions       (IDA search export; aibl_mrimeta + aibl_mri3meta
                              as fallback)
  §3  PET acquisitions       (IDA search export; aibl_pibmeta + aibl_av45meta
                              + aibl_flutemeta + aibl_navmeta as fallback)
  §4  Neuropsychological     (aibl_neurobat + aibl_mmse + aibl_cdr)
  §5  CSF biomarkers         (aibl_labdata)
  §6  Demographics + APOE    (aibl_ptdemog + aibl_apoeres)

Output:
  <root>/Tables_AIBL/Integrated_AIBL/AIBL_longitudinal.csv

Author  : VRS_inbrain project
Version : 2  —  2026-11-16  (IDA search export: real Image IDs, dates, scanner)

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("AIBL", root_path)
"""
from __future__ import annotations

import gc
import numpy as np
import os
import pandas as pd
import re
from collections import Counter
from functools import reduce

from ..datasets.base import DatasetContext
from ..imaging.ida import (
    ida_build_acquisitions,
    ida_classify,
    ida_combine_blocks,
    ida_load,
    ida_match_index,
    ida_report_unclassified,
)
from ..imaging.scanner import add_scanner_channels
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

DATASET = "AIBL"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
AIBL = OUT_DIR = OUT_LONG = ""


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, AIBL, OUT_DIR, OUT_LONG
    _CTX = ctx
    AIBL = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write the per-subject table of AIBL; returns it."""
    configure(ctx)
    return main()



def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns of one measure (see series.build_scored_columns)."""
    return build_scored_columns(scores_by_subj, test_name, missing_date_as_text=True)
"""
AIBL-SPECIFIC NOTES
-------------------
Subject ID  : RID  (integer roster ID — normalised to "PTID" internally)
Visit key   : VISCODE  (integer: 0=baseline, 1=m18, 2=m36, 3=m54 …)
              Date anchor: EXAMDATE when present; otherwise the visit-to-
              approximate-date lookup below is used.
Diagnosis   : aibl_pdxconv — NO EXAMDATE column. Date must be joined from
              another table (aibl_mmse used as anchor).
IDA search  : idaSearch_*.csv in Tables_AIBL/ (LONI IDA export downloaded
              2026-11-16; found by ida_load): one row per image series with
              Visit, Study Date, Description, Imaging Protocol (vendor,
              model, field strength, radiopharmaceutical) and Image ID.
              PRIMARY source of §2 / §3: real LONI Image IDs and scans up to
              2021 (the *meta tables stop in 2018). Subject ID = RID.
              Subjects that exist only in the IDA export (no row in the 2018
              clinical tables) are kept, with imaging columns only.
MRI         : IDA types MRI_T1w (3 T), MRI_T1w_1.5T, MRI_T2w, MRI_T2w_FLAIR,
              MRI_T2star (SWI). Fallback tables, used only for a subject x
              type missing from IDA:
                aibl_mrimeta   — 1.5T acquisitions (MMSMPRAGE = T1w)
                aibl_mri3meta  — 3T acquisitions   (MMSMPRAGE = T1w)
              Both have EXAMDATE + MMCONDCT (1=conducted) + MMSMPRAGE flag.
              Their image IDs are synthetic: RID_VISCODE_EXAMDATE.
PET         : IDA tracer from the Radiopharmaceutical field. Fallback: four
              separate tables, one per tracer:
                aibl_pibmeta   → ABETA_PIB      (11C-PIB)
                aibl_av45meta  → ABETA_AV45     (18F-AV45)
                aibl_flutemeta → ABETA_Flutemetamol  (18F-Flutemetamol)
                aibl_navmeta   → ABETA_NAV4694  (18F-NAV4694)
              All share the same column schema (PBCONDCT, EXAMDATE, …).
              PBCONDCT=1 means scan was conducted.
NPS         : aibl_neurobat   — LIMMTOTAL (LM Immediate), LDELTOTAL (LM Delayed)
              aibl_mmse       — MMSCORE (MMSE total)
              aibl_cdr        — CDGLOBAL (CDR global)
CSF         : aibl_labdata    — lab chemistry; no standard CSF proteins
              (ABETA/TAU/PTAU not present in AIBL public release)
Demographics: aibl_ptdemog    — RID, PTGENDER (1=M 2=F), PTDOB (YYYY-MM-DD)
APOE        : aibl_apoeres    — APGEN1, APGEN2 (allele integers, e.g. 3,4)
"""




# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================









# IDA search: the idaSearch_*.csv of Tables_AIBL/ is discovered by ida_load()
# (shared-helpers cell) — no download date hard-coded here.
# True  -> 1.5 T T1w goes to its own family MRI_T1w_1.5T and MRI_T1w holds the
#          3 T scans only (same meaning as before: only T1w_3T fed MRI_T1w).
# False -> every T1w is MRI_T1w (the field is in all_MRI_T1w_field_strength).
AIBL_SPLIT_T1W_15T = True

# VISCODE → approximate calendar offset (months from baseline)
# Used when EXAMDATE is missing to order visits.
VISCODE_MONTHS = {0: 0, 1: 18, 2: 36, 3: 54, 4: 72, 5: 90, 6: 108}

# =============================================================================
# ── HELPERS (identical to ADNI v2 pattern) ────────────────────────────────────
# =============================================================================










def load(key: str) -> pd.DataFrame:
    """
    Multi-strategy CSV loader.  Tries progressively more permissive options:
      1. C engine, strict
      2. C engine, skip bad lines
      3. Python engine, skip bad lines, UTF-8
      4. Python engine, skip bad lines, latin-1   ← handles non-UTF-8 chars
      5. Python engine, skip bad lines, cp1252    ← Windows-1252 encoding
      6. Python engine, sep=None                  ← auto-detect separator
    Prints [WARN] if a fallback was needed, [ERROR] if all strategies fail.
    """
    path = PATHS.get(key, "")
    if not path or not os.path.exists(path):
        _print(f"  [SKIP] {key}: not found → {path}")
        return pd.DataFrame()

    STRATEGIES = [
        dict(engine="c",      on_bad_lines="error", encoding="utf-8"),
        dict(engine="c",      on_bad_lines="skip",  encoding="utf-8"),
        dict(engine="python", on_bad_lines="skip",  encoding="utf-8"),
        dict(engine="python", on_bad_lines="skip",  encoding="latin-1"),
        dict(engine="python", on_bad_lines="skip",  encoding="cp1252"),
        dict(engine="python", on_bad_lines="skip",  encoding="latin-1", sep=None),
    ]

    last_err = None
    for i, kwargs in enumerate(STRATEGIES):
        try:
            df = pd.read_csv(path, low_memory=False, **kwargs)
            if i > 0:
                enc = kwargs.get("encoding", "?")
                sep = "auto-sep" if kwargs.get("sep") is None else ""
                _print(f"  [WARN] {key}: loaded with fallback "
                      f"(strategy {i+1}: {enc} {sep})")
            _print(f"  [OK]   {key}: {len(df):>6,} rows  {df.shape[1]} cols")
            return df
        except Exception as e:
            last_err = e
            continue

    _print(f"  [ERROR] {key}: all parse strategies failed → {path}")
    _print(f"          Last error: {last_err}")
    return pd.DataFrame()


def normalise_id(df: pd.DataFrame) -> pd.DataFrame:
    """Rename RID (or subject_id) → PTID for consistency with shared builders."""
    if "PTID" in df.columns:
        return df
    col_lower = {c.lower(): c for c in df.columns}
    for alt in ("rid", "ptid", "subject_id"):
        if alt in col_lower:
            df = df.rename(columns={col_lower[alt]: "PTID"})
            return df
    return df


def sort_by_viscode_or_date(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort a per-visit DataFrame. Prefer EXAMDATE; fall back to VISCODE integer.
    Adds _sort_key column for stable ordering.
    """
    if "EXAMDATE" in df.columns:
        safe_date(df, "EXAMDATE")
        df["_sort_key"] = df["EXAMDATE"]
    elif "VISCODE" in df.columns:
        df["_sort_key"] = pd.to_numeric(df["VISCODE"], errors="coerce")
    else:
        df["_sort_key"] = range(len(df))
    return df.sort_values(["PTID", "_sort_key"])


# =============================================================================
# §1  DIAGNOSIS TRAJECTORY
# =============================================================================
# Source: aibl_pdxconv_01-Jun-2018.csv
#   Columns: RID, SITEID, VISCODE, DXCURREN, DXNORM, DXNODEP, DXMCI, DXAD,
#            DXPARK, DXOTHDEM
#   NOTE: NO EXAMDATE in this table.
#
# Diagnosis coding (same boolean-flag scheme as ADNI-1 / legacy):
#   DXNORM=1   → CN
#   DXMCI=1    → MCI
#   DXAD=1     → AD
#   DXPARK=1   → PD  (Parkinson's — rare, kept separate)
#   DXOTHDEM=1 → OtherDem
#
# DXCURREN (current diagnosis integer) is also available:
#   1=CN  2=MCI  3=AD  4=PD  5=OtherDem
#   Used as the primary signal when non-null; boolean flags as fallback.
#
# Date strategy: join EXAMDATE from aibl_mmse on (RID, VISCODE).
# =============================================================================

def _classify_aibl(row: pd.Series) -> str:
    """
    Map one aibl_pdxconv row to a compound diagnosis label.

    Design principle — collect ALL active flags, join with '+':
      DXNORM=1    → CN
      DXMCI=1     → MCI
      DXAD=1      → AD
      DXPARK=1    → PD
      DXOTHDEM=1  → OtherDem

    DXCURREN (integer summary) is used ONLY as fallback when ALL boolean
    flags are missing (e.g. older data exports), never as a priority gate
    that would suppress co-occurring flags.

    Examples:
      DXMCI=1, DXAD=1               → "AD+MCI"
      DXNORM=1                      → "CN"
      DXCURREN=2 (flags absent)     → "MCI"
      DXMCI=1, DXPARK=1             → "MCI+PD"
    """
    FLAG_MAP = [
        ("DXAD",     "AD"),
        ("DXMCI",    "MCI"),
        ("DXNORM",   "CN"),
        ("DXPARK",   "PD"),
        ("DXOTHDEM", "OtherDem"),
    ]

    # Collect every active boolean flag
    active = [label for col, label in FLAG_MAP
              if row.get(col, 0) == 1]

    if active:
        return "+".join(sorted(active))   # sorted → canonical order

    # Fallback: DXCURREN integer (when flags absent / legacy export)
    dc = row.get("DXCURREN", np.nan)
    if pd.notna(dc):
        dc = int(dc)
        mapping = {1: "CN", 2: "MCI", 3: "AD", 4: "PD", 5: "OtherDem"}
        return mapping.get(dc, "Unknown")

    return "Unknown"


def build_diagnosis(df_dx: pd.DataFrame,
                    df_date_src: pd.DataFrame) -> pd.DataFrame:
    """
    §1  Diagnosis trajectory.
    df_date_src is aibl_mmse, used to supply EXAMDATE keyed on (RID, VISCODE).
    """
    _print("  §1  Diagnosis trajectory...")
    if df_dx.empty:
        return pd.DataFrame()

    dx = normalise_id(df_dx.copy())

    # ── Attach dates ─────────────────────────────────────────────────────────
    if not df_date_src.empty:
        date_src = normalise_id(df_date_src[["RID", "VISCODE", "EXAMDATE"]].copy()
                                if "RID" in df_date_src.columns
                                else df_date_src[["PTID", "VISCODE", "EXAMDATE"]].copy())
        date_src = date_src.rename(columns={"RID": "PTID"}) if "RID" in date_src.columns else date_src
        date_src = date_src.drop_duplicates(["PTID", "VISCODE"])
        safe_date(date_src, "EXAMDATE")
        dx = dx.merge(date_src[["PTID", "VISCODE", "EXAMDATE"]],
                      on=["PTID", "VISCODE"], how="left")

    # If still no EXAMDATE, synthesise from VISCODE month offset
    if "EXAMDATE" not in dx.columns:
        dx["EXAMDATE"] = pd.NaT
    dx["_approx_order"] = pd.to_numeric(dx.get("VISCODE", 0), errors="coerce").fillna(0)

    dx["DX_LABEL"] = dx.apply(_classify_aibl, axis=1)

    records = []
    for ptid, grp in dx.sort_values(["PTID", "_approx_order"]).groupby("PTID"):
        grp = grp.sort_values("_approx_order")
        labels = grp["DX_LABEL"].tolist()
        dates  = grp["EXAMDATE"].tolist()
        first_date = dates[0]

        changes_dates, changes_days = [], []
        prev = labels[0]
        for lbl, dt in zip(labels[1:], dates[1:]):
            if lbl != prev:
                changes_dates.append(dt)
                changes_days.append(days_from(first_date, dt))
                prev = lbl

        d_fmt = lambda d: str(d.date()) if pd.notna(d) and hasattr(d, "date") else np.nan

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

    df = pd.DataFrame(records)
    return df.set_index("ID") if not df.empty else df


# =============================================================================
# §2  MRI ACQUISITIONS
# =============================================================================
# Two source tables:
#   aibl_mrimeta_01-Jun-2018.csv   (1.5T)
#   aibl_mri3meta_01-Jun-2018.csv  (3T)
#   Columns: RID, SITEID, VISCODE, MMCONDCT, EXAMDATE,
#            MMSCOUT, MMSMPRAGE, MMRMPRAGE, MMB1HEAD, MMB1BODY, MMECHO
#
# MMCONDCT=1: MRI session conducted (filter before processing)
# MMSMPRAGE=1: MPRAGE (T1w structural) was acquired
# MMRMPRAGE=1: Repeat MPRAGE acquired
#
# No explicit series_type column — sequence type inferred from flag columns:
#   MMSMPRAGE=1  → T1w
#   MMECHO=1     → T2w (echo sequence)
#
# Image ID: synthetic = "RID_VISCODE_DATE"  (no LONI image ID in AIBL public)
# Both tables are concatenated and deduplicated.
# =============================================================================

def build_mri(df_15t: pd.DataFrame, df_3t: pd.DataFrame) -> pd.DataFrame:
    _print("  §2  MRI acquisitions...")

    parts = []
    for label, df in [("1.5T", df_15t), ("3T", df_3t)]:
        if df.empty:
            continue
        df = normalise_id(df.copy())
        safe_date(df, "EXAMDATE")
        # Keep only conducted sessions
        if "MMCONDCT" in df.columns:
            df = df[df["MMCONDCT"] == 1]
        if df.empty:
            continue
        df["_field_strength"] = label
        parts.append(df)

    if not parts:
        return pd.DataFrame()

    mri = pd.concat(parts, ignore_index=True)

    # Infer sequence types from flag columns
    records = {}
    avail_key = "MRI_type_disponible"

    for ptid, grp in mri.sort_values("EXAMDATE").groupby("PTID"):
        row = records.setdefault(ptid, {"ID": ptid, avail_key: set()})

        # Build synthetic image ID
        def synth_id(r):
            date_str = str(r["EXAMDATE"].date()) if pd.notna(r.get("EXAMDATE")) else "nodate"
            return f"{ptid}_{r.get('VISCODE', 'v?')}_{date_str}"

        grp = grp.copy()
        grp["_img_id"] = grp.apply(synth_id, axis=1)

        for _, r in grp.iterrows():
            date = r.get("EXAMDATE", pd.NaT)
            img  = r["_img_id"]
            field = r.get("_field_strength", "")

            # Assign canonical type per flag
            types_this_row = []
            if r.get("MMSMPRAGE", 0) == 1:
                types_this_row.append(f"T1w_{field}" if field else "T1w")
            if r.get("MMRMPRAGE", 0) == 1:
                types_this_row.append(f"T1w_{field}_repeat" if field else "T1w_repeat")
            if r.get("MMECHO", 0) == 1:
                types_this_row.append("T2w")
            if not types_this_row:
                types_this_row = [f"MRI_{field}" if field else "MRI_unknown"]

            for mtype in types_this_row:
                row[avail_key].add(mtype)
                date_str = str(date.date()) if pd.notna(date) and hasattr(date, "date") else np.nan

                # Accumulate lists
                row.setdefault(f"all_{mtype}_date_list", [])
                row.setdefault(f"all_{mtype}_id_list",   [])
                if date_str and pd.notna(date_str):
                    row[f"all_{mtype}_date_list"].append(date_str)
                row[f"all_{mtype}_id_list"].append(img)

    # Collapse list buffers into final columns
    for ptid, row in records.items():
        row[avail_key] = "|".join(sorted(row[avail_key]))
        for mtype in list(row[avail_key].split("|")):
            dl = row.pop(f"all_{mtype}_date_list", [])
            il = row.pop(f"all_{mtype}_id_list",   [])
            uniq_dates = list(dict.fromkeys(dl))
            row[f"all_{mtype}_date"]  = list_str(uniq_dates)
            row[f"nb_total_{mtype}"]  = len(uniq_dates)
            row[f"{mtype}_IDs"]       = list_str(il)
            row[f"1st_{mtype}_ID"]    = il[0]  if il else np.nan
            row[f"Lst_{mtype}_ID"]    = il[-1] if il else np.nan

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(list(records.values()))
    return result.set_index("ID") if "ID" in result.columns else pd.DataFrame()


# =============================================================================
# §3  PET ACQUISITIONS
# =============================================================================
# Four separate tables, one per radiotracer.  All share the same schema:
#   RID, SITEID, VISCODE, PBCONDCT, EXAMDATE, PBGE, PBGEMODEL, PBSIEMENS,
#   PBSIEMMODEL, PBPHILLIPS, PBPHMODEL, PBPIBTIME, PBPIBDOS, PBINJDOS,
#   PBINJTIME, PBSCTIME, PBSTATDYN, PBNMFRAM1, PBDURATN1, ... PBTRSNFR
#
# PBCONDCT=1: scan was conducted.
# Image ID: synthetic = "RID_VISCODE_TRACER_DATE"
# Scanner manufacturer: derived from PBGE / PBSIEMENS / PBPHILLIPS flags.
# =============================================================================

PET_TABLES = {
    "pet_pib"   : "ABETA_PIB",
    "pet_av45"  : "ABETA_AV45",
    "pet_flute" : "ABETA_Flutemetamol",
    "pet_nav"   : "ABETA_NAV4694",
}


def _pet_manufacturer(row: pd.Series) -> str:
    """Infer scanner manufacturer from boolean flags."""
    if row.get("PBGE",       0) == 1: return "GE"
    if row.get("PBSIEMENS",  0) == 1: return "Siemens"
    if row.get("PBPHILLIPS", 0) == 1: return "Philips"
    # Try model string columns as fallback
    for mc in ("PBGEMODEL", "PBSIEMODEL", "PBSIEMMODEL", "PBPHMODEL"):
        val = row.get(mc, None)
        if pd.notna(val) and str(val).strip() not in ("", "nan", "0"):
            return str(val).strip()
    return "Unknown"


def build_pet(dfs: dict) -> pd.DataFrame:
    """
    §3  PET acquisitions.
    dfs: the full loaded-files dict from main(); keys like 'pet_pib', 'pet_av45'…
    """
    _print("  §3  PET acquisitions...")
    available_tracers = {k: v for k, v in PET_TABLES.items()
                         if not dfs.get(k, pd.DataFrame()).empty}
    missing_tracers   = {k: v for k, v in PET_TABLES.items()
                         if dfs.get(k, pd.DataFrame()).empty}
    if missing_tracers:
        _print(f"    [INFO] Tracers skipped (file empty/unreadable): "
              f"{list(missing_tracers.values())}")
    records = {}

    for key, tracer_label in PET_TABLES.items():
        df = dfs.get(key, pd.DataFrame())
        if df.empty:
            continue

        df = normalise_id(df.copy())
        safe_date(df, "EXAMDATE")

        # Keep conducted scans only
        if "PBCONDCT" in df.columns:
            df = df[df["PBCONDCT"] == 1]
        if df.empty:
            continue

        for ptid, grp in df.sort_values("EXAMDATE").groupby("PTID"):
            row = records.setdefault(ptid, {"ID": ptid,
                                            "PET_tracer_disponible": set()})
            row["PET_tracer_disponible"].add(tracer_label)

            dates = grp["EXAMDATE"].dropna().tolist()
            uniq_dates = list(dict.fromkeys(
                str(d.date()) if hasattr(d, "date") else str(d) for d in dates
            ))

            # Synthetic image IDs
            ids = []
            for _, r in grp.iterrows():
                ds = str(r["EXAMDATE"].date()) if pd.notna(r.get("EXAMDATE")) \
                     and hasattr(r["EXAMDATE"], "date") else "nodate"
                ids.append(f"{ptid}_{r.get('VISCODE', 'v?')}_{tracer_label}_{ds}")

            row[f"all_{tracer_label}_date"]  = list_str(uniq_dates)
            row[f"nb_total_{tracer_label}"]  = len(uniq_dates)
            row[f"{tracer_label}_IDs"]       = list_str(ids)
            row[f"1st_{tracer_label}_ID"]    = ids[0]  if ids else np.nan
            row[f"Lst_{tracer_label}_ID"]    = ids[-1] if ids else np.nan

            # Scanner metadata from first row
            first = grp.iloc[0]
            row[f"{tracer_label}_manufacturer"] = _pet_manufacturer(first)
            if "PBINJDOS" in grp.columns:
                row[f"{tracer_label}_injected_dose"] = list_str(
                    grp["PBINJDOS"].dropna().tolist())

    for r in records.values():
        r["PET_tracer_disponible"] = "|".join(sorted(r["PET_tracer_disponible"]))

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(list(records.values()))
    if "ID" not in result.columns:
        return pd.DataFrame()
    result = result.set_index("ID")

    # manufacturer per scan (flags), aligned with all_{tracer}_date
    sc = []
    for key, tracer_label in PET_TABLES.items():
        df = dfs.get(key, pd.DataFrame())
        if df.empty:
            continue
        df = normalise_id(df.copy())
        safe_date(df, "EXAMDATE")
        if "PBCONDCT" in df.columns:
            df = df[df["PBCONDCT"] == 1]
        if df.empty:
            continue
        df["_T"] = tracer_label
        df["_mfr"] = df.apply(_pet_manufacturer, axis=1)
        model = None
        for mc in ("PBGEMODEL", "PBSIEMODEL", "PBSIEMMODEL", "PBPHMODEL"):
            if mc in df.columns:
                v = df[mc].where(df[mc].astype(str).str.strip().ne("") & df[mc].notna())
                model = v if model is None else model.fillna(v)
        df["_model"] = model if model is not None else np.nan
        sc.append(df[["PTID", "_T", "EXAMDATE", "_mfr", "_model"]])
    if sc:
        result = add_scanner_channels(result, pd.concat(sc, ignore_index=True),
                                      "PTID", "_T", "EXAMDATE", "_mfr", "_model")
    return result


# =============================================================================
# §2b / §3b  MRI + PET FROM THE IDA SEARCH EXPORT  (added 2026-11-16)
# =============================================================================
# One row per image series. A SESSION = subject x type x Study Date.
# Columns per canonical type (every list has nb_total_{TYPE} tokens):
#   all_{TYPE}_date        study date of each session (YYYY-MM-DD)
#   all_{TYPE}_visit       IDA visit of each session (bl, m18, m36 …)
#   {TYPE}_IDs             LONI Image ID of each session
#   all_{TYPE}_image_ids   every Image ID of the session ('/'-joined)
#   all_{TYPE}_manufacturer / _model / _field_strength, {TYPE}_manufacturer
#   nb_total_{TYPE}, nb_images_{TYPE}, 1st_/Lst_{TYPE}_visit, 1st_/Lst_{TYPE}_ID
# MRI types : MRI_T1w, MRI_T1w_1.5T, MRI_T2w, MRI_T2w_FLAIR, MRI_T2star
# PET types : ABETA_PIB, ABETA_AV45, ABETA_Flutemetamol, ABETA_NAV4694
# =============================================================================

AIBL_PET_DESC_RULES = [          # only when the Radiopharmaceutical field is empty
    ("pib", "ABETA_PIB"), ("av45", "ABETA_AV45"), ("fbp", "ABETA_AV45"),
    ("flute", "ABETA_Flutemetamol"), ("nav", "ABETA_NAV4694"),
]

# names written by the *meta-table builders -> canonical families
AIBL_META_TYPE_MAP = {           # longest first
    "T1w_3T_repeat"  : "MRI_T1w_repeat",
    "T1w_1.5T_repeat": "MRI_T1w_1.5T_repeat",
    "T1w_3T"         : "MRI_T1w",
    "T1w_1.5T"       : "MRI_T1w_1.5T",
    "T2w"            : "MRI_T2w",
}


def classify_aibl_description(desc) -> str | None:
    """IDA Description (MRI rows) -> canonical MRI type, or None to skip."""
    if not isinstance(desc, str) or not desc.strip():
        return None
    dl = desc.lower()
    if any(x in dl for x in ("localis", "localiz", "scout")):
        return None
    if any(x in dl for x in ("mpr", "spgr")):            # MPRAGE …, t1_mpr_…, 3D_SPGR
        return "MRI_T1w"
    if "flair" in dl:
        return "MRI_T2w_FLAIR"
    if any(x in dl for x in ("swi", "t2_star", "t2star")):
        return "MRI_T2star"
    if "t2" in dl or "tse" in dl:                        # Axial PD-T2 TSE
        return "MRI_T2w"
    return None


def _aibl_pet_from_description(desc) -> str | None:
    dl = str(desc).lower()
    for key, canon in AIBL_PET_DESC_RULES:
        if key in dl:
            return canon
    return None


def _aibl_visit_code(v) -> str | float:
    """'Baseline' -> 'bl' ; '18 Month follow-up' -> 'm18' (AIBL VISCODE style)."""
    t = str(v).strip().lower()
    if t.startswith("baseline"):
        return "bl"
    m = re.search(r"(\d+)\s*month", t)
    return f"m{int(m.group(1))}" if m else np.nan


def canonise_meta_types(df: pd.DataFrame) -> pd.DataFrame:
    """Rename the families of build_mri() (T1w_3T, T1w_1.5T, T2w …) to the
    canonical names used by the IDA block (MRI_T1w, MRI_T1w_1.5T, MRI_T2w)."""
    if df is None or df.empty:
        return df
    ren = {}
    for col in df.columns:
        for old, new in AIBL_META_TYPE_MAP.items():
            rx = rf"(^|_){re.escape(old)}(_|$)"
            if re.search(rx, col):
                ren[col] = re.sub(rx, lambda m: m.group(1) + new + m.group(2), col, count=1)
                break
    df = df.rename(columns=ren)
    if "MRI_type_disponible" in df.columns:
        df["MRI_type_disponible"] = df["MRI_type_disponible"].map(
            lambda v: "|".join(sorted(AIBL_META_TYPE_MAP.get(t, t) for t in v.split("|")))
            if isinstance(v, str) else v)
    return df


def build_imaging_ida(ref_ids=None) -> pd.DataFrame:
    """
    §2b/§3b  MRI + PET sessions from the IDA search export (index cast to the
    RID dtype of the clinical tables). Empty frame when there is no export.
    Requires the shared-helpers cell (ida_load, ida_classify …).
    """
    _print("  §2b MRI + PET acquisitions (IDA search export)...")
    ida = ida_load(AIBL, label="AIBL")
    if ida.empty:
        return pd.DataFrame()

    is_pet = ida["modality"].eq("PET")
    mri_canon = ida_classify(ida, classify_aibl_description)
    pet_canon = ida["tracer"].where(ida["tracer"].notna(),
                                    ida_classify(ida, _aibl_pet_from_description))
    ida["canon"] = mri_canon.where(~is_pet, pet_canon)
    # a PET tracer never labels an MRI row and vice-versa
    ida.loc[~is_pet & ~ida["canon"].astype(str).str.startswith("MRI_"), "canon"] = np.nan
    ida.loc[is_pet & ida["canon"].astype(str).str.startswith("MRI_"), "canon"] = np.nan
    if AIBL_SPLIT_T1W_15T:
        low = ida["canon"].eq("MRI_T1w") & ida["field_strength"].eq("1.5T")
        ida.loc[low, "canon"] = "MRI_T1w_1.5T"

    ida["visit_tok"] = ida["visit"].map(_aibl_visit_code)
    ida["visit_ord"] = pd.to_numeric(ida["visit_tok"].astype(str).str.extract(r"(\d+)")[0],
                                     errors="coerce").fillna(0)
    ida_report_unclassified(ida, "canon",
                            skip_fn=lambda d: any(x in str(d).lower()
                                                  for x in ("localis", "localiz", "scout")))
    result = ida_build_acquisitions(ida, "canon")
    if not result.empty and ref_ids is not None:
        result = ida_match_index(result, ref_ids)
    return result


# =============================================================================
# §4  NEUROPSYCHOLOGICAL + CDR
# =============================================================================
# aibl_neurobat: RID, SITEID, VISCODE, EXAMDATE, LIMMTOTAL, LDELTOTAL
# aibl_mmse:     RID, SITEID, VISCODE, EXAMDATE, MMSCORE
# aibl_cdr:      RID, SITEID, VISCODE, EXAMDATE, CDGLOBAL
#
# Uses the same build_nps_columns pattern from ADNI v2.
# =============================================================================




def _nps_from_df(df: pd.DataFrame, test_name: str, score_col: str) -> pd.DataFrame:
    if df.empty or score_col not in df.columns:
        _print(f"    [WARN] '{score_col}' not found → {test_name} skipped")
        return pd.DataFrame()
    df = normalise_id(df.copy())
    safe_date(df, "EXAMDATE")
    df = df.sort_values(["PTID", "EXAMDATE"])
    scores_by_subj = {}
    for subj, grp in df.groupby("PTID"):
        scores_by_subj[subj] = list(zip(grp["EXAMDATE"], grp[score_col]))
    return build_nps_columns(scores_by_subj, test_name)


def build_all_clinical(dfs: dict) -> pd.DataFrame:
    """§4  Neuropsychological + CDR series.

    Variables and source columns come from ``variables.yaml`` (``sources: AIBL:``);
    each series is dated with the table's ``EXAMDATE``.
    """
    _print("  §4  Clinical / NPS...")
    parts = []
    for name, src in variable_sources(DATASET, config=_CTX.config):
        df = dfs.get(src.table, pd.DataFrame())
        if df.empty:
            continue
        frame, col = materialise(df, name, src)
        if col is None:
            continue
        parts.append(_nps_from_df(frame, name, col))

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()
    return reduce(lambda a, b: a.join(b, how="outer"), parts)


# =============================================================================
# §5  LAB DATA (proxy for bio-related metrics — no standard CSF proteins)
# =============================================================================
# aibl_labdata: RID, SITEID, VISCODE, AXT117 (TSH), BAT126, HMT3, HMT7,
#               HMT13, HMT40, HMT100, HMT102, RCT6, RCT11, RCT20, RCT392
#
# AIBL's public release does NOT include ABETA/TAU/PTAU CSF concentrations.
# Labdata contains blood/serum markers. We collect what's available as
# longitudinal summaries using the same pattern.
#
# Stored as all_{col}_values | nb_of_{col}_values per marker.
# Glucose = RCT11 (AIBL reports SI units, mmol/L). The merge notebook renames
# all_RCT11_values → all_Lab_Glucose_AIBL_score (kept apart from ADNI mg/dL).
# =============================================================================

LAB_COLS = ["AXT117", "BAT126", "HMT3", "HMT7", "HMT13",
            "HMT40", "HMT100", "HMT102", "RCT6", "RCT11", "RCT20", "RCT392"]


AIBL_LAB_NAMES = {   # code → analyte (ADNI/AIBL LABDATA dictionary)
    "AXT117": "TSH", "BAT126": "Vitamin B12", "HMT3": "RBC", "HMT7": "WBC",
    "HMT13": "Platelets", "HMT40": "Hemoglobin", "HMT100": "MCH", "HMT102": "MCHC",
    "RCT6": "Urea nitrogen", "RCT11": "Serum glucose", "RCT20": "Cholesterol",
    "RCT392": "Creatinine",
}


def build_lab(df_lab: pd.DataFrame, df_date_src: pd.DataFrame = None) -> pd.DataFrame:
    """
    all_{code}_values + all_{code}_date (aligned). aibl_labdata has no date:
    EXAMDATE is joined from aibl_mmse on (RID, VISCODE), like the diagnosis.
    """
    _print("  §5  Lab data...")
    if df_lab.empty:
        return pd.DataFrame()

    lab = normalise_id(df_lab.copy())
    if df_date_src is not None and not df_date_src.empty and "EXAMDATE" not in lab.columns:
        ds = normalise_id(df_date_src.copy())
        if {"PTID", "VISCODE", "EXAMDATE"} <= set(ds.columns):
            ds = ds[["PTID", "VISCODE", "EXAMDATE"]].drop_duplicates(["PTID", "VISCODE"])
            lab["VISCODE"] = lab["VISCODE"].astype(str)
            ds["VISCODE"]  = ds["VISCODE"].astype(str)
            lab = lab.merge(ds, on=["PTID", "VISCODE"], how="left")
    if "EXAMDATE" not in lab.columns:
        lab["EXAMDATE"] = pd.NaT
    safe_date(lab, "EXAMDATE")
    lab["_order"] = pd.to_numeric(lab.get("VISCODE", 0), errors="coerce").fillna(0)
    lab = lab.sort_values(["PTID", "EXAMDATE", "_order"])

    records = []
    for ptid, grp in lab.groupby("PTID"):
        row = {"ID": ptid}
        for col in LAB_COLS:
            if col not in grp.columns:
                continue
            g = grp.assign(_v=pd.to_numeric(grp[col], errors="coerce"))
            g = g[g["_v"].notna() & (g["_v"] >= 0)]
            if g.empty:
                continue
            dts = ["NA" if pd.isna(d) else str(d.date()) for d in g["EXAMDATE"]]
            row[f"all_{col}_values"]   = "|".join(map(str, g["_v"]))
            row[f"all_{col}_date"]     = "|".join(dts)
            row[f"nb_of_{col}_values"] = len(g)
        records.append(row)

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(records)
    return result.set_index("ID") if "ID" in result.columns else pd.DataFrame()


# =============================================================================
# §5b MEDICAL HISTORY — neurologic history proxy for seizures (added 2026-11-08)
# =============================================================================
# aibl_medhist: RID, SITEID, VISCODE, MHPSYCH, MH2NEURL, MH4CARD, MH6HEPAT,
#               MH8MUSCL, MH9ENDO, MH10GAST, MH12RENA, MH16SMOK, MH17MALI
# AIBL has no seizure-specific field. MH2NEURL (1=yes, 0=no; "neurologic other
# than AD") is the closest proxy → Neuro_history. MH9ENDO (endocrine-metabolic,
# includes diabetes) → Endo_history.
# =============================================================================

MEDHIST_COLS = {"MH2NEURL": "Neuro_history", "MH9ENDO": "Endo_history"}


def build_medhist(df_mh: pd.DataFrame, df_date_src: pd.DataFrame = None) -> pd.DataFrame:
    _print("  §5b Medical history (MH2NEURL / MH9ENDO)...")
    if df_mh is None or df_mh.empty:
        _print("    [INFO] aibl_medhist not found → skipped")
        return pd.DataFrame()
    mh = normalise_id(df_mh.copy())
    if df_date_src is not None and not df_date_src.empty and "EXAMDATE" not in mh.columns:
        ds = normalise_id(df_date_src.copy())
        if {"PTID", "VISCODE", "EXAMDATE"} <= set(ds.columns):
            ds = ds[["PTID", "VISCODE", "EXAMDATE"]].drop_duplicates(["PTID", "VISCODE"])
            mh["VISCODE"] = mh["VISCODE"].astype(str)
            ds["VISCODE"] = ds["VISCODE"].astype(str)
            mh = mh.merge(ds, on=["PTID", "VISCODE"], how="left")
    if "EXAMDATE" not in mh.columns:
        mh["EXAMDATE"] = pd.NaT
    safe_date(mh, "EXAMDATE")
    mh["_order"] = pd.to_numeric(mh.get("VISCODE", 0), errors="coerce").fillna(0)
    mh = mh.sort_values(["PTID", "EXAMDATE", "_order"])
    records = []
    for ptid, grp in mh.groupby("PTID"):
        row = {"ID": ptid}
        for col, nm in MEDHIST_COLS.items():
            if col not in grp.columns:
                continue
            g = grp.assign(_v=pd.to_numeric(grp[col], errors="coerce"))
            g = g[g["_v"].isin([0, 1])]
            if g.empty:
                continue
            dts  = ["NA" if pd.isna(d) else str(d.date()) for d in g["EXAMDATE"]]
            vals = g["_v"].astype(int).tolist()
            row[f"all_{nm}_date"]     = "|".join(dts)
            row[f"all_{nm}_score"]    = "|".join(map(str, vals))
            row[f"1st_{nm}_score"]    = vals[0]
            row[f"1st_{nm}_date"]     = dts[0]
            row[f"Lst_{nm}_score"]    = vals[-1]
            row[f"Lst_{nm}_date"]     = dts[-1]
            row[f"nb_of_{nm}_visits"] = len(vals)
            row[f"{nm}_ever"]         = int(max(vals))
        records.append(row)
    return pd.DataFrame(records).set_index("ID") if records else pd.DataFrame()


# =============================================================================
# §6  DEMOGRAPHICS + APOE
# =============================================================================
# aibl_ptdemog: RID, SITEID, VISCODE, PTGENDER (1=M 2=F), PTDOB (YYYY-MM-DD)
#   One row per visit (take first per subject).
#
# aibl_apoeres: RID, SITEID, VISCODE, APTESTDT, APGEN1, APGEN2, APRECEIVE, APAMBTEMP
#   APGEN1 / APGEN2 are allele integers (2, 3, or 4).
#   APOE_genotype = f"{APGEN1}/{APGEN2}"
#   APOE_e4_count = (APGEN1==4) + (APGEN2==4)
# =============================================================================

def build_demographics(df_demog: pd.DataFrame,
                       df_apoe:  pd.DataFrame) -> pd.DataFrame:
    _print("  §6  Demographics + APOE...")
    records: dict = {}

    # ── ptdemog ──────────────────────────────────────────────────────────────
    if not df_demog.empty:
        demog = normalise_id(df_demog.copy())
        # First visit per subject (lowest VISCODE)
        demog["_vc"] = pd.to_numeric(demog.get("VISCODE", 0), errors="coerce").fillna(0)
        first = demog.sort_values("_vc").drop_duplicates("PTID")

        for _, row in first.iterrows():
            ptid = row["PTID"]
            # DOB: PTDOB column is "YYYY-MM-DD" string
            dob = pd.NaT
            ptdob = row.get("PTDOB", None)
            if pd.notna(ptdob) and str(ptdob).strip() not in ("", "nan"):
                dob = pd.to_datetime(str(ptdob).strip(), errors="coerce")

            records[ptid] = {
                "ID"         : ptid,
                "gender"     : row.get("PTGENDER", np.nan),   # 1=Male 2=Female
                "birth_date" : dob,
            }

    # ── apoeres ──────────────────────────────────────────────────────────────
    if not df_apoe.empty:
        apoe = normalise_id(df_apoe.copy())
        # Keep first non-null result per subject
        apoe = apoe.dropna(subset=["APGEN1", "APGEN2"], how="all")
        apoe["_vc"] = pd.to_numeric(apoe.get("VISCODE", 0), errors="coerce").fillna(0)
        apoe = apoe.sort_values("_vc").drop_duplicates("PTID")

        for _, row in apoe.iterrows():
            ptid = row["PTID"]
            rec  = records.setdefault(ptid, {"ID": ptid})
            g1   = int(row.get("APGEN1", 0) or 0)
            g2   = int(row.get("APGEN2", 0) or 0)
            rec["APOE_genotype"] = f"{g1}/{g2}" if g1 and g2 else np.nan
            rec["APOE_e4_count"] = int(g1 == 4) + int(g2 == 4)

    df = pd.DataFrame(list(records.values()))
    return df.set_index("ID") if not df.empty else df


# =============================================================================
# §7  AGE AT FIRST DIAGNOSIS (post-merge)
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
    _print("AIBL LONGITUDINAL INTEGRATION PIPELINE  v2")
    _print("=" * 62)

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load ──────────────────────────────────────────────────────────────────
    _print("\nLoading raw tables...")
    dfs = {k: load(k) for k in PATHS}

    # ── Build domain tables ───────────────────────────────────────────────────
    _print("\nBuilding domain tables...")
    diag  = build_diagnosis(dfs["dx"], dfs["dx_date_src"])
    mri   = canonise_meta_types(build_mri(dfs["mri_15t"], dfs["mri_3t"]))
    pet   = build_pet(dfs)

    # IDA search export = primary imaging source; the *meta tables only cover a
    # subject x type that the export does not have (whole blocks, never mixed).
    ref_ids = next((normalise_id(dfs[k])["PTID"] for k in ("dx", "demog", "mmse")
                    if not dfs[k].empty and "PTID" in normalise_id(dfs[k]).columns), None)
    ida_img = build_imaging_ida(ref_ids)
    if not ida_img.empty:
        mri_cols = [c for c in ida_img.columns if "MRI_" in c or c == "IDA_sex"]
        pet_cols = [c for c in ida_img.columns
                    if "ABETA_" in c or "TAU_" in c or "PET_" in c]
        n_meta = len(set(mri.index) | set(pet.index))
        mri = ida_combine_blocks(ida_img[mri_cols], mri)
        pet = ida_combine_blocks(ida_img[pet_cols], pet)
        if ref_ids is not None:
            only_ida = len(set(ida_img.index) - set(ref_ids))
            _print(f"    [INFO] imaging subjects: {len(ida_img):,} in the IDA export "
                  f"({only_ida:,} of them absent from the clinical tables), "
                  f"{n_meta:,} in the *meta tables")
    else:
        _print("    [INFO] no IDA search export — MRI / PET from the *meta tables only")
    clin  = build_all_clinical(dfs)
    lab   = build_lab(dfs["labdata"], dfs["dx_date_src"])
    medh  = build_medhist(dfs.get("medhist", pd.DataFrame()), dfs["dx_date_src"])
    demog = build_demographics(dfs["demog"], dfs["apoe"])

    # ── Merge ─────────────────────────────────────────────────────────────────
    _print("\nMerging all domains on subject ID...")
    parts = [p for p in [diag, mri, pet, clin, lab, medh, demog]
             if isinstance(p, pd.DataFrame) and not p.empty]

    if not parts:
        _print("ERROR: no data loaded — check paths.")
        return pd.DataFrame()

    master = reduce(lambda a, b: a.join(b, how="outer"), parts)
    master.index.name = "ID"
    master = master.reset_index()
    gc.collect()

    # ── Post-processing ───────────────────────────────────────────────────────
    master = add_age_at_1st_diag(master)

    master["Plasma_available"] = "not_available_in_AIBL_public (routine lab only)"

    LEAD = [
        "ID",
        "gender", "birth_date", "age_at_1st_diag",
        "APOE_genotype", "APOE_e4_count",
        "1st_diag", "Lst_diag", "1st_diag_date", "Lst_diag_date",
        "all_diag", "all_diag_date", "nb_of_visits",
        "1st_time_of_change", "1st_time_of_change_days_to",
        "all_time_of_change", "all_time_of_change_days_to",
        "Lst_time_of_change", "nb_of_change", "Diag_total",
        "MRI_type_disponible", "PET_tracer_disponible",
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
        _print("\nMRI types:")
        cnt: Counter = Counter()
        for v in master["MRI_type_disponible"].dropna():
            cnt.update(v.split("|"))
        for t, n in cnt.most_common():
            _print(f"  {t:<35s} {n:>5,} subjects")

    if "PET_tracer_disponible" in master.columns:
        _print("\nPET tracers:")
        cnt = Counter()
        for v in master["PET_tracer_disponible"].dropna():
            cnt.update(v.split("|"))
        for t, n in cnt.most_common():
            _print(f"  {t:<35s} {n:>5,} subjects")

    return master


# =============================================================================
