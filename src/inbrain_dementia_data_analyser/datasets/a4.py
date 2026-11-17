"""A4 LONGITUDINAL INTEGRATION PIPELINE  v1
Applies the same domain-builder logic as ADNI v2, adapted for A4's
multi-CSV structure with BID-based subject IDs and EXAMDAY time axis.

A4-SPECIFIC ARCHITECTURE
--------------------------
A4 ships multiple CSV files (one per instrument form), all keyed on:
  - BID      : Blinded Participant ID  → normalised to PTID
  - VISCODE  : Visit code (scv1=screening1, scv2=screening2,
                           scv3=baseline, scv4=follow-up1, scv5=follow-up2)
  - EXAMDAY  : Days since consent (integer) — primary time axis
               No real calendar dates for most tables.

VISIT CODE MAPPING:
  scv1 → Screening 1          (~day -60 to -30)
  scv2 → Screening 2          (~day -30 to 0)
  scv3 → Baseline             (day 0)
  scv4 → Follow-up 1          (~18 months)
  scv5 → Follow-up 2          (~36 months)

VISCODE ORDER:  scv1 < scv2 < scv3 < scv4 < scv5
Ordering uses EXAMDAY when available; otherwise VISCODE integer position.

DIAGNOSIS:
  A4 is a PREVENTION trial — all subjects are cognitively unimpaired (CU)
  at enrollment. Diagnosis is STATIC: "CU" for all subjects.
  The trial enrolled amyloid-positive CU adults (confirmed via PET).
  No DXSUM-style longitudinal diagnosis table exists.

SOURCE FILES (all in Tables_A4/):
  REGISTRY         : BID, VISCODE, EXAMDAY, VISTYPE    (visit registry)
  PTDEMOG          : demographics + EXAMDAY
  SUBJINFO         : BID, AGEYR, APOEGN               (age at consent, APOE)
  CDR              : CDGLOBAL, CDSOB
  MMSE             : MMSCORE (calculated from item responses)
  COGDIGIT         : DIGITTOTAL (Digit Symbol Substitution, 0-93)
  COGLOGIC         : LIMMTOTAL (LM immediate), LDELTOTAL (LM delayed)
  COGFCSR16        : FCTOTF (free recall), FCTOTC (cued recall), FCTOTAL
  PSYCHWELL        : GDTOTAL (GDS-15), STAITOTAL (STAI-6)
  CFI              : CFIPTTOTAL (Cognitive Function Index, patient)
  CFISP            : CFSPTTOTAL (Cognitive Function Index, study partner)
  ADLPQ            : Activities of Daily Living (17 items, scored 0-3)
  VITALS           : VSWEIGHT, VSHEIGHT, VSBPSYS, VSBPDIA, VSPULSE
  HABITS           : SMOKE, ALCOHOL, CAFFEINE, AEROBIC, WALKING, SLEEP
  IES              : IESCORE (Impact of Event Scale total)
  FTPSCALE         : Future Time Perspective Scale (items)
  CONCERNS         : Alzheimer's disease concerns (items)
  SPPACC           : PACC composite (Preclinical Alzheimer Cognitive Composite)

Imaging-derived tables (were in Tables_A4/Matias Files/ — that folder no
longer exists; the files are read from Tables_A4/, any export date):
  A4_VMRI_PRV2_09Jun2026.csv   : BID, VISCODE + 50 volumetric MRI cols
  TAUSUVR_09Jun2026.csv        : ID + ~120 regional tau SUVR cols
  TAUSUVR_PETSURFER_09Jun2026.csv : ID + ~120 regional tau SUVR (PETsurfer)

PETSUVR table:
  BID, VISCODE, ligand, scan_analyzed, brain_region, suvr_cer, centiloid
  Ligands: Florbetapir (AV45) for amyloid; Flortaucipir (AV1451) for tau.

DWI/fMRI params: BID, VISCODE, series_id, SeriesDescription,
  Manufacturer, imageid (LONI ImageID) — used for MRI section.

IDA search export (idaSearch_*.csv in Tables_A4/, single file with
  'Visit', 'Study Date', 'Image ID' …; found by ida_load): every MRI / PET
  series with its LONI Image ID. 'Study Date' is year-only (01/01/YYYY) in
  A4, so the time axis stays EXAMDAY: 'SCV4 (Screening MRI)' -> registry
  VISCODE 4 -> EXAMDAY. 'Imaging Protocol' is empty in the A4 export (no
  scanner there; the vendor still comes from the DWI / fMRI params tables).

Output:
  <root>/Tables_A4/Integrated_A4/A4_longitudinal.csv

Author  : VRS_inbrain project
Version : 2  —  2026-11-16  (single-file IDA search: Image IDs, registry EXAMDAY)

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("A4", root_path)
"""
from __future__ import annotations

import gc
import numpy as np
import os
import pandas as pd
import re
from functools import reduce

from ..datasets.base import DatasetContext
from ..fluid import (
    plasma_add_summary,
    plasma_canon_analyte,
    plasma_clean_numeric,
    plasma_detect_long,
    plasma_from_spec,
    plasma_join,
    plasma_read_csv,
    plasma_sanitize,
)
from ..imaging.ida import (
    ida_build_acquisitions,
    ida_classify,
    ida_combine_blocks,
    ida_load,
    ida_report_unclassified,
)
from ..imaging.scanner import add_scanner_channels, norm_vendor
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

DATASET = "A4"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
A4 = OUT_DIR = OUT_LONG = ""


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, A4, OUT_DIR, OUT_LONG, A4_PLASMA_DIR, A4_PLASMA_FILE_RX
    _CTX = ctx
    A4 = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()
    A4_PLASMA_DIR = str(ctx.folder)
    A4_PLASMA_FILE_RX = re.compile(ctx.settings.get(
        "plasma_file_regex", "plasma|ptau|p-tau|pt217|biomarker|roche|c2n|nfl|gfap"), re.I)


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write the per-subject table of A4; returns it."""
    configure(ctx)
    return main()



def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns of one measure (see series.build_scored_columns)."""
    return build_scored_columns(scores_by_subj, test_name)

A4_PLASMA_DIR = ""
A4_PLASMA_FILE_RX = re.compile("plasma", re.I)




# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================













# VISCODE → approximate month offset for ordering when EXAMDAY is missing
VISCODE_ORDER = {
    "scv1": 0, "scv2": 1, "scv3": 2, "scv4": 20, "scv5": 38,
}

# ADLPQ item columns (14 items, scored 0-3 each, max 42)


# VMRI key ROI columns
VMRI_KEY_ROIS = [
    "LeftHippocampus", "RightHippocampus",
    "LeftEntorhinal",  "RightEntorhinal",
    "LeftAmygdala",    "RightAmygdala",
    "IntraCranialVolume", "HOC",
    "LeftCorticalGrayMatter", "RightCorticalGrayMatter",
    "LeftWMHypo", "RightWMHypo",
]

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================










def load(key: str) -> pd.DataFrame:
    """Multi-strategy CSV loader with 5 fallback strategies."""
    path = PATHS.get(key, "")
    if not path or not os.path.exists(path):
        _print(f"  [SKIP] {key}: not found")
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
                _print(f"  [WARN] {key}: loaded with fallback (strategy {i+1})")
            _print(f"  [OK]   {key}: {len(df):>6,} rows  {df.shape[1]} cols")
            return df
        except Exception as e:
            last_err = e
    _print(f"  [ERROR] {key}: all strategies failed — {last_err}")
    return pd.DataFrame()


def normalise_id(df: pd.DataFrame) -> pd.DataFrame:
    """Rename BID (or ID) → PTID for consistency."""
    if "PTID" in df.columns:
        return df
    col_lower = {c.lower(): c for c in df.columns}
    for alt in ("bid", "ptid", "id", "subject_id"):
        if alt in col_lower:
            df = df.rename(columns={col_lower[alt]: "PTID"})
            return df
    return df


def attach_examday(df: pd.DataFrame,
                   registry: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure df has an _examday column for sorting.
    Priority: EXAMDAY column → join from registry → VISCODE_ORDER approximation.
    """
    df = normalise_id(df)

    if "EXAMDAY" in df.columns:
        df["_examday"] = pd.to_numeric(df["EXAMDAY"], errors="coerce")
        return df

    if not registry.empty and "VISCODE" in df.columns:
        reg = normalise_id(registry.copy())
        reg = reg.rename(columns={"BID": "PTID"}) if "BID" in reg.columns else reg
        reg = reg.drop_duplicates(["PTID", "VISCODE"])
        reg["EXAMDAY"] = pd.to_numeric(reg.get("EXAMDAY", pd.Series()),
                                        errors="coerce")
        # Coerce both VISCODE columns to lowercase str to prevent
        # int64/object type mismatch error during merge
        reg["VISCODE"] = reg["VISCODE"].astype(str).str.strip().str.lower()
        df["VISCODE"]  = df["VISCODE"].astype(str).str.strip().str.lower()
        df = df.merge(reg[["PTID", "VISCODE", "EXAMDAY"]],
                      on=["PTID", "VISCODE"], how="left")

    if "EXAMDAY" not in df.columns:
        df["EXAMDAY"] = np.nan

    if "VISCODE" in df.columns:
        vc_order = df["VISCODE"].map(VISCODE_ORDER).fillna(99)
        df["_examday"] = df["EXAMDAY"].fillna(vc_order * 30)
    else:
        df["_examday"] = df["EXAMDAY"].fillna(0)

    df["_examday"] = pd.to_numeric(df["_examday"], errors="coerce").fillna(0)
    return df





def _nps_from_df(df: pd.DataFrame, registry: pd.DataFrame,
                 test_name: str, score_col: str) -> pd.DataFrame:
    """Extract one NPS test longitudinally from a per-visit A4 form table."""
    if df.empty or score_col not in df.columns:
        _print(f"    [WARN] '{score_col}' not found → {test_name} skipped")
        return pd.DataFrame()
    df = normalise_id(df.copy())
    df = attach_examday(df, registry)
    df = df.sort_values(["PTID", "_examday"])
    scores_by_subj = {}
    for subj, grp in df.groupby("PTID"):
        scores_by_subj[subj] = list(zip(grp["_examday"], grp[score_col]))
    return build_nps_columns(scores_by_subj, test_name)


# =============================================================================
# §1  DIAGNOSIS — STATIC CU
# =============================================================================

def build_diagnosis(registry: pd.DataFrame) -> pd.DataFrame:
    _print("  §1  Diagnosis (static CU — A4 prevention trial)...")
    if registry.empty:
        return pd.DataFrame()

    reg = normalise_id(registry.copy())
    reg["_examday"] = pd.to_numeric(reg.get("EXAMDAY", 0),
                                     errors="coerce").fillna(0)
    reg = reg.sort_values(["PTID", "_examday"])

    records = []
    for ptid, grp in reg.groupby("PTID"):
        baseline = grp[grp.get("VISCODE", pd.Series()).eq("scv3")] \
                   if "VISCODE" in grp.columns else pd.DataFrame()
        baseline_day = (baseline["_examday"].iloc[0]
                        if not baseline.empty
                        else grp["_examday"].min())
        records.append({
            "ID"                         : ptid,
            "1st_diag"                   : "CU",
            "Lst_diag"                   : "CU",
            "1st_diag_date"              : baseline_day,
            "Lst_diag_date"              : grp["_examday"].max(),
            "all_diag"                   : "CU",
            "all_diag_date"              : list_str(grp["_examday"].tolist()),
            "nb_of_visits"               : len(grp),
            "1st_time_of_change"         : np.nan,
            "1st_time_of_change_days_to" : np.nan,
            "all_time_of_change"         : np.nan,
            "all_time_of_change_days_to" : np.nan,
            "Lst_time_of_change"         : np.nan,
            "nb_of_change"               : 0,
            "Diag_total"                 : "CU",
        })

    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §2  MRI — DWI Params + fMRI Params + VMRI volumetric
# =============================================================================

MRI_TYPE_MAP_A4 = {
    "mprage": "MRI_T1w", "t1_mprage": "MRI_T1w", "t1mprage": "MRI_T1w",
    "mp-rage": "MRI_T1w", "fl3d": "MRI_T1w",
    "t2_flair": "MRI_T2w_FLAIR", "flair": "MRI_T2w_FLAIR",
    "t2_spc": "MRI_T2w", "t2_tse": "MRI_T2w",
    "bold": "MRI_fMRI", "resting": "MRI_fMRI", "rest": "MRI_fMRI",
    "fieldmap": None, "field_map": None, "gre_field": None,
    "scout": None, "localizer": None, "mocoseries": None,
}


_A4_SCANNER_ROWS = None   # PTID / _examday / _T / Manufacturer rows (set by build_mri)


def build_mri(dfs: dict, registry: pd.DataFrame) -> pd.DataFrame:
    _print("  §2  MRI acquisitions...")
    global _A4_SCANNER_ROWS
    _A4_SCANNER_ROWS = None
    records = {}

    # ── DWI params ────────────────────────────────────────────────────────────
    dwi = dfs.get("dwi_params", pd.DataFrame())
    if not dwi.empty:
        dwi = normalise_id(dwi)
        dwi = attach_examday(dwi, registry)
        for _, row in dwi.sort_values(["PTID", "_examday"]).iterrows():
            ptid = row["PTID"]
            r = records.setdefault(ptid, {"ID": ptid,
                                           "MRI_type_disponible": set(),
                                           "_dwi_ids": [], "_dwi_days": []})
            img_id = row.get("imageid", None) or \
                     f"{ptid}_{row.get('VISCODE','?')}_{row.get('series_id','?')}"
            r["_dwi_ids"].append(str(img_id))
            r["_dwi_days"].append(row["_examday"])
            r["MRI_type_disponible"].add("MRI_DTI")

    # ── fMRI params ───────────────────────────────────────────────────────────
    fmri = dfs.get("fmri_params", pd.DataFrame())
    if not fmri.empty:
        fmri = normalise_id(fmri)
        fmri = attach_examday(fmri, registry)
        for _, row in fmri.sort_values(["PTID", "_examday"]).iterrows():
            ptid = row["PTID"]
            r = records.setdefault(ptid, {"ID": ptid,
                                           "MRI_type_disponible": set()})
            img_id = row.get("imageid", None) or \
                     f"{ptid}_{row.get('VISCODE','?')}_{row.get('series_id','?')}"
            r.setdefault("_fmri_ids", []).append(str(img_id))
            r.setdefault("_fmri_days", []).append(row["_examday"])
            r["MRI_type_disponible"].add("MRI_fMRI")

    # ── VMRI volumetric ───────────────────────────────────────────────────────
    vmri = dfs.get("vmri", pd.DataFrame())
    if not vmri.empty:
        vmri = normalise_id(vmri)
        vmri = attach_examday(vmri, registry)
        for ptid, grp in vmri.sort_values(["PTID", "_examday"]).groupby("PTID"):
            r = records.setdefault(ptid, {"ID": ptid,
                                           "MRI_type_disponible": set()})
            r["MRI_type_disponible"].add("MRI_T1w")
            days = grp["_examday"].tolist()
            ids  = [f"{ptid}_vmri_v{row.get('VISCODE','?')}"
                    for _, row in grp.iterrows()]
            r["all_MRI_T1w_visit"]  = list_str(days)
            r["nb_total_MRI_T1w"]   = len(days)
            r["MRI_T1w_IDs"]        = list_str(ids)
            r["1st_MRI_T1w_visit"]  = days[0]  if days else np.nan
            r["Lst_MRI_T1w_visit"]  = days[-1] if days else np.nan
            # Key ROI metrics
            for roi in VMRI_KEY_ROIS:
                if roi in vmri.columns:
                    cn = f"VMRI_{roi}"
                    vals = grp[roi].dropna().tolist()
                    r[f"all_{cn}"] = list_str(vals)
                    r[f"1st_{cn}"] = vals[0]  if vals else np.nan
                    r[f"Lst_{cn}"] = vals[-1] if vals else np.nan

    # Collapse accumulation lists
    for ptid, r in records.items():
        for mtype, id_key, day_key in [
            ("MRI_DTI",  "_dwi_ids",  "_dwi_days"),
            ("MRI_fMRI", "_fmri_ids", "_fmri_days"),
        ]:
            il = r.pop(id_key,  [])
            dl = r.pop(day_key, [])
            if il:
                r[f"all_{mtype}_visit"] = list_str(dl)
                r[f"nb_total_{mtype}"]  = len(il)
                r[f"{mtype}_IDs"]       = list_str(il)
                r[f"1st_{mtype}_ID"]    = il[0]
                r[f"Lst_{mtype}_ID"]    = il[-1]

        r["MRI_type_disponible"] = (
            "|".join(sorted(r["MRI_type_disponible"])) or np.nan
        )

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(list(records.values()))
    if "ID" not in result.columns:
        return pd.DataFrame()
    result = result.set_index("ID")

    # Scanner per acquisition: DWI / fMRI params carry 'Manufacturer';
    # the T1w (VMRI) of the same PTID+VISCODE was acquired on that scanner.
    sc = []
    for df_, T in ((dwi, "MRI_DTI"), (fmri, "MRI_fMRI")):
        if not df_.empty and "Manufacturer" in df_.columns:
            sc.append(df_.assign(_T=T)[["PTID", "VISCODE", "_examday", "_T", "Manufacturer"]])
    if sc:
        sc = pd.concat(sc, ignore_index=True)
        if not vmri.empty and "VISCODE" in vmri.columns:
            vm = sc.dropna(subset=["Manufacturer"]).assign(
                VISCODE=lambda x: x["VISCODE"].astype(str).str.lower()) \
                .groupby(["PTID", "VISCODE"])["Manufacturer"].first()
            t1 = vmri[["PTID", "VISCODE", "_examday"]].copy()
            t1["VISCODE"] = t1["VISCODE"].astype(str).str.lower()
            t1["Manufacturer"] = [vm.get((p, v), np.nan) for p, v in zip(t1["PTID"], t1["VISCODE"])]
            t1["_T"] = "MRI_T1w"
            sc = pd.concat([sc, t1], ignore_index=True)
        result = add_scanner_channels(result, sc, "PTID", "_T", "_examday",
                                      "Manufacturer", None, list_suffix="visit")
        _A4_SCANNER_ROWS = sc            # reused by main() after the IDA merge
    return result


# =============================================================================
# §3  PET — amyloid (AV45) + tau (AV1451) + regional tau SUVR tables
# =============================================================================

PET_LIGAND_MAP = {
    "florbetapir"  : "ABETA_AV45",
    "av45"         : "ABETA_AV45",
    "18f-av45"     : "ABETA_AV45",
    "av-45"        : "ABETA_AV45",
    "flortaucipir" : "TAU_Flortaucipir",
    "av1451"       : "TAU_Flortaucipir",
    "18f-av1451"   : "TAU_Flortaucipir",
    "av-1451"      : "TAU_Flortaucipir",
}

TAU_KEY_ROIS_FSL = [
    "Mean_Left_Hippocampus", "Mean_Right_Hippocampus",
    "Mean_Left_Amygdala",    "Mean_Right_Amygdala",
    "Mean_Left_Entorhinal",  "Mean_Right_Entorhinal",
    "Mean_ctx_lh_inferiortemporal", "Mean_ctx_rh_inferiortemporal",
]

TAU_KEY_ROIS_PS = [
    "bi_Hippocampus", "bi_Amygdala", "bi_entorhinal",
    "bi_inferiortemporal", "bi_middletemporal", "bi_parahippocampal",
]


def _norm_ligand(raw) -> str | None:
    if pd.isna(raw):
        return None
    t = str(raw).strip().lower()
    for key, canon in PET_LIGAND_MAP.items():
        if key in t:
            return canon
    return None


def build_pet(dfs: dict, registry: pd.DataFrame) -> pd.DataFrame:
    _print("  §3  PET acquisitions (amyloid + tau)...")
    records = {}

    # ── PETSUVR per-visit table ───────────────────────────────────────────────
    pet = dfs.get("petsuvr", pd.DataFrame())
    if not pet.empty:
        pet = normalise_id(pet)
        pet = attach_examday(pet, registry)
        pet["_ligand"] = pet.get("ligand", pd.Series()).apply(_norm_ligand)
        pet = pet[pet["_ligand"].notna()]

        suvr_col = next((c for c in pet.columns if "suvr" in c.lower()), None)
        centiloid_col = next((c for c in pet.columns
                              if "centiloid" in c.lower()), None)

        for (ptid, vc, ligand), grp in pet.groupby(
                ["PTID", "_examday", "_ligand"]):
            r = records.setdefault(ptid, {
                "ID": ptid, "PET_tracer_disponible": set()
            })
            r["PET_tracer_disponible"].add(ligand)
            mean_suvr = grp[suvr_col].mean() if suvr_col else np.nan
            mean_cl   = grp[centiloid_col].mean() if centiloid_col else np.nan
            r.setdefault(f"_pet_{ligand}_days", []).append(vc)
            r.setdefault(f"_pet_{ligand}_suvr", []).append(mean_suvr)
            r.setdefault(f"_pet_{ligand}_cl",   []).append(mean_cl)

    # ── TAUSUVR regional tables ──────────────────────────────────────────────────
    for tau_key, key_rois, label in [
        ("tausuvr",    TAU_KEY_ROIS_FSL, "TAU_FSL"),
        ("tausuvr_ps", TAU_KEY_ROIS_PS,  "TAU_PS"),
    ]:
        tau = dfs.get(tau_key, pd.DataFrame())
        if tau.empty:
            continue
        tau = normalise_id(tau)
        for _, row in tau.iterrows():
            ptid = row["PTID"]
            r = records.setdefault(ptid, {
                "ID": ptid, "PET_tracer_disponible": set()
            })
            r["PET_tracer_disponible"].add("TAU_Flortaucipir")
            for roi in key_rois:
                if roi in tau.columns:
                    clean = roi.replace("Mean_", "").replace("ctx_", "")\
                               .replace("lh_", "L_").replace("rh_", "R_")
                    r[f"{label}_{clean}"] = row.get(roi, np.nan)

    # Collapse accumulation lists
    for ptid, r in records.items():
        r["PET_tracer_disponible"] = (
            "|".join(sorted(r["PET_tracer_disponible"])) or np.nan
        )
        for ligand in set(PET_LIGAND_MAP.values()):
            days_k = f"_pet_{ligand}_days"
            suvr_k = f"_pet_{ligand}_suvr"
            cl_k   = f"_pet_{ligand}_cl"
            if days_k in r:
                days = r.pop(days_k)
                suvr = r.pop(suvr_k, [])
                cl   = r.pop(cl_k, [])
                r[f"all_{ligand}_visit"]     = list_str(days)
                r[f"nb_total_{ligand}"]      = len(days)
                r[f"1st_{ligand}_visit"]     = days[0]  if days else np.nan
                r[f"Lst_{ligand}_visit"]     = days[-1] if days else np.nan
                r[f"all_{ligand}_suvr"]      = list_str(suvr)
                r[f"all_{ligand}_centiloid"] = list_str(cl)

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(list(records.values()))
    return result.set_index("ID") if "ID" in result.columns else pd.DataFrame()


# =============================================================================
# §4  CLINICAL / NPS / FUNCTIONAL / VITALS / HABITS
# =============================================================================
#
# NPS tests: (table_key, score_col, canonical_name)
#
# CDR          → CDGLOBAL (CDR), CDSOB (CDR_SB), MEMORY, ORIENT
# MMSE         → MMSCORE
# COGDIGIT     → DIGITTOTAL (Digit Symbol, 0-93)
# COGLOGIC     → LIMMTOTAL (LM_I), LDELTOTAL (LM_II)
# COGFCSR16    → FCTOTF (free recall), FCTOTC (cued), FCTOTAL
# PSYCHWELL    → GDTOTAL (GDS-15), STAITOTAL (STAI-6)
# CFI          → CFIPTTOTAL (patient CFI)
# CFISP        → CFSPTTOTAL (study partner CFI)
# IES          → IESCORE
# VITALS       → VSWEIGHT, VSHEIGHT, VSBPSYS, VSBPDIA, VSPULSE
# HABITS       → SLEEP, AEROBIC, WALKING, ALCOHOL, CAFFEINE, SMOKE
# ADLPQ        → sum of 14 items → ADLPQ_total
# SPPACC       → PACCRN (normalised PACC composite)
# =============================================================================




_LEGACY_DOUBLE_ATTACH = {"ADLPQ", "PACC"}


def build_all_clinical(dfs: dict, registry: pd.DataFrame) -> pd.DataFrame:
    """Clinical / NPS / functional / vitals / habits series.

    Variables come from ``variables.yaml`` (``sources: A4:``); every series is
    timed with the registry ``EXAMDAY`` of its visit (days since consent).
    """
    _print("  §4  Clinical / NPS / Functional / Vitals / Habits...")
    parts = []
    for name, src in variable_sources(DATASET, config=_CTX.config):
        df = dfs.get(src.table, pd.DataFrame())
        if df.empty:
            _print(f"    [WARN] table '{src.table}' empty or missing -> {name} skipped")
            continue
        frame, col = materialise(df, name, src)
        if col is None:
            continue
        if name in _LEGACY_DOUBLE_ATTACH:
            # The notebook attached EXAMDAY to these tables before calling
            # _nps_from_df, which attaches it again. On the second pass the
            # EXAMDAY column already exists, so visits whose VISCODE has no
            # registry day keep a missing time instead of the 99 x 30 = 2970-day
            # placeholder used for every other measure. Kept to reproduce the
            # notebook output; see the report (A4 placeholder days).
            frame = attach_examday(normalise_id(frame.copy()), registry)
        part = _nps_from_df(frame, registry, name, col)
        if not part.empty:
            parts.append(part)

    # Initial health history (added 2026-11-08). A4 has no seizure field and
    # the public INITHEALTH file carries no free-text description, so only the
    # body-system category is usable:
    #   IHSYMPTOM 2 = Neurologic (other than AD) → Neuro_history_*
    #   IHSYMPTOM 9 = Endocrine-Metabolic        → Endo_history_*
    # Subject-level: *_ever (1/0), *_n (records), *_current (IHCUR=1 on any).
    ih = dfs.get("inithealth", pd.DataFrame())
    if not ih.empty and "IHSYMPTOM" in ih.columns:
        ih = normalise_id(ih.copy())
        ih["_cat"] = pd.to_numeric(ih["IHSYMPTOM"], errors="coerce")
        ih["_cur"] = pd.to_numeric(ih.get("IHCUR", np.nan), errors="coerce")
        recs = []
        for subj, g in ih.groupby("PTID"):
            r = {"ID": subj}
            for code, nm in ((2, "Neuro_history"), (9, "Endo_history")):
                h = g[g["_cat"] == code]
                r[f"{nm}_ever"]    = int(len(h) > 0)
                r[f"{nm}_n"]       = int(len(h))
                r[f"{nm}_current"] = int((h["_cur"] == 1).any())
            recs.append(r)
        if recs:
            parts.append(pd.DataFrame(recs).set_index("ID"))

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()
    return reduce(lambda a, b: a.join(b, how="outer"), parts)


# =============================================================================
# §5  DEMOGRAPHICS + APOE
# =============================================================================
# PTDEMOG: BID, VISCODE, EXAMDAY, PTGENDER (1=M 2=F), PTAGE, PTETHNIC,
#          PTRACE, PTEDUCAT, PTMARRY, PTNOTRT, PTHOME, PTLANG
# SUBJINFO: BID, AGEYR (age at consent), APOEGN (e.g. "3/4")

def build_demographics(dfs: dict, registry: pd.DataFrame) -> pd.DataFrame:
    _print("  §5  Demographics + APOE...")
    records: dict = {}

    dem = dfs.get("ptdemog", pd.DataFrame())
    if not dem.empty:
        dem = normalise_id(dem.copy())
        dem = attach_examday(dem, registry)
        first = dem.sort_values("_examday").drop_duplicates("PTID", keep="first")
        for _, row in first.iterrows():
            ptid = row["PTID"]
            records[ptid] = {
                "ID"             : ptid,
                "gender"         : row.get("PTGENDER",  np.nan),
                "race"           : row.get("PTRACE",    np.nan),
                "hispanic"       : row.get("PTETHNIC",  np.nan),
                "education"      : row.get("PTEDUCAT",  np.nan),
                "marital_status" : row.get("PTMARRY",   np.nan),
                "retired"        : row.get("PTNOTRT",   np.nan),
                "residence_type" : row.get("PTHOME",    np.nan),
                "language"       : row.get("PTLANG",    np.nan),
                "age_at_screen"  : row.get("PTAGE",     np.nan),
            }

    subj = dfs.get("subjinfo", pd.DataFrame())
    if not subj.empty:
        subj = normalise_id(subj.copy())
        for _, row in subj.drop_duplicates("PTID").iterrows():
            ptid = row["PTID"]
            rec  = records.setdefault(ptid, {"ID": ptid})
            rec["age_at_consent"] = row.get("AGEYR", np.nan)
            geno = str(row.get("APOEGN", "") or "").strip()
            rec["APOE_genotype"] = geno if geno and geno.lower() != "nan" else np.nan
            rec["APOE_e4_count"] = geno.count("4") if geno else np.nan

    df = pd.DataFrame(list(records.values()))
    return df.set_index("ID") if not df.empty else df


# =============================================================================
# §6  AGE AT FIRST DIAGNOSIS
# =============================================================================

def add_age_at_1st_diag(master: pd.DataFrame) -> pd.DataFrame:
    if "age_at_screen" in master.columns:
        master["age_at_1st_diag"] = master["age_at_screen"]
    elif "age_at_consent" in master.columns:
        master["age_at_1st_diag"] = master["age_at_consent"]
    return master



# =============================================================================
# ── IDA IMAGING SECTION (LONI IDA search export) ─────────────────────────────
# =============================================================================
# File: the idaSearch_*.csv of Tables_A4/ (discovered by ida_load, shared-
# helpers cell). Columns: 'Subject ID', 'Sex', 'Research Group', 'Visit',
#   'Study Date', 'Age', 'Modality', 'Description', 'Type',
#   'Imaging Protocol', 'Image ID'
#
# DESCRIPTION COLUMN FORMAT: "TYPE_PREFIX <- Free text description"
#   The Type column contains only "Pre-processed" (not useful).
#   The Description column encodes both type and scan name:
#     "Florbetapir <- PET Scan"            → ABETA_AV45
#     "Flortaucipir <- PET Scan"           → TAU_Flortaucipir
#     "DWI <- Axial DWI"                   → MRI_DTI
#     "fMRI_rest <- Axial fcMRI (EYES OPEN)"→ MRI_fMRI
#     "T1; GradWarp; DeFaced <- MPRAGE"    → MRI_T1w
#     "T1; DeFaced <- MPRAGE"              → MRI_T1w
#     "T2_SE; DeFaced <- Axial T2 TSE"     → MRI_T2w
#     "FLAIR; DeFaced <- Axial T2-FLAIR"   → MRI_T2w_FLAIR
#     "T2_star <- AXIAL_T2_STAR"           → MRI_T2star
#
#   Strategy: extract the part BEFORE " <- " and apply prefix match.
#   If no " <- " separator, fall back to matching the whole string.
#
# VISIT COLUMN → registry VISCODE (numeric in A4_REGISTRY):
#   'SCV4 (Screening MRI)' → VISCODE 4
#   'SCV2 (Screening PET)' → VISCODE 2
#   EXAMDAY is then read from the registry row of that subject + VISCODE.
#   (Before 2026-11-16 the visit was matched as the text 'scv4', which never
#    equals the numeric VISCODE, so every IDA acquisition fell back to the
#    constants 600 / 30 / 2970 days.)
#
# STUDY DATE: year-only in A4 (01/01/YYYY) → not used; time axis = EXAMDAY.
# ARCHIVE DATE (export of 2026-11-16): date the series was UPLOADED to LONI
#   (2020-2021), one batch per sequence type (all T2* on one day, all DWI on
#   another …) → spreads 1-241 days inside ONE visit; not a scan date, not used.
# AGE: 0.0 in ~700 rows (mostly LEARN) = missing → NaN (ida_load).
# IMAGE ID  : real LONI Image ID of every series → {TYPE}_IDs.
# =============================================================================

# Type prefix (from Description before " <- ") → canonical
IDA_TYPE_MAP = [
    ("florbetapir",  "ABETA_AV45"),
    ("flortaucipir", "TAU_Flortaucipir"),
    ("dwi",          "MRI_DTI"),
    ("fmri_rest",    "MRI_fMRI"),
    ("t1",           "MRI_T1w"),   # T1; GradWarp; DeFaced, T1; DeFaced
    ("t2_se",        "MRI_T2w"),   # T2_SE; DeFaced
    ("flair",        "MRI_T2w_FLAIR"),
    ("t2_star",      "MRI_T2star"),
]


A4_DESC_SEP = re.compile(r"\s*(?:<-|\u21d0|\u2190)\s*")


def classify_a4_type(desc_str: str) -> str | None:
    """
    Map an A4 Description string to a canonical type.

    Description format: "TYPE_PREFIX <- Free text scan name"
    Strategy: extract TYPE_PREFIX (part before ' <- '), lowercase,
    then apply prefix match against IDA_TYPE_MAP.

    Examples:
      "Florbetapir <- PET Scan"             → ABETA_AV45
      "T1; GradWarp; DeFaced <- MPRAGE"     → MRI_T1w
      "T2_SE; DeFaced <- Axial T2 TSE..."   → MRI_T2w
      "FLAIR; DeFaced <- Axial T2-FLAIR"    → MRI_T2w_FLAIR
      "T2_star <- AXIAL_T2_STAR"            → MRI_T2star
    """
    if not isinstance(desc_str, str) or not desc_str.strip():
        return None

    # Prefix before the separator. The CSV export writes " <- "; the IDA web
    # table shows the same thing as "⇐" — accept both (and "←").
    prefix = A4_DESC_SEP.split(desc_str, maxsplit=1)[0].strip().lower()

    for key, canon in IDA_TYPE_MAP:
        if prefix.startswith(key):
            return canon
    return None


def normalise_a4_visit(visit_str: str) -> str | None:
    """
    Registry VISCODE of an IDA Visit label (the registry stores it as a number).
    'SCV4 (Screening MRI)' → '4'
    'SCV2 (Screening PET)' → '2'
    """
    if not isinstance(visit_str, str) or not visit_str.strip():
        return None
    m = re.search(r"(\d+)", visit_str.strip().split()[0])
    return m.group(1) if m else None


def build_imaging_ida_a4(registry: pd.DataFrame) -> pd.DataFrame:
    """
    MRI + PET acquisition table from the IDA search export.

    Uses:
      - 'Description' prefix for the canonical type
      - 'Visit' → registry VISCODE → EXAMDAY (time axis)
      - 'Subject ID' as PTID (BID), 'Image ID' as LONI Image ID

    A SESSION = subject x type x EXAMDAY. Columns per canonical type
    (every list has nb_total_{TYPE} tokens):
      all_{TYPE}_visit        EXAMDAY of each session
      {TYPE}_IDs              LONI Image ID of each session
      all_{TYPE}_image_ids    every Image ID of the session ('/'-joined)
      nb_total_{TYPE}, nb_images_{TYPE}
      1st_/Lst_{TYPE}_visit, 1st_/Lst_{TYPE}_ID, 1st_{TYPE}_age
      MRI_type_disponible, PET_tracer_disponible
    """
    ida = ida_load(A4, label="A4")
    if ida.empty:
        _print(f"    [SKIP] no IDA search export in → {A4}")
        return pd.DataFrame()

    ida["canon"] = ida_classify(ida, classify_a4_type)
    _print(f"    [INFO] Classified {int(ida['canon'].notna().sum()):,}/{len(ida):,} rows → "
          f"{ida['canon'].value_counts().to_dict()}")
    ida_report_unclassified(ida, "canon")

    # Visit → registry VISCODE → EXAMDAY of that subject's visit
    ida["VISCODE"] = ida["visit"].map(normalise_a4_visit)
    ida["EXAMDAY"] = np.nan
    if not registry.empty and {"PTID", "VISCODE", "EXAMDAY"}.issubset(registry.columns):
        reg = registry[["PTID", "VISCODE", "EXAMDAY"]].copy()
        reg["PTID"]    = reg["PTID"].astype(str).str.strip()
        reg["VISCODE"] = (reg["VISCODE"].astype(str).str.strip().str.lower()
                          .str.replace(r"\.0$", "", regex=True))
        reg["EXAMDAY"] = pd.to_numeric(reg["EXAMDAY"], errors="coerce")
        reg = reg.dropna(subset=["EXAMDAY"]).drop_duplicates(["PTID", "VISCODE"])
        lut = reg.set_index(["PTID", "VISCODE"])["EXAMDAY"]
        key = pd.MultiIndex.from_arrays([ida["subject"], ida["VISCODE"].fillna("")])
        ida["EXAMDAY"] = lut.reindex(key).values
        miss = ida["EXAMDAY"].isna() & ida["canon"].notna()
        if miss.any():
            # no registry row for that subject+visit → typical day of the visit
            med = reg.groupby("VISCODE")["EXAMDAY"].median()
            ida.loc[miss, "EXAMDAY"] = ida.loc[miss, "VISCODE"].map(med)
            _print(f"    [INFO] {int(miss.sum()):,} IDA rows "
                  f"({ida.loc[miss, 'subject'].nunique():,} subjects) without a registry "
                  f"row for their visit → EXAMDAY = median day of that VISCODE")
    else:
        _print("    [WARN] registry without PTID/VISCODE/EXAMDAY → IDA sessions untimed")

    ida["visit_tok"] = ida["EXAMDAY"]
    ida["visit_ord"] = ida["EXAMDAY"]
    return ida_build_acquisitions(ida, "canon", use_dates=False, first_age=True)


# =============================================================================
# MAIN
# =============================================================================


# =============================================================================
# =============================================================================
# §5  PLASMA BIOMARKERS  (auto-discovered — not part of the PRV2 DATADIC)
# =============================================================================
# The A4/LEARN plasma results (e.g. Lilly p-tau217, Roche Elecsys plasma
# panel, C2N) are distributed as separate "biomarker_*.csv" files outside the
# PRV2 data dictionary. Any CSV under Tables_A4/ whose name matches
# A4_PLASMA_FILE_RX is read:
#   • long format  (TESTCD/LBTESTCD + ORRES/LBORRES…)  → one series per test
#   • wide format  numeric columns whose name matches A4_PLASMA_COL_RX
# Series name: Plasma_{test_or_column}_{file_tag}; time axis = EXAMDAY.
# Requires the "Plasma / lab helpers" cell.
# =============================================================================



A4_PLASMA_COL_RX  = re.compile(r"p-?tau|pt217|ab_?4[02]|abeta|a[bβ]42|a[bβ]40|nfl|nf_?l|gfap|t-?tau|aps", re.I)
A4_PLASMA_SKIP_COL_RX = re.compile(r"^(bid|ptid|id|viscode|examday|visit|sample|spid|subjid|.*date.*)$", re.I)


def _a4_plasma_files() -> list:
    known = {os.path.abspath(p) for p in PATHS.values()}
    out = []
    for root, _, files in os.walk(A4_PLASMA_DIR):
        if "Integrated_A4" in root:
            continue
        for fn in files:
            p = os.path.abspath(os.path.join(root, fn))
            if fn.lower().endswith(".csv") and A4_PLASMA_FILE_RX.search(fn) and p not in known:
                out.append(p)
    return sorted(out)


def build_plasma_a4(registry: pd.DataFrame) -> pd.DataFrame:
    _print("  §5  Plasma biomarkers (auto-discovery)...")
    files = _a4_plasma_files()
    if not files:
        _print(f"    [INFO] no plasma files found in {A4_PLASMA_DIR} "
              f"(pattern: {A4_PLASMA_FILE_RX.pattern})")
        return pd.DataFrame()

    parts = []
    for path in files:
        raw = plasma_read_csv(path)
        if raw.empty:
            continue
        raw = normalise_id(raw)
        if "PTID" not in raw.columns:
            _print(f"    [SKIP] {os.path.basename(path)}: no BID/ID column")
            continue
        raw = attach_examday(raw, registry)
        tag = re.sub(r"(?i)\.csv$|^biomarker_?|plasma_?|_?results?", "", os.path.basename(path))
        tag = plasma_sanitize(tag, 20) if tag.strip("_ ") else "A4"
        visit_col = "VISCODE" if "VISCODE" in raw.columns else None

        tcol, vcol = plasma_detect_long(raw)
        if tcol and vcol:
            # avoid 'Plasma_pTau217_pTau217' when the file is named after the analyte
            tag_ = "A4" if plasma_canon_analyte(tag) in set(raw[tcol].map(plasma_canon_analyte)) else tag
            spec = dict(long=True, test_col=tcol, value_col=vcol,
                        name_tpl="Plasma_{test}_" + tag_)
        else:
            cols = [c for c in raw.columns
                    if A4_PLASMA_COL_RX.search(c) and not A4_PLASMA_SKIP_COL_RX.match(c)
                    and plasma_clean_numeric(raw[c]).notna().any()]
            spec = dict(markers=[(c, f"Plasma_{plasma_canon_analyte(c)}_{tag}") for c in cols])
        got = plasma_from_spec(raw, spec, "PTID", "_examday", visit_col)
        parts.extend(got)
        _print(f"    [OK]   {os.path.basename(path)}: {len(raw):,} rows → {len(got)} marker(s)")

    return plasma_add_summary(plasma_join(parts))


def main() -> pd.DataFrame:
    _print("=" * 62)
    _print("A4 LONGITUDINAL INTEGRATION PIPELINE  v2")
    _print("=" * 62)
    _print("\nNOTE: A4 is a prevention trial — all subjects are CU (amyloid+).")
    _print("      Time axis: EXAMDAY (days since consent), not calendar dates.")
    _print("      APOE from SUBJINFO.APOEGN; age from PTDEMOG.PTAGE.")

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load all tables ───────────────────────────────────────────────────────
    _print("\nLoading raw tables...")
    dfs = {k: load(k) for k in PATHS}

    registry = normalise_id(dfs.get("registry", pd.DataFrame()).copy()) \
               if not dfs.get("registry", pd.DataFrame()).empty \
               else pd.DataFrame()

    n_subjects = registry["PTID"].nunique() if not registry.empty else "unknown"
    _print(f"\n  Subjects (from registry): {n_subjects}")

    # ── Build domains ─────────────────────────────────────────────────────────
    _print("\nBuilding domain tables...")

    # ── IDA imaging (LONI IDA search export) ─────────────────────────────────
    _print("  Loading IDA search export...")
    ida_imaging = build_imaging_ida_a4(registry)

    diag  = build_diagnosis(registry)
    mri   = build_mri(dfs, registry)
    pet   = build_pet(dfs, registry)
    clin  = build_all_clinical(dfs, registry)
    demog = build_demographics(dfs, registry)
    plasma = build_plasma_a4(registry)

    # Merge IDA results at the subject x type level (lists of two sources are
    # never mixed inside one type):
    #   MRI: IDA first (real LONI Image IDs, one session per EXAMDAY); the
    #        DWI / fMRI / VMRI tables cover subjects or types absent from IDA.
    #   PET: PETSUVR first (its SUVR / centiloid lists are aligned with its own
    #        visits); IDA fills missing tracers and lends the Image IDs when
    #        both have the same number of sessions.
    if not ida_imaging.empty:
        mri_cols_ida = [c for c in ida_imaging.columns if "MRI_" in c or c == "IDA_sex"]
        pet_cols_ida = [c for c in ida_imaging.columns
                        if "ABETA_" in c or "TAU_" in c or c == "PET_tracer_disponible"]
        if mri_cols_ida:
            mri = ida_combine_blocks(ida_imaging[mri_cols_ida], mri)
        if pet_cols_ida:
            pet = ida_combine_blocks(pet, ida_imaging[pet_cols_ida], borrow_ids=True)
    else:
        _print("  [INFO] IDA file not found — using DWI/fMRI/VMRI tables only")

    # Scanner vendor (DWI / fMRI params) for every MRI type of the session.
    # A4 has one MRI visit per subject, so a subject with a single vendor gets
    # it on all of its MRI sessions (robust to EXAMDAY differences between the
    # tables and the registry); subjects with >1 vendor are matched by EXAMDAY.
    sc = globals().get("_A4_SCANNER_ROWS")
    if isinstance(sc, pd.DataFrame) and not sc.empty and not mri.empty:
        sc = sc.dropna(subset=["Manufacturer"]).copy()
        sc["_vend"] = sc["Manufacturer"].map(norm_vendor)
        n_vend = sc.dropna(subset=["_vend"]).groupby("PTID")["_vend"].nunique()
        single = (sc.dropna(subset=["_vend"]).drop_duplicates("PTID")
                    .set_index("PTID")["_vend"])
        single = single[single.index.isin(n_vend[n_vend == 1].index)]
        mri_types = [c[len("nb_total_"):] for c in mri.columns if c.startswith("nb_total_MRI_")]
        rows = []
        for T in mri_types:
            lc = f"all_{T}_visit"
            if lc not in mri.columns:
                continue
            cells = mri[lc].dropna()
            cells = cells[cells.index.isin(single.index)]
            for sid, cell in cells.items():
                rows += [(sid, T, tok, single[sid]) for tok in str(cell).split("|")]
        sc_all = pd.DataFrame(rows, columns=["PTID", "_T", "_examday", "Manufacturer"])
        multi = sc[~sc["PTID"].isin(single.index)]
        if not multi.empty and mri_types:
            multi = multi.drop_duplicates(["PTID", "_examday"])[["PTID", "_examday", "Manufacturer"]]
            sc_all = pd.concat([sc_all] + [multi.assign(_T=T) for T in mri_types],
                               ignore_index=True)
        if not sc_all.empty:
            mri = add_scanner_channels(mri, sc_all, "PTID", "_T", "_examday",
                                       "Manufacturer", None, list_suffix="visit")

    # ── Merge ─────────────────────────────────────────────────────────────────
    _print("\nMerging all domains on subject ID...")
    parts = [p for p in [diag, mri, pet, clin, plasma, demog]
             if isinstance(p, pd.DataFrame) and not p.empty]

    if not parts:
        _print("ERROR: no domain tables produced — check paths.")
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
    master["CSF_available"] = "not_available_in_A4"

    master["Plasma_available"] = (
        master["Plasma_panels_disponible"].notna().map({True: "yes", False: "no"})
        if "Plasma_panels_disponible" in master.columns else "no")

    LEAD = [
        "ID",
        "gender", "race", "hispanic", "education",
        "marital_status", "retired", "residence_type", "language",
        "age_at_screen", "age_at_consent", "age_at_1st_diag",
        "APOE_genotype", "APOE_e4_count",
        "1st_diag", "Lst_diag", "1st_diag_date", "Lst_diag_date",
        "all_diag", "all_diag_date", "nb_of_visits",
        "nb_of_change", "Diag_total",
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
    _print(f"\nDiagnosis: all subjects CU (n={len(master):,})")
    n_ida = len(ida_imaging) if not ida_imaging.empty else 0
    _print(f"Subjects with IDA imaging data: {n_ida:,}")

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
                _print(f"  {name:<25s} {n:>5,} subjects")

    return master
