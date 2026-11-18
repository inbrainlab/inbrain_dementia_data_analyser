"""ABVIB LONGITUDINAL INTEGRATION PIPELINE  v2
MRI acquisitions come from the LONI IDA search export. Since the 2026-11-16
download this is ONE file (any idaSearch_*.csv in Tables_ABVIB/ that has
'Study Date' + 'Image ID'; found automatically by ida_load) with columns:
  'Subject ID', 'Sex', 'Visit', 'Study Date', 'Archive Date', 'Age',
  'Modality', 'Description', 'Type', 'Imaging Protocol', 'Image ID'
The five old per-visit files (idaSearch_6_09_2026_visit_{1..5}.csv, only
'Subject ID', 'Sex', 'Age', 'Description', 'Type') now live in
old_ida_search/ and are read only to fill rows missing from the new export.

Description → canonical MRI type mapping (from full inventory analysis):
  MRI_T1w       : MPRAGE*, t1_mpr_sag*, T1 SPGR, SAG T1 SPGR
  MRI_T2w       : t2_tse3d*, HighResCoronalHippo*, T2 TSE*, cor HIPPO,
                  axial Hippo, tse2d*, t2w_3d*
  MRI_T2w_FLAIR : t2_tirm3d*, Axial T2-FLAIR, 3D T2 FLAIR, 3D_T2_FLAIR*,
                  flair_1mm, Axial FLAIR
  MRI_DTI       : DTI*, ep2d_diff*  (excludes _ADC, _FA, _TRACEW)
  MRI_fMRI      : rs-fmri*, rs-fMRI, rs_fmri*
  MRI_ASL       : 3DASL, ASL CBF, ep2d_tra_pasl, MN-ep2d_pcasl*,
                  fbirn_pasl*, Perfusion_Weighted, relCBF
  MRI_T2star    : AXIAL_T2_STAR, T2Star GRE

SKIP (infra): MoCoSeries, scout*, loc, Circle Scout, ASSET cal*,
              Sag Cal*, Localizer, localizer
SKIP (derived): Fractional Ansio., T2-weighted trace, Average DC*,
                *_ADC

Session rule: one session per subject x type x Study Date (real dates, same
calendar as the UDS visit dates); the IDA visit number is kept as a token.
Subject IDs: the IDA 'Subject ID' is cast to the dtype of the clinical PTID
(integer) before the join — as text it never matched and every IDA subject
appeared twice in the output (one clinical row + one imaging row).

Author  : VRS_inbrain project
Version : 3  —  2026-11-16  (single-file IDA search: dates, Image IDs, scanner)

Run through the package::

    from inbrain_dementia_data_analyser import preprocess_dataset
    preprocess_dataset("ABVIB", root_path)
"""
from __future__ import annotations

import gc
import numpy as np
import os
import pandas as pd
import re
from functools import reduce

from ..clinical import AED_TEXT_RX, text_flag_series
from ..datasets.base import DatasetContext
from ..diagnosis import add_dx_columns, trajectory_record, uds_d1_stage_etiology
from ..imaging.ida import (
    ida_build_acquisitions,
    ida_classify,
    ida_combine_blocks,
    ida_load,
    ida_match_index,
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

DATASET = "ABVIB"

# Module state set by configure(); the functions below read these names.
_CTX: DatasetContext | None = None
PATHS: dict = {}
ABVIB = OUT_DIR = OUT_LONG = ""


def configure(ctx: DatasetContext) -> None:
    """Point the module at one root folder / configuration."""
    global _CTX, PATHS, ABVIB, OUT_DIR, OUT_LONG
    _CTX = ctx
    ABVIB = str(ctx.folder)
    OUT_DIR, OUT_LONG = str(ctx.output_dir), str(ctx.output_path)
    PATHS = ctx.paths()


def run(ctx: DatasetContext) -> pd.DataFrame:
    """Build and write the per-subject table of ABVIB; returns it."""
    configure(ctx)
    return main()



def build_nps_columns(scores_by_subj: dict, test_name: str) -> pd.DataFrame:
    """Scored-family columns of one measure (see series.build_scored_columns)."""
    return build_scored_columns(scores_by_subj, test_name)




# =============================================================================
# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# =============================================================================











# IDA search: the idaSearch_*.csv of Tables_ABVIB/ is discovered by ida_load()
# (shared-helpers cell) — no download date hard-coded here.
def _abvib_visit_num(v):
    """IDA Visit label 'Visit 3' -> 3 (visit token)."""
    m = re.search(r"(\d+)", str(v))
    return int(m.group(1)) if m else np.nan


def _abvib_legacy_visit(fname: str) -> str:
    """Old per-visit file name -> Visit label of the new export."""
    m = re.search(r"_visit_(\d+)", fname, re.I)
    return f"Visit {int(m.group(1))}" if m else ""


MRI_TYPES_ABVIB = [
    "MRI_T1w", "MRI_T2w", "MRI_T2w_FLAIR",
    "MRI_DTI", "MRI_fMRI", "MRI_ASL", "MRI_T2star",
]

# =============================================================================
# ── CANONICAL TYPE RULES ──────────────────────────────────────────────────────
# =============================================================================

SKIP_INFRA = [
    "scout", "loc", "moco", "asset cal", "asset cal.", "sag cal",
    "circle scout", "localizer", "t2 tse axial hp scout",
    "t2_tse_axial_hp_scout", "t2 scout",
]
SKIP_DERIVED_EXACT = [
    "fractional ansio.", "t2-weighted trace",
]
SKIP_DERIVED_CONTAINS = ["average dc", "average dcs", "average dcls"]
SKIP_DERIVED_SUFFIX   = ["_adc", "_fa", "_tracew"]


def classify_abvib_description(desc: str) -> str | None:
    """
    Map a Description string to a canonical MRI type, or None to skip.
    Rules derived from full ABVIB IDA description inventory.
    """
    if not isinstance(desc, str) or not desc.strip():
        return None

    dl = desc.lower().strip()
    du = desc.upper().strip()

    # ── SKIP: infrastructure ─────────────────────────────────────────────────
    if any(p in dl for p in SKIP_INFRA):
        return None
    # ── SKIP: derived maps ───────────────────────────────────────────────────
    if dl in SKIP_DERIVED_EXACT:
        return None
    if any(x in dl for x in SKIP_DERIVED_CONTAINS):
        return None
    if any(dl.endswith(s) for s in SKIP_DERIVED_SUFFIX):
        return None

    # ── DTI — primary only ────────────────────────────────────────────────────
    if any(x in dl for x in ["dti", "ep2d_diff"]):
        return "MRI_DTI"

    # ── T1w ──────────────────────────────────────────────────────────────────
    if any(x in dl for x in ["mprage", "t1 spgr", "t1_mpr", "t1_spgr",
                               "sag t1 spgr", "spgr"]):
        return "MRI_T1w"

    # ── T2w FLAIR (before T2w) ────────────────────────────────────────────────
    if any(x in dl for x in ["flair", "tirm", "t2_tirm"]):
        return "MRI_T2w_FLAIR"

    # ── T2w structural ────────────────────────────────────────────────────────
    if any(x in dl for x in ["t2_tse", "t2 tse", "t2w_3d", "tse2d",
                               "cor hippo", "highresco", "axial hippo",
                               "t2w_3d"]):
        return "MRI_T2w"

    # ── T2* / SWI ─────────────────────────────────────────────────────────────
    if any(x in dl for x in ["t2star", "t2_star", "axial_t2_star",
                               "t2star gre"]):
        return "MRI_T2star"

    # ── fMRI / resting-state ──────────────────────────────────────────────────
    if any(x in dl for x in ["rs-fmri", "rs_fmri", "rsfmri", "rs-fmri"]):
        return "MRI_fMRI"

    # ── ASL / perfusion ───────────────────────────────────────────────────────
    if any(x in dl for x in ["asl", "pasl", "pcasl", "perfusion",
                               "relcbf", "cbf", "3dasl", "fbirn_pasl",
                               "ep2d_tra_pasl", "mn-ep2d"]):
        return "MRI_ASL"

    return None   # unknown — skip silently


def is_infra_abvib(desc) -> bool:
    """True for scouts / calibration / derived DTI maps (skipped on purpose)."""
    if not isinstance(desc, str) or not desc.strip():
        return True
    dl = desc.lower().strip()
    return (any(p in dl for p in SKIP_INFRA) or dl in SKIP_DERIVED_EXACT
            or any(x in dl for x in SKIP_DERIVED_CONTAINS)
            or any(dl.endswith(x) for x in SKIP_DERIVED_SUFFIX))


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


def load_ida() -> pd.DataFrame:
    """IDA search of ABVIB, classified (canon) with visit-number tokens.
    Requires the shared-helpers cell (ida_load, ida_classify …)."""
    ida = ida_load(ABVIB, legacy_visit=_abvib_legacy_visit, label="ABVIB")
    if ida.empty:
        return ida
    ida["canon"]     = ida_classify(ida, classify_abvib_description)
    ida["visit_tok"] = ida["visit"].map(_abvib_visit_num)
    ida["visit_ord"] = ida["visit_tok"]
    return ida


def normalise_id(df: pd.DataFrame) -> pd.DataFrame:
    if "PTID" in df.columns:
        return df
    col_lower = {c.lower(): c for c in df.columns}
    for alt in ("ptid", "subject_id", "subid", "rid"):
        if alt in col_lower:
            df = df.rename(columns={col_lower[alt]: "PTID"})
            return df
    return df


def make_visit_date(df: pd.DataFrame) -> pd.DataFrame:
    if all(c in df.columns for c in ["VISITYR", "VISITMO"]):
        yr = pd.to_numeric(df["VISITYR"],  errors="coerce").fillna(0).astype(int)
        mo = pd.to_numeric(df["VISITMO"],  errors="coerce").fillna(1).astype(int).clip(1, 12)
        dy = pd.to_numeric(df.get("VISITDAY", pd.Series(1, index=df.index)),
                           errors="coerce").fillna(1).astype(int).clip(1, 28)
        df["_date"] = pd.to_datetime(dict(year=yr, month=mo, day=dy),
                                      errors="coerce")
    else:
        df["_date"] = pd.NaT
    return df





def _nps_from_df(df: pd.DataFrame, test_name: str, score_col: str,
                 date_col: str = "_date") -> pd.DataFrame:
    if score_col not in df.columns:
        _print(f"    [WARN] '{score_col}' not found → {test_name} skipped")
        return pd.DataFrame()
    actual_date = date_col if date_col in df.columns else "_date"
    sub = df[["PTID", actual_date, score_col]].copy()
    safe_date(sub, actual_date)
    sub = sub.sort_values(["PTID", actual_date])
    scores_by_subj = {}
    for subj, grp in sub.groupby("PTID"):
        scores_by_subj[subj] = list(zip(grp[actual_date], grp[score_col]))
    return build_nps_columns(scores_by_subj, test_name)


# =============================================================================
# §1  DIAGNOSIS TRAJECTORY (unchanged from v1)
# =============================================================================

MCI_COLS = [
    "MCIAMEM", "MCIAPLUS", "MCIAPLAN", "MCIAPATT", "MCIAPEX", "MCIAPVIS",
    "MCINON1", "MCIN1LAN", "MCIN1ATT", "MCIN1EX",  "MCIN1VIS",
    "MCINON2", "MCIN2LAN", "MCIN2ATT", "MCIN2EX",  "MCIN2VIS", "IMPNOMCI",
]
AD_COLS  = ["PROBAD", "POSSAD"]
LBD_COLS = ["DLB"]
VAD_COLS = ["VASC", "VASCPS"]
FTD_COLS = ["FTD", "PSP", "CORT", "PPAPH", "PNAPH", "SEMDEMAN", "SEMDEMAG",
            "PPAOTHR"]
HD_COLS  = ["HUNT"]
PD_COLS  = ["PARK"]


def _classify_abvib(row: pd.Series) -> str:
    def active(cols):
        return any(
            pd.notna(row.get(c)) and float(row.get(c, 0) or 0) == 1
            for c in cols
        )
    if active(["NORMCOG"]):
        return "CN"
    parts = []
    if active(MCI_COLS):  parts.append("MCI")
    if active(AD_COLS):   parts.append("AD")
    if active(LBD_COLS):  parts.append("LBD")
    if active(VAD_COLS):  parts.append("VaD")
    if active(FTD_COLS):  parts.append("FTD")
    if active(HD_COLS):   parts.append("HD")
    if active(PD_COLS):   parts.append("PD")
    if not parts:
        if pd.notna(row.get("DEMENTED")) and float(row.get("DEMENTED", 0) or 0) == 1:
            return "Dementia-NOS"
        return "Unknown"
    return "+".join(sorted(parts))


def build_diagnosis(df_raw: pd.DataFrame) -> pd.DataFrame:
    _print("  §1  Diagnosis trajectory...")
    dx = normalise_id(df_raw.copy())
    dx = make_visit_date(dx)
    dx["DX_LABEL"] = dx.apply(_classify_abvib, axis=1)   # pre-0.4.0 label -> *_diag_legacy
    dx = add_dx_columns(dx, uds_d1_stage_etiology)  # UDS D1 v1.2-v3 (upper/lower case)
    dx = dx.sort_values(["PTID", "_date",
                          "VISITNUM"] if "VISITNUM" in dx.columns
                        else ["PTID", "_date"])
    records = []
    for ptid, grp in dx.groupby("PTID"):
        grp    = grp.sort_values("_date")
        dates  = grp["_date"].tolist()
        records.append(trajectory_record(ptid, grp, dates))
    result = pd.DataFrame(records)
    return result.set_index("ID") if not result.empty else result


# =============================================================================
# §2  MRI — IDA search export (primary) + UDS_E1 flag fallback
# =============================================================================
# PRIMARY: IDA search export — one row per image series.
#   A SESSION = subject x type x Study Date (ida_build_acquisitions).
#
# FALLBACK: UDS_E1_IMAGING_LABS — MRI1FLM/MRI1DIG presence flags.
#   Used for subjects not in the IDA export (older visits, different format).
#
# Columns produced per canonical type (every list has nb_total_{TYPE} tokens):
#   all_{TYPE}_date        study date of each session (YYYY-MM-DD)
#   all_{TYPE}_visit       IDA visit number of each session
#   {TYPE}_IDs             LONI Image ID of each session
#   all_{TYPE}_image_ids   every Image ID of the session ('/'-joined)
#   all_{TYPE}_manufacturer / _model / _field_strength, {TYPE}_manufacturer
#   nb_total_{TYPE}, nb_images_{TYPE}, 1st_/Lst_{TYPE}_visit, 1st_/Lst_{TYPE}_ID
#   MRI_type_disponible — pipe-sorted types present
# =============================================================================

def build_mri_ida(ida: pd.DataFrame, ref_ids=None) -> pd.DataFrame:
    """MRI sessions from the IDA search export, index cast to the PTID dtype
    of the clinical tables (ref_ids) so both sides join on the same key."""
    if ida is None or ida.empty:
        return pd.DataFrame()
    ida_report_unclassified(ida, "canon", skip_fn=is_infra_abvib)
    result = ida_build_acquisitions(ida, "canon")
    if not result.empty and ref_ids is not None:
        result = ida_match_index(result, ref_ids)
    return result


def build_mri_uds_flags(df_raw: pd.DataFrame,
                         ida_subjects: set) -> pd.DataFrame:
    """
    Fallback: build MRI T1w presence from UDS_E1 flags for subjects
    NOT already covered by IDA files.
    """
    if df_raw.empty:
        return pd.DataFrame()

    mri = normalise_id(df_raw.copy())
    mri = make_visit_date(mri)

    mri_cols = [c for c in ["MRI1FLM", "MRI1DIG", "MRI2FLM", "MRI2DIG"]
                if c in mri.columns]
    if not mri_cols:
        return pd.DataFrame()

    for c in mri_cols:
        mri[c] = pd.to_numeric(mri[c], errors="coerce")

    mri["_mri_done"] = mri[mri_cols].eq(1).any(axis=1)
    mri = mri[mri["_mri_done"] & ~mri["PTID"].isin(ida_subjects)].copy()

    if mri.empty:
        return pd.DataFrame()

    def synth_id(r):
        ds = str(r["_date"].date()) if pd.notna(r["_date"]) \
             and hasattr(r["_date"], "date") else "nodate"
        return f"{r['PTID']}_{ds}"

    mri["_img_id"] = mri.apply(synth_id, axis=1)
    mri = mri.sort_values(["PTID", "_date"])

    records = {}
    for ptid, grp in mri.groupby("PTID"):
        row = records.setdefault(ptid, {
            "ID": ptid, "MRI_type_disponible": "MRI_T1w",
            "_dates": [], "_ids": [],
        })
        for _, r in grp.iterrows():
            ds = d_fmt(r["_date"])
            if ds and pd.notna(ds):
                row["_dates"].append(ds)
            row["_ids"].append(r["_img_id"])

    for ptid, row in records.items():
        dl = row.pop("_dates", [])
        il = row.pop("_ids",   [])
        uniq = list(dict.fromkeys(dl))
        row["all_MRI_T1w_date"]  = list_str(uniq)
        row["nb_total_MRI_T1w"]  = len(uniq)
        row["MRI_T1w_IDs"]       = list_str(il)
        row["1st_MRI_T1w_visit"] = 1   # visit number unknown from flags
        row["Lst_MRI_T1w_visit"] = len(uniq)

    if not records:
        return pd.DataFrame()
    result = pd.DataFrame(list(records.values()))
    return result.set_index("ID") if "ID" in result.columns else pd.DataFrame()


# =============================================================================
# §3–§5  CLINICAL / NPS / DEMOGRAPHICS (unchanged from v1)
# =============================================================================








# SUBCAT (enrolment categorisation; added 2026-11-08)
#   FASTING : fasting status at the lab draw (pairs with LABORATORY GLUCOSE)
#   DIABETS : diabetes at enrolment



def _prepare_table(key: str, df: pd.DataFrame) -> pd.DataFrame:
    """Normalise the subject ID and build the visit date ``_date`` of one table.

    UDS forms: VISITYR / VISITMO / VISITDAY. Tables listed in the dataset
    setting ``table_date_columns`` (COGSCOR, LABORATORY, SUBCAT) carry their
    own date column.
    """
    df = normalise_id(df.copy())
    date_cols = _CTX.settings.get("table_date_columns", {})
    if key in date_cols:
        col = date_cols[key]
        safe_date(df, col)
        df["_date"] = df[col] if col in df.columns else pd.NaT
    else:
        df = make_visit_date(df)
    return df


def build_all_clinical(dfs: dict) -> pd.DataFrame:
    """Clinical / NPS / lab series (``variables.yaml``, ``sources: ABVIB:``)."""
    _print("  §3  Clinical / NPS / Functional / Vitals...")
    parts = []
    prepared: dict = {}
    for name, src in variable_sources(DATASET, config=_CTX.config):
        if src.table not in prepared:
            raw = dfs.get(src.table, pd.DataFrame())
            prepared[src.table] = raw if raw.empty else _prepare_table(src.table, raw)
        df = prepared[src.table]
        if df.empty:
            continue
        frame, col = materialise(df, name, src)
        if col is None:
            continue
        parts.append(_nps_from_df(frame, name, col, src.date_column or "_date"))
    # ── UDS A4 medications: antiepileptic text flag (added 2026-11-08) ──────
    # DRUGID may hold NACC drug codes (dXXXXX) rather than names — then only
    # OTHERMEDICINE (free text) can match. Check AED_text_n after the run.
    meds = dfs.get("meds", pd.DataFrame())
    if not meds.empty:
        meds = normalise_id(meds.copy())
        meds = make_visit_date(meds)
        parts.append(text_flag_series(meds, "PTID", ["DRUGID", "OTHERMEDICINE"],
                                      AED_TEXT_RX, "AED_text", "_date"))

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()
    return reduce(lambda a, b: a.join(b, how="outer"), parts)


def build_demographics(df_demog, df_subcat, df_apoe) -> pd.DataFrame:
    _print("  §4  Demographics + APOE...")
    records: dict = {}
    if not df_demog.empty:
        dem = normalise_id(df_demog.copy())
        dem = make_visit_date(dem)
        first = dem.sort_values("_date").drop_duplicates("PTID", keep="first")
        for _, row in first.iterrows():
            ptid = row["PTID"]
            dob = pd.NaT
            try:
                yr = int(row.get("BIRTHYR", 0) or 0)
                mo = max(1, int(row.get("BIRTHMO", 1) or 1))
                if yr > 1800:
                    dob = pd.Timestamp(year=yr, month=mo, day=1)
            except Exception:
                dob = pd.NaT
            records[ptid] = {
                "ID"             : ptid,
                "gender"         : row.get("SEX",       np.nan),
                "race"           : row.get("RACE",       np.nan),
                "education"      : row.get("EDUC",       np.nan),
                "hispanic"       : row.get("HISPANIC",   np.nan),
                "marital_status" : row.get("MARISTAT",   np.nan),
                "handedness"     : row.get("HANDED",     np.nan),
                "birth_date"     : dob,
            }
    if not df_subcat.empty:
        sub = normalise_id(df_subcat.copy())
        safe_date(sub, "DATE")
        safe_date(sub, "DOB")
        first_sub = sub.sort_values("DATE").drop_duplicates("PTID", keep="first")
        for _, row in first_sub.iterrows():
            ptid = row["PTID"]
            rec  = records.setdefault(ptid, {"ID": ptid})
            if pd.notna(row.get("DOB")):
                rec["birth_date"] = row["DOB"]
            if "entry_age" not in rec:
                rec["entry_age"] = row.get("AGE", np.nan)
            if "gender" not in rec:
                rec["gender"] = row.get("GENDER", np.nan)
            if "education" not in rec or pd.isna(rec.get("education", np.nan)):
                rec["education"] = row.get("EDUC", np.nan)
    if not df_apoe.empty:
        apoe = normalise_id(df_apoe.copy())
        apoe = apoe.drop_duplicates("PTID")
        for _, row in apoe.iterrows():
            ptid = row["PTID"]
            rec  = records.setdefault(ptid, {"ID": ptid})
            geno = str(row.get("ApoE", "") or "").strip()
            rec["APOE_genotype"] = geno if geno and geno != "nan" else np.nan
            rec["APOE_e4_count"] = geno.count("4") if geno else np.nan
            rec["rs429358"]      = row.get("rs429358", np.nan)
            rec["rs7412"]        = row.get("rs7412",   np.nan)
    df = pd.DataFrame(list(records.values()))
    return df.set_index("ID") if not df.empty else df


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
        master["age_at_1st_diag"] = np.nan
    return master


# =============================================================================
# MAIN
# =============================================================================

def main() -> pd.DataFrame:
    _print("=" * 62)
    _print("ABVIB LONGITUDINAL INTEGRATION PIPELINE  v3")
    _print("=" * 62)
    _print("\nNOTE: MRI from the IDA search export (primary) + UDS_E1 flags (fallback).")
    _print("      No PET or CSF protein data available for ABVIB.")

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Load clinical tables ──────────────────────────────────────────────────
    _print("\nLoading raw tables...")
    dfs = {k: load(k) for k in PATHS}

    # ── Load IDA search export ───────────────────────────────────────────────
    _print("\nLoading IDA search export...")
    ida = load_ida()

    # ── Build domains ─────────────────────────────────────────────────────────
    _print("\nBuilding domain tables...")
    diag = build_diagnosis(dfs["dx"])

    # MRI: IDA primary + UDS_E1 fallback for uncovered subjects
    ref_ids = normalise_id(dfs["dx"]).get("PTID") if not dfs["dx"].empty else None
    mri_ida = build_mri_ida(ida, ref_ids)
    ida_subjects = set(mri_ida.index.tolist()) if not mri_ida.empty else set()
    mri_flags = build_mri_uds_flags(dfs["imaging"], ida_subjects)
    if not mri_flags.empty:
        _print(f"    [INFO] {len(mri_flags)} subjects not in the IDA export -> "
              f"MRI_T1w from UDS_E1 flags")

    # Combine IDA + flag-based MRI (whole subject x type blocks, IDA first)
    mri = ida_combine_blocks(mri_ida, mri_flags)

    clin  = build_all_clinical(dfs)
    demog = build_demographics(dfs["demog"], dfs["subcat"], dfs["apoe"])

    # ── Merge ─────────────────────────────────────────────────────────────────
    _print("\nMerging all domains on subject ID...")
    parts = [p for p in [diag, mri, clin, demog]
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
    master["PET_tracer_disponible"] = "not_available_in_ABVIB"
    master["CSF_available"]          = "not_available_in_ABVIB"

    master["Plasma_available"] = "not_available_in_ABVIB"

    LEAD = [
        "ID",
        "gender", "race", "education", "hispanic", "marital_status",
        "handedness", "birth_date", "age_at_1st_diag", "entry_age",
        "APOE_genotype", "APOE_e4_count", "rs429358", "rs7412",
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
        cnt = Counter()
        for v in master["MRI_type_disponible"].dropna():
            cnt.update(str(v).split("|"))
        if cnt:
            _print("\nMRI types (IDA + flags):")
            for t, n in cnt.most_common():
                _print(f"  {t:<30s} {n:>5,} subjects")

    _print("\nTest coverage (subjects with ≥1 valid score):")
    for name in list(dict.fromkeys(n for n, _ in variable_sources(DATASET, config=_CTX.config))):
        col = f"nb_of_{name}_visits"
        if col in master.columns:
            n = master[col].gt(0).sum()
            if n > 0:
                _print(f"  {name:<20s} {n:>5,} subjects")

    return master


# =============================================================================
