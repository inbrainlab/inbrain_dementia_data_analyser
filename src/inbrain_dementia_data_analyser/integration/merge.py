"""Merge of the per-dataset tables into one master table.

Column renames, demographic codebooks and the diagnosis map are read from
``harmonization.yaml``, ``demographics.yaml`` and ``diagnosis.yaml``
(:func:`configure`); :func:`merge_cohorts` runs the merge.
"""
from __future__ import annotations

import ast
import gc
import numpy as np
import os
import pandas as pd
import re

from ..config import ProjectConfig, load_config
from ..provenance import write_provenance

from ..config import ProjectConfig, load_config

import logging as _logging
from ..logging_utils import make_print as _make_print
from ..logging_utils import quiet_pandas_warnings

log = _logging.getLogger(__name__)
_print = _make_print(log)

_unmapped_demog: dict = {}     # (cohort, column) -> set of raw codes without a rule
_unmapped_diag: set = set()    # diagnosis tokens without a rule

# Filled by configure() from harmonization.yaml / demographics.yaml / diagnosis.yaml
GLOBAL_RENAME: dict = {}
COLUMN_RENAME: dict = {}
FAMILY_RENAME: dict = {}
AGE_PRIORITY: list = []
_SEX_STD: dict = {}
SEX_CODES: dict = {}
RACE_CODES: dict = {}
RACE_TEXT: dict = {}
HISPANIC_CODES: dict = {}
_HISP_01: dict = {}
HAND_CODES: dict = {}
_HAND_TEXT: dict = {}
EDU_SENTINELS: set = set()
EDU_MAX_YEARS: float = 40
NACC_OLD_APOE_FIX: dict = {}
_NACC_APOE_CODE: dict = {}
DIAG_MAP: dict = {}
DIAG_LABEL_COLS: list = []


def configure(cfg: ProjectConfig) -> None:
    """Load renames, codebooks and the diagnosis map from a configuration."""
    global GLOBAL_RENAME, COLUMN_RENAME, FAMILY_RENAME, AGE_PRIORITY, _SEX_STD, SEX_CODES
    global RACE_CODES, RACE_TEXT, HISPANIC_CODES, _HISP_01, HAND_CODES, _HAND_TEXT
    global EDU_SENTINELS, EDU_MAX_YEARS, NACC_OLD_APOE_FIX, _NACC_APOE_CODE, DIAG_MAP, DIAG_LABEL_COLS
    h, d, g = cfg.harmonization, cfg.demographics, cfg.diagnosis
    GLOBAL_RENAME = dict(h.get("global_rename") or {})
    COLUMN_RENAME = {k: dict(v) for k, v in (h.get("column_rename") or {}).items()}
    FAMILY_RENAME = {k: dict(v) for k, v in (h.get("family_rename") or {}).items()}
    AGE_PRIORITY = list(h.get("age_priority") or [])
    _SEX_STD = dict(d.get("sex_standard") or {})
    SEX_CODES = {k: dict(v) for k, v in (d.get("sex_codes") or {}).items()}
    RACE_CODES = {k: dict(v) for k, v in (d.get("race_codes") or {}).items()}
    RACE_TEXT = dict(d.get("race_text") or {})
    HISPANIC_CODES = {k: dict(v) for k, v in (d.get("hispanic_codes") or {}).items()}
    _HISP_01 = dict(d.get("hispanic_default") or {})
    HAND_CODES = {k: dict(v) for k, v in (d.get("handedness_codes") or {}).items()}
    _HAND_TEXT = dict(d.get("handedness_text") or {})
    EDU_SENTINELS = set(d.get("education_sentinels") or [])
    EDU_MAX_YEARS = d.get("education_max_years", 40)
    NACC_OLD_APOE_FIX = dict(d.get("nacc_old_apoe_fix") or {})
    _NACC_APOE_CODE = dict(d.get("nacc_apoe_codes") or {})
    DIAG_MAP = dict(g.get("diag_map") or {})
    DIAG_LABEL_COLS = list(g.get("label_columns") or [])

ID_COL_RX = re.compile(r"(^ID$|_IDs$|_image_ids$|_ID$)")


def _code_key(v):
    """Raw cell -> lookup key: 1.0 -> '1', ' White ' -> 'white', NaN -> None."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    s = str(v).strip()
    if s == "" or s.lower() in ("nan", "none", "nat"):
        return None
    s = re.sub(r"\.0+$", "", s)
    return s.lower()


def _map_codes(series, table, cohort, col):
    """Map through `table` once per DISTINCT raw value; unknown codes -> NaN + report."""
    lut = {}
    for raw in pd.unique(series):
        k = _code_key(raw)
        if k is None:
            lut[raw] = np.nan
        elif k in table:
            lut[raw] = table[k] if table[k] is not None else np.nan
        else:
            lut[raw] = np.nan
            _unmapped_demog.setdefault((cohort, col), set()).add(str(raw))
    return series.map(lut)


def _race_one(raw, table, cohort):
    """One race cell (may hold several codes: '1|5', '1:4:5') -> (label, detail)."""
    k = _code_key(raw)
    if k is None:
        return (np.nan, np.nan)
    if k in RACE_TEXT and k not in table:
        return (RACE_TEXT[k] if RACE_TEXT[k] else np.nan, np.nan)
    labels = []
    for tok in re.split(r"[|:;,]", k):
        tok = re.sub(r"\.0+$", "", tok.strip())
        if not tok:
            continue
        if tok in table:
            lab = table[tok]
        elif tok in RACE_TEXT:
            lab = RACE_TEXT[tok]
        else:
            _unmapped_demog.setdefault((cohort, "race"), set()).add(str(raw))
            lab = None
        if lab and lab not in labels:
            labels.append(lab)
    if not labels:
        return (np.nan, np.nan)
    if len(labels) == 1:
        return (labels[0], np.nan)
    return ("Multiracial", "+".join(sorted(labels)))


def canon_apoe(v):
    """'3/4' · 'E3/E4' · 'e3,e4' · '4 3' · 34.0 -> 'e3/e4' (lower allele first)."""
    k = _code_key(v)
    if k is None:
        return np.nan
    m = re.fullmatch(r"([234])([234])", k)                     # OASIS: 34
    al = list(m.groups()) if m else re.findall(r"(?<![\d-])([234])(?!\d)", k)
    if len(al) != 2:
        return np.nan                                           # e.g. AIBL '-4/-4'
    a, b = sorted(al)
    return f"e{a}/e{b}"


def harmonise_demographics(df: pd.DataFrame, cohort: str) -> pd.DataFrame:
    """Rewrite sex / race / hispanic / education / handedness / APOE of one
    cohort with the shared vocabulary (see 2d). Originals are kept in *_raw."""
    def keep_raw(col):
        if col in df.columns and f"{col}_raw" not in df.columns:
            df[f"{col}_raw"] = df[col]

    # ── sex ──────────────────────────────────────────────────────────────────
    if "gender" in df.columns:
        keep_raw("gender")
        df["gender"] = _map_codes(df["gender"], SEX_CODES.get(cohort, _SEX_STD), cohort, "gender")
    if "IDA_sex" in df.columns:                      # sex recorded in the IDA export
        ida_sex = df["IDA_sex"].where(df["IDA_sex"].isin(["M", "F"]))
        df["gender"] = df["gender"].where(df["gender"].notna(), ida_sex) \
            if "gender" in df.columns else ida_sex

    # ── race ─────────────────────────────────────────────────────────────────
    if "race" in df.columns or "race_label" in df.columns:
        src = df["race"] if "race" in df.columns else df["race_label"]
        if "race" in df.columns:
            keep_raw("race")
        table = RACE_CODES.get(cohort, {})
        lut = {raw: _race_one(raw, table, cohort) for raw in pd.unique(src.dropna())}
        df["race"] = src.map({k: v[0] for k, v in lut.items()})
        detail = src.map({k: v[1] for k, v in lut.items()})
        if "race_label" in df.columns:               # text label fills what the code left empty
            txt = df["race_label"].map(lambda x: RACE_TEXT.get(_code_key(x) or "", None))
            df["race"] = df["race"].where(df["race"].notna(), txt)
        if detail.notna().any():
            df["race_detail"] = detail

    # ── hispanic / latino ────────────────────────────────────────────────────
    if "hispanic" in df.columns:
        keep_raw("hispanic")
        df["hispanic"] = pd.to_numeric(
            _map_codes(df["hispanic"], HISPANIC_CODES.get(cohort, _HISP_01), cohort, "hispanic"),
            errors="coerce")

    # ── education (years) ────────────────────────────────────────────────────
    if "education" in df.columns:
        edu = pd.to_numeric(df["education"], errors="coerce")
        bad = edu.isin(EDU_SENTINELS) | (edu < 0) | (edu > EDU_MAX_YEARS)
        if bad.any():
            keep_raw("education")
        df["education"] = edu.where(~bad)

    # ── handedness ───────────────────────────────────────────────────────────
    if "handedness" in df.columns:
        keep_raw("handedness")
        df["handedness"] = _map_codes(df["handedness"],
                                      {**_HAND_TEXT, **HAND_CODES.get(cohort, {})},
                                      cohort, "handedness")

    # ── APOE ─────────────────────────────────────────────────────────────────
    if "APOE_genotype" in df.columns:
        keep_raw("APOE_genotype")
        lut = {raw: canon_apoe(raw) for raw in pd.unique(df["APOE_genotype"])}
        geno = df["APOE_genotype"].map(lut)
        if cohort == "NACC":
            if "APOE_code" in df.columns:            # pipeline >= 2026-11-16: raw NACCAPOE
                code = _map_codes(df["APOE_code"], _NACC_APOE_CODE, cohort, "APOE_code")
                geno = code.where(code.notna(), geno)
            elif geno.notna().any() and (geno == "e2/e2").sum() > (geno == "e3/e3").sum():
                # e2/e2 more frequent than e3/e3 is impossible in a real cohort:
                # this file still carries the old wrong NACCAPOE labels
                geno = geno.map(NACC_OLD_APOE_FIX)
                _print("  [WARN] NACC: APOE labels of an OLD NACC_longitudinal.csv "
                      "(wrong NACCAPOE map) were corrected here; re-run the NACC "
                      "pipeline to fix the file itself.")
        df["APOE_genotype"] = geno
        n4 = geno.str.count("e4")
        old = pd.to_numeric(df["APOE_e4_count"], errors="coerce") \
            if "APOE_e4_count" in df.columns else pd.Series(np.nan, index=df.index)
        df["APOE_e4_count"] = n4.where(geno.notna(), old if cohort != "NACC" else np.nan)
    if "APOE_e4_count" in df.columns:
        cnt = pd.to_numeric(df["APOE_e4_count"], errors="coerce")
        any_old = pd.to_numeric(df["APOE_e4_any"], errors="coerce") \
            if "APOE_e4_any" in df.columns else pd.Series(np.nan, index=df.index)
        df["APOE_e4_any"] = (cnt > 0).astype(float).where(cnt.notna(), any_old)
    return df


def _norm_token(tok: str) -> str:
    """Map one raw diagnosis token to its canonical label."""
    key = str(tok).strip()
    up = key.upper()
    if up in DIAG_MAP:
        return DIAG_MAP[up]
    if key:                       # remember anything we did not recognise
        _unmapped_diag.add(key)
    return key                    # keep unknown tokens as-is


def normalize_diag(value):
    """
    Normalise a diagnosis cell. Handles:
      - plain labels:        "CON"            -> "CN"
      - compound labels:     "AD+VaD"         -> "AD+VaD"
      - list-strings:        "['CON','MCI']"  -> "['CN', 'MCI']"
      - pipe lists (all_diag / Diag_total as written by the pipelines):
                             "CON|CON|MCI"    -> "CN|CN|MCI"
    NaN/None passes through unchanged.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return np.nan
    s = str(value).strip()
    if s == "" or s.lower() in ("nan", "none", "nat"):
        return np.nan

    if s.startswith("[") and s.endswith("]"):
        try:
            items = ast.literal_eval(s)
            if isinstance(items, (list, tuple)):
                out = ["+".join(_norm_token(t) for t in str(it).split("+"))
                       for it in items]
                return str(out)
        except Exception:
            pass

    # pipe-separated visit lists: normalise every visit label on its own
    # (before 2026-11-16 the whole "CON|CON|MCI" string was looked up as ONE
    # token, so all_diag kept raw labels such as CON / DEMENTIA)
    return "|".join("+".join(_norm_token(t) for t in part.split("+"))
                    for part in s.split("|"))


def load_csv(path: str) -> pd.DataFrame | None:
    if not os.path.exists(path):
        _print(f"  [WARN] missing file, skipping: {path}")
        return None
    attempts = [
        dict(),
        dict(engine="python"),
        dict(encoding="utf-8"),
        dict(encoding="latin-1"),
        dict(encoding="cp1252", engine="python"),
        dict(encoding="latin-1", engine="python", on_bad_lines="skip"),
    ]
    last_err = None
    for kw in attempts:
        try:
            hdr = pd.read_csv(path, nrows=0, **{k: v for k, v in kw.items()
                                                if k != "on_bad_lines"}).columns
            as_text = {c: str for c in hdr if ID_COL_RX.search(c)}
            if kw.get("engine") == "python":
                return pd.read_csv(path, dtype=as_text, **kw)
            return pd.read_csv(path, low_memory=False, dtype=as_text, **kw)
        except Exception as e:
            last_err = e
    _print(f"  [ERROR] could not read {path}: {last_err}")
    return None


def apply_family_rename(df: pd.DataFrame, cohort: str) -> pd.DataFrame:
    """Rename whole test blocks per FAMILY_RENAME[cohort].
    Safe with already-canonical columns:
      * a stem that already carries the canonical prefix is left alone
        ('T2w' -> 'MRI_T2w' never turns 'all_MRI_T2w_FLAIR_date' into
        'all_MRI_MRI_T2w_FLAIR_date');
      * a rename that would create a column name that already exists is
        skipped and reported (duplicate names break the outer concat)."""
    mapping = FAMILY_RENAME.get(cohort, {})
    if not mapping:
        return df
    rename, taken, clash = {}, set(df.columns), []
    for raw, canon in mapping.items():
        rx = re.compile(rf"(^|_){re.escape(raw)}(_|$)")
        prefix = canon[:-len(raw)] if canon.endswith(raw) and canon != raw else None
        for col in df.columns:
            # match all_/1st_/Lst_/nb_of_ <stem> ... patterns
            m = rx.search(col)
            if not m:
                continue
            stem_at = m.start() + len(m.group(1))
            if prefix and col[:stem_at].endswith(prefix):
                continue                      # already canonical
            new = col[:stem_at] + canon + col[stem_at + len(raw):]
            if new in taken and new != col:
                clash.append(col)
                continue
            taken.add(new)
            rename[col] = new
    if clash:
        _print(f"  [WARN] {cohort}: {len(clash)} columns NOT renamed (canonical name "
              f"already present): {clash[:6]}{' …' if len(clash) > 6 else ''}")
    return df.rename(columns=rename)


def add_age_at_baseline(df: pd.DataFrame) -> pd.DataFrame:
    """Coalesce the first available baseline-age column into age_at_baseline."""
    present = [c for c in AGE_PRIORITY if c in df.columns]
    if not present:
        return df
    base = df[present[0]].copy()
    for c in present[1:]:
        base = base.combine_first(df[c])
    df["age_at_baseline"] = base
    return df

@quiet_pandas_warnings
def merge_cohorts(root_path: str | os.PathLike, datasets: list[str] | None = None, *,
                  out_dir: str | os.PathLike | None = None,
                  config: ProjectConfig | None = None) -> pd.DataFrame:
    """Stack the per-dataset tables into one harmonised master table.

    For every dataset (in ``integration.yaml -> merge_order``) the table
    written by :func:`~inbrain_dementia_data_analyser.preprocess.preprocess_dataset`
    is read, its columns are renamed to the canonical names
    (``harmonization.yaml``), demographics are recoded (``demographics.yaml``),
    ``Databank`` and ``uid`` (``<Databank>_<ID>``) are added, and the tables are
    concatenated on identical column names. Diagnosis labels are then mapped
    with ``diagnosis.yaml``.

    Writes ``<root>/Integrated_ALL/MASTER_all_cohorts.csv`` and
    ``MASTER_column_coverage.csv`` (non-missing values per column and dataset).

    Parameters
    ----------
    root_path:
        Root folder of the datasets.
    datasets:
        Datasets to include (default: all, in ``merge_order``). Missing
        per-dataset tables are skipped with a warning.
    out_dir:
        Output folder (default ``<root>/Integrated_ALL``).
    config:
        Project configuration.

    Returns
    -------
    pandas.DataFrame
        The master table (one row per subject).
    """
    cfg = config or load_config()
    configure(cfg)
    _unmapped_demog.clear()
    _unmapped_diag.clear()
    order = [cfg.dataset(n).name for n in cfg.integration.get("merge_order", list(cfg.datasets))]
    wanted = [cfg.dataset(n).name for n in (datasets or order)]
    names = [n for n in order if n in wanted] + [n for n in wanted if n not in order]
    COHORTS = {n: str(cfg.dataset(n).output_path(root_path)) for n in names}
    OUT_DIR = str(out_dir or os.path.join(root_path, cfg.integration.get("output_folder", "Integrated_ALL")))
    OUT_MASTER = os.path.join(OUT_DIR, "MASTER_all_cohorts.csv")
    OUT_REPORT = os.path.join(OUT_DIR, "MASTER_column_coverage.csv")
    os.makedirs(OUT_DIR, exist_ok=True)

    frames = []
    _print("=" * 70)
    _print("LOADING & HARMONISING COHORTS")
    _print("=" * 70)
    
    for cohort, path in COHORTS.items():
        df = load_csv(path)
        if df is None:
            continue
    
        # 5a. column-name harmonisation
        df = df.rename(columns={k: v for k, v in GLOBAL_RENAME.items() if k in df.columns})
        df = df.rename(columns={k: v for k, v in COLUMN_RENAME.get(cohort, {}).items()
                                if k in df.columns})
        df = apply_family_rename(df, cohort)
        df = harmonise_demographics(df, cohort)      # sex / race / education / APOE … (2d)
        df = add_age_at_baseline(df)
    
        # 5b. bookkeeping columns
        if "ID" not in df.columns:
            _print(f"  [WARN] {cohort}: no 'ID' column found — using row index as ID")
            df["ID"] = df.index.astype(str)
        df["ID"] = df["ID"].astype(str).str.strip()
        df["Databank"] = cohort
        df["uid"] = cohort + "_" + df["ID"]
    
        _print(f"  {cohort:<10} rows={len(df):>6}  cols={df.shape[1]:>4}")
        frames.append(df)
    
    if not frames:
        raise FileNotFoundError("No per-dataset table could be loaded - run preprocess_dataset first "
                                f"(looked for: {list(COHORTS.values())})")
    
    # -----------------------------------------------------------------------------
    # 6. OUTER CONCAT  (join on identical column names; union of all columns)
    # -----------------------------------------------------------------------------
    master = pd.concat(frames, axis=0, ignore_index=True, sort=False)
    del frames                       # the per-cohort copies are no longer needed:
    gc.collect()                     # keeping them doubled the RAM held by this cell
    _print("\n" + "=" * 70)
    _print(f"MERGED: {len(master)} subjects  x  {master.shape[1]} columns")
    _print("=" * 70)
    
    # -----------------------------------------------------------------------------
    # 7. NORMALISE DIAGNOSIS LABELS  (applied to the combined table)
    # -----------------------------------------------------------------------------
    for col in DIAG_LABEL_COLS:
        if col in master.columns:
            master[col] = master[col].apply(normalize_diag)
    
    # -----------------------------------------------------------------------------
    # 8. COLUMN ORDER  — identifiers + key demographics + diagnosis first
    # -----------------------------------------------------------------------------
    PRIORITY = [
        "uid", "Databank", "ID",
        "gender", "race", "hispanic", "education", "education_level",
        "race_detail", "marital_status", "retired", "employment", "handedness",
        "living_situation", "residence_type",
        "age_at_baseline", "age_at_1st_diag", "birth_date",
        "APOE_genotype", "APOE_e4_count", "APOE_e4_any",
        "1st_diag", "Lst_diag", "1st_diag_date", "Lst_diag_date",
        "all_diag", "all_diag_date", "Diag_total",
        "1st_stage", "Lst_stage", "max_stage", "all_stage",
        "1st_etiology", "Lst_etiology", "all_etiology",
        "epilepsy_any", "nb_of_epilepsy_visits",
        "1st_diag_legacy", "Lst_diag_legacy", "all_diag_legacy",
        "nb_of_visits", "nb_of_change", "nb_of_label_change",
        "1st_time_of_change", "Lst_time_of_change", "all_time_of_change",
    ]
    front = [c for c in PRIORITY if c in master.columns]
    rest  = sorted(c for c in master.columns if c not in front)
    master = master[front + rest]
    
    # -----------------------------------------------------------------------------
    # 9. SAVE + REPORTS
    # -----------------------------------------------------------------------------
    master.to_csv(OUT_MASTER, index=False)
    _print(f"\n[saved] master  -> {OUT_MASTER}")
    
    # 9a. per-column non-null coverage by cohort
    # (boolean mask + groupby-sum: no full copy of the master, unlike assign/apply)
    cov = master.notna().groupby(master["Databank"]).sum().T
    cov.loc["_n"] = master.groupby("Databank").size()
    cov.to_csv(OUT_REPORT)
    _print(f"[saved] coverage report -> {OUT_REPORT}")
    
    # 9b. console diagnostics
    
    # harmonised demographics (subjects per category, by cohort)
    for _col in ("gender", "race", "hispanic", "handedness", "APOE_genotype"):
        if _col in master.columns:
            _print(f"\nHarmonised {_col} (subjects):")
            _print(pd.crosstab(master["Databank"], master[_col].fillna("(missing)")).to_string())
    if "education" in master.columns:
        _print("\nEducation in years (after removing sentinels):")
        _print(master.groupby("Databank")["education"].agg(["count", "min", "median", "max"])
              .round(1).to_string())
    if _unmapped_demog:
        _print("\n[ACTION] Demographic codes WITHOUT a rule (set to missing — add them to "
              "the codebooks in section 2d):")
        for (_c, _k), _v in sorted(_unmapped_demog.items()):
            _print(f"   - {_c:<10} {_k:<14} {sorted(_v)[:12]}")
    
    if "1st_stage" in master.columns:
        _print("\nCognitive stage at the first visit with a stage (subjects):")
        _print(pd.crosstab(master["Databank"], master["1st_stage"].fillna("(missing)")).to_string())
    if "1st_diag" in master.columns:
        _print("\nHarmonised 1st_diag (primary label only):")
        prim = master["1st_diag"].astype(str).str.split("+").str[0]
        _print(prim.value_counts(dropna=False).to_string())
    
    # 9c. seizure / glucose / BMI coverage (subjects with >= 1 value)  — added 2026-11-08
    EXT_MEASURES = {
        "SEIZURES": "nb_of_SEIZURES_visits", "EPILEP": "nb_of_EPILEP_visits",
        "Antiepileptic_med": "nb_of_Antiepileptic_med_visits",
        "SEIZURE_text_flag": "SEIZURE_text_flag", "AED_text_flag": "AED_text_flag",
        "Neuro_history": "nb_of_Neuro_history_visits",
        "Neuro_history_ever": "Neuro_history_ever",
        "Glucose(ABVIB)": "nb_of_Glucose_visits", "Lab_Glucose": "nb_of_Lab_Glucose_visits",
        "Lab_Glucose_preFDG": "nb_of_Lab_Glucose_preFDG_visits",
        "Lab_Glucose_NMR": "nb_of_Lab_Glucose_NMR_visits",
        "Lab_Glucose_AIBL": "nb_of_Lab_Glucose_AIBL_visits",
        "HbA1c": "nb_of_HbA1c_visits", "Fasting": "nb_of_Fasting_visits",
        "DIABETES": "nb_of_DIABETES_visits", "Diabetes": "nb_of_Diabetes_visits",
        "Diabetes_treatment": "nb_of_Diabetes_treatment_visits",
        "Diabetes_med": "nb_of_Diabetes_med_visits",
        "BMI": "nb_of_BMI_visits", "CAIDE_BMI": "nb_of_CAIDE_BMI_visits",
    }
    _ext = {}
    for lab, col in EXT_MEASURES.items():
        if col in master.columns:
            v = pd.to_numeric(master[col], errors="coerce")
            ok = v.gt(0) if col.startswith("nb_of_") else v.notna()
            _ext[lab] = ok.groupby(master["Databank"]).sum().astype(int)
    if _ext:
        _print("\nSeizure / glucose / BMI coverage (subjects with >= 1 value; *_flag/_ever = subjects with any record):")
        _print(pd.DataFrame(_ext).T.fillna(0).astype(int).to_string())
    
    if _unmapped_diag:
        _print("\n[ACTION] Diagnosis tokens NOT in DIAG_MAP (kept as-is — review/extend):")
        for t in sorted(_unmapped_diag):
            _print(f"   - {t!r}")
    else:
        _print("\nAll diagnosis tokens were mapped successfully.")
    
    _print("\nDone.")

    write_provenance(OUT_MASTER, step="merge", config=cfg,
                     inputs=[p for p in COHORTS.values() if os.path.exists(p)], table=master,
                     parameters={"root_path": str(root_path), "datasets": names})
    return master
