"""Fluid biomarker (plasma, CSF) and routine blood-lab series.

Used by the ADNI, A4 and PREVENT-AD pipelines. Output per marker follows the
scored-family schema (``all_{name}_date`` / ``all_{name}_score`` / ``1st_`` /
``Lst_`` / ``nb_of_{name}_visits``) plus ``nb_of_{name}_censored``,
``all_{name}_batch`` / ``all_{name}_rundate`` and ``{name}_lab`` / ``{name}_kit``.
Canonical names: ``Plasma_{Analyte}_{Platform}``, ``CSF_{Analyte}...``,
``Lab_{Analyte}``. Every ``all_*`` list is positionally aligned.
"""
from __future__ import annotations

import numpy as np
import os
import pandas as pd
import re
from functools import reduce

import logging as _logging
from .logging_utils import make_print as _make_print

log = _logging.getLogger(__name__)
_print = _make_print(log)

# ── value cleaning ────────────────────────────────────────────────────────────
def plasma_clean_numeric(s: pd.Series) -> pd.Series:
    """
    '<0.11' / '>1000' / '1,23' → numeric (censored values kept at the limit).
    Negative values are sentinel codes in ADNI/A4 (-1, -4, -5, -9 …) → NaN.
    """
    if s.dtype.kind in "if":
        out = s.astype(float)
    else:
        txt = (s.astype(str).str.strip()
                .str.replace(r"^[<>≤≥=]+", "", regex=True)
                .str.replace(",", ".", regex=False))
        out = pd.to_numeric(txt, errors="coerce")
    return out.where(out >= 0)


def plasma_sanitize(token, maxlen: int = 40) -> str:
    t = re.sub(r"[^0-9A-Za-z]+", "_", str(token)).strip("_")
    return t[:maxlen] if t else "NA"


_ANALYTE_RULES = [            # (regex on raw test/column name, canonical analyte)
    (r"p-?t(au)?_?217.*(ratio|%|npt)",            "pTau217_ratio"),
    (r"np-?t(au)?_?217",                          "npTau217"),
    (r"p-?t(au)?_?217",                           "pTau217"),
    (r"p-?t(au)?_?181",                           "pTau181"),
    (r"p-?t(au)?_?231",                           "pTau231"),
    (r"nt1",                                      "NT1tau"),
    (r"(a[bβ]|abeta|amyloid)\D*(1\D?)?42\D*(1\D?)?40|4240", "AB42_40"),
    (r"(a[bβ]|abeta|amyloid)\D*(1\D?)?42",        "AB42"),
    (r"(a[bβ]|abeta|amyloid)\D*(1\D?)?40",        "AB40"),
    (r"nf_?l|neurofilament",                      "NfL"),
    (r"gfap",                                     "GFAP"),
    (r"t-?tau|total.?tau|^tau$",                  "tTau"),
    (r"aps",                                      "APS"),
]


def plasma_canon_analyte(raw) -> str:
    """'PTAU217' / 'p-tau217' / 'Plasma pTau217 (Lilly)' → 'pTau217', etc."""
    r = str(raw).strip().lower()
    for rx, canon in _ANALYTE_RULES:
        if re.search(rx, r):
            return canon
    return plasma_sanitize(raw)


def _viscode_month(v) -> float:
    """'bl'/'sc' → 0 ; 'm24' → 24 ; anything else → 999 (sorted last)."""
    if pd.isna(v):
        return 999.0
    v = str(v).strip().lower()
    if v in ("bl", "sc", "scmri", "f", "init", "scv3", "baseline"):
        return 0.0
    m = re.search(r"m(\d+)", v)
    return float(m.group(1)) if m else 999.0


# ── file discovery ────────────────────────────────────────────────────────────
_PREFIX_RE = re.compile(r"^(copy of\s+)?(all_subjects_)?", re.I)


def plasma_find_tables(folder: str, stems, recursive: bool = True) -> dict:
    """
    Map each table stem (e.g. 'UPENN_PLASMA_FUJIREBIO_QUANTERIX') to the most
    recently modified CSV whose file name starts with that stem
    (ignoring 'Copy of ' / 'All_Subjects_' prefixes). A file is assigned to the
    LONGEST matching stem, so 'BLENNOWPLASMANFLLONG_…' never feeds
    'BLENNOWPLASMANFL', and 'BATEMANLAB_ADNI_Plasma_…' never feeds 'BATEMANLAB'.
    """
    if not os.path.isdir(folder):
        _print(f"  [SKIP] plasma folder not found → {folder}")
        return {}
    stems_sorted = sorted(stems, key=len, reverse=True)
    found = {}
    walker = os.walk(folder) if recursive else [(folder, [], os.listdir(folder))]
    for root, _, files in walker:
        for fn in files:
            if not fn.lower().endswith(".csv"):
                continue
            base = _PREFIX_RE.sub("", fn).lower()
            for stem in stems_sorted:
                st = stem.lower()
                if base.startswith(st) and (len(base) == len(st) or
                                            base[len(st)] in "_.( -"):
                    found.setdefault(stem, []).append(os.path.join(root, fn))
                    break
    out = {}
    for stem, paths in found.items():
        paths.sort(key=os.path.getmtime, reverse=True)
        out[stem] = paths[0]
        if len(paths) > 1:
            _print(f"    [INFO] {stem}: {len(paths)} files, using newest → "
                  f"{os.path.basename(paths[0])}")
    return out


def plasma_read_csv(path: str) -> pd.DataFrame:
    for enc in ("utf-8", "latin-1", "cp1252"):
        try:
            return pd.read_csv(path, low_memory=False, encoding=enc)
        except UnicodeDecodeError:
            continue
        except Exception as e:
            _print(f"    [ERROR] {os.path.basename(path)}: {e}")
            return pd.DataFrame()
    return pd.DataFrame()


# ── core longitudinal builder ─────────────────────────────────────────────────
# Columns describing WHO measured and WITH WHAT (auto-detected when present)
META_LAB_COLS   = ("NAM", "LBNAM", "PERFORMINGLAB", "LAB", "LABORATORY")
META_KIT_COLS   = ("KIT", "ASSAYPLATFORM", "ASSAYVERSION", "IMMUNOASSAY", "ASSAY",
                   "LBMETHOD", "METHOD", "PLATFORM", "INSTRUMENT")
META_BATCH_COLS = ("BATCH", "BATCH_N", "RUN", "RUN_NUM", "RUN_NO", "RUN_NUMBER",
                   "RBARCODE", "PLATE", "PLATEID")
META_RUNDATE_COLS = ("RUNDATE", "MS_RUN_DATE", "ASSAY_START_DATE", "DATE_OF_ANALYSIS")


def _first_col(df: pd.DataFrame, candidates, exclude=()):
    cols = {c.upper(): c for c in df.columns}
    for c in candidates:
        if c in cols and cols[c] not in exclude:
            return cols[c]
    return None


def plasma_series(df: pd.DataFrame, id_col: str, value_col: str, name: str,
                  time_col: str | None = None,
                  visit_col: str | None = None,
                  method_col: str | None = None,
                  lab: str | None = None,
                  kit: str | None = None,
                  exclude=()) -> pd.DataFrame:
    """
    One canonical marker → one block of per-subject columns.
    time_col  : collection date (parsed) or numeric time axis (e.g. EXAMDAY)
    visit_col : VISCODE-like label, used for ordering/labels when time is missing
    lab / kit : static laboratory / assay-kit description (from the DATADIC);
                columns present in the table (NAM, KIT, ASSAYPLATFORM, …) win.
    Replicates at the same subject/visit are averaged.

    Columns (all_* lists are POSITIONALLY ALIGNED — one entry per measurement):
      all_{n}_date, all_{n}_score, 1st/Lst_{n}_score, 1st/Lst_{n}_date,
      nb_of_{n}_visits, nb_of_{n}_censored  (values reported as <x / >x),
      all_{n}_batch, all_{n}_rundate        (when the table has them),
      {n}_lab, {n}_kit
    """
    if df.empty or value_col not in df.columns or id_col not in df.columns:
        return pd.DataFrame()

    raw = df[value_col]
    d = pd.DataFrame({"ID": df[id_col], "_v": plasma_clean_numeric(raw)})
    d["_cens"] = raw.astype(str).str.strip().str.match(r"^[<>≤≥]").astype(int)
    if time_col and time_col in df.columns:
        t = df[time_col]
        if t.dtype.kind not in "iufM":
            t = pd.to_datetime(t, errors="coerce")
        d["_t"] = t
    else:
        d["_t"] = pd.NaT
    d["_vis"] = df[visit_col].astype(str) if visit_col and visit_col in df.columns else np.nan
    d["_vm"]  = d["_vis"].map(_viscode_month)

    excl      = (id_col, value_col) + tuple(exclude)
    lab_col   = _first_col(df, META_LAB_COLS, exclude=excl)
    kit_col   = (method_col if method_col in df.columns
                 else _first_col(df, META_KIT_COLS, exclude=excl))
    batch_col = _first_col(df, META_BATCH_COLS)
    run_col   = _first_col(df, META_RUNDATE_COLS)
    for tag, col in (("_lab", lab_col), ("_kit", kit_col), ("_batch", batch_col), ("_run", run_col)):
        d[tag] = df[col].astype(str).where(df[col].notna()) if col else np.nan
    if run_col:
        d["_run"] = pd.to_datetime(d["_run"], errors="coerce").dt.strftime("%Y-%m-%d")

    d = d.dropna(subset=["ID", "_v"])
    if d.empty:
        return pd.DataFrame()

    def _lab(r):
        t = r["_t"]
        if pd.notna(t):
            return str(t.date()) if hasattr(t, "date") else str(t)
        return r["_vis"] if pd.notna(r["_vis"]) else "NA"
    d["_label"] = d.apply(_lab, axis=1)

    _j = lambda s: "/".join(sorted(set(s.dropna().astype(str)))) or "NA"
    d = (d.groupby(["ID", "_label"], dropna=False)
           .agg(_v=("_v", "mean"), _t=("_t", "first"), _vm=("_vm", "first"),
                _cens=("_cens", "max"), _lab=("_lab", _j), _kit=("_kit", _j),
                _batch=("_batch", _j), _run=("_run", _j))
           .reset_index())
    d = d.sort_values(["ID", "_t", "_vm"], na_position="last")

    rows = []
    for subj, g in d.groupby("ID", sort=False):
        vals, labs = g["_v"].round(4).tolist(), g["_label"].tolist()
        row = {
            "ID"                    : subj,
            f"all_{name}_date"      : "|".join(map(str, labs)),
            f"all_{name}_score"     : "|".join(map(str, vals)),
            f"1st_{name}_score"     : vals[0],
            f"1st_{name}_date"      : labs[0],
            f"Lst_{name}_score"     : vals[-1],
            f"Lst_{name}_date"      : labs[-1],
            f"nb_of_{name}_visits"  : len(vals),
            f"nb_of_{name}_censored": int(g["_cens"].sum()),
        }
        if batch_col:
            row[f"all_{name}_batch"] = "|".join(g["_batch"])
        if run_col:
            row[f"all_{name}_rundate"] = "|".join(g["_run"])
        labs_found = sorted({x for x in g["_lab"] if x != "NA"})
        kits_found = sorted({x for x in g["_kit"] if x != "NA"})
        def _combine(static, found):
            if found and static:
                return f"{static} [{'|'.join(found)}]"
            return "|".join(found) if found else (static or np.nan)
        row[f"{name}_lab"] = _combine(lab, labs_found)
        row[f"{name}_kit"] = _combine(kit, kits_found)
        rows.append(row)
    return pd.DataFrame(rows).set_index("ID")


# ── spec-driven extraction ────────────────────────────────────────────────────
_LONG_TEST_COLS  = ("TESTCD", "LBTESTCD", "TESTNAME", "LBTEST", "TEST",
                    "ANALYTE", "ASSAY", "PLASMA_BIOMARKER", "BIOMARKER")
_LONG_VALUE_COLS = ("ORRES", "LBORRESN", "LBSTRESN", "LBORRES", "LBSTRESC",
                    "RESULT", "CONCENTRATION", "TESTVALUE", "VALUE")


def plasma_detect_long(df: pd.DataFrame):
    cols = {c.upper(): c for c in df.columns}
    tcol = next((cols[c] for c in _LONG_TEST_COLS if c in cols), None)
    vcol = next((cols[c] for c in _LONG_VALUE_COLS if c in cols), None)
    return (tcol, vcol) if tcol and vcol else (None, None)


def plasma_from_spec(df: pd.DataFrame, spec: dict, id_col: str,
                     time_col: str | None, visit_col: str | None) -> list:
    """
    spec keys:
      markers   : [(raw_col, canonical_name[, lab[, kit]]), …]   wide tables
      lab / kit : static laboratory / kit for the whole table (DATADIC)
      long      : True → test_col / value_col / name_tpl ('{test}' placeholder)
      test_map  : long tables only - {raw test name: canonical name}; when given,
                  only the listed tests are read and name_tpl is not used
      split_col : run every marker separately per value of this column
      method    : column describing assay / batch / platform (optional)
      filter    : {col: regex} — keep rows whose col matches (optional)
    """
    parts = []
    if df.empty:
        return parts
    for col, rx in spec.get("filter", {}).items():
        if col in df.columns:
            df = df[df[col].astype(str).str.contains(rx, case=False, regex=True, na=False)]
    method = spec.get("method")

    if spec.get("long"):
        tcol = spec.get("test_col") or plasma_detect_long(df)[0]
        vcol = spec.get("value_col") or plasma_detect_long(df)[1]
        if not tcol or not vcol or tcol not in df.columns or vcol not in df.columns:
            _print(f"    [WARN] long-format columns not found ({tcol}, {vcol})")
            return parts
        test_map = spec.get("test_map") or {}
        for test, g in df.groupby(tcol):
            if test_map:
                if test not in test_map:
                    continue
                name = test_map[test]
            else:
                name = spec["name_tpl"].format(test=plasma_canon_analyte(test))
            parts.append(plasma_series(g, id_col, vcol, name, time_col, visit_col, method,
                                       spec.get("lab"), spec.get("kit"), exclude=(tcol,)))
        return parts

    groups = [(None, df)]
    if spec.get("split_col") and spec["split_col"] in df.columns:
        groups = list(df.groupby(spec["split_col"]))
    for split_val, g in groups:
        for mk in spec.get("markers", []):
            raw_col, canonical = mk[0], mk[1]
            m_lab = mk[2] if len(mk) > 2 else spec.get("lab")
            m_kit = mk[3] if len(mk) > 3 else spec.get("kit")
            if raw_col not in g.columns:
                continue
            name = canonical if split_val is None else f"{canonical}_{plasma_sanitize(split_val, 20)}"
            parts.append(plasma_series(g, id_col, raw_col, name, time_col, visit_col, method,
                                       m_lab, m_kit))
    return [p for p in parts if not p.empty]


def plasma_join(parts: list) -> pd.DataFrame:
    parts = [p for p in parts if isinstance(p, pd.DataFrame) and not p.empty]
    if not parts:
        return pd.DataFrame()
    out = reduce(lambda a, b: a.join(b, how="outer"), parts)
    out.index.name = "ID"
    return out


def plasma_add_summary(df: pd.DataFrame, prefix: str = "Plasma") -> pd.DataFrame:
    """Add {prefix}_panels_disponible / {prefix}_analytes_disponible."""
    if df.empty:
        return df
    nb_cols = [c for c in df.columns
               if c.startswith(f"nb_of_{prefix}_") and c.endswith("_visits")]
    if not nb_cols:
        return df
    names = [c[len("nb_of_"):-len("_visits")] for c in nb_cols]
    analyte = {n: n.split("_")[1] if n.count("_") >= 1 else n for n in names}
    # keep composite analyte tokens readable (AB42_40, pTau217_AB42 …)
    for n in names:
        m = re.match(rf"{prefix}_(AB42_40|pTau217_AB42|pTau217_ratio|AB_ratio)", n)
        if m:
            analyte[n] = m.group(1)
    nb = df[nb_cols].fillna(0).gt(0)
    df[f"{prefix}_panels_disponible"] = nb.apply(
        lambda r: "|".join(n for n, c in zip(names, nb_cols) if r[c]) or np.nan, axis=1)
    df[f"{prefix}_analytes_disponible"] = nb.apply(
        lambda r: "|".join(sorted({analyte[n] for n, c in zip(names, nb_cols) if r[c]})) or np.nan,
        axis=1)
    return df


def specs_from_config(specs: dict) -> dict:
    """Fluid-table specs read from YAML: marker lists become tuples."""
    out = {}
    for stem, spec in (specs or {}).items():
        spec = dict(spec)
        if "markers" in spec:
            spec["markers"] = [tuple(m) for m in spec["markers"]]
        out[stem] = spec
    return out
