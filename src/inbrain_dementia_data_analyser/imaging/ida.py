"""LONI IDA search export: discovery, loading, classification and per-session lists.

Since the 2026-11-16 download each cohort has ONE ``idaSearch_<m>_<d>_<yyyy>_*.csv``
export with ``Subject ID, Sex, [Research Group], Visit, Study Date,
[Archive Date], Age, [Weight], Modality, Description, Type, Imaging Protocol,
Image ID``. ``Imaging Protocol`` holds ``key=value;key=value`` pairs
(Manufacturer, Mfg Model, Field Strength, Weighting, Acquisition Type,
Radiopharmaceutical, ...).

Main entry points:

* :func:`ida_find` - new-format exports in a cohort folder, plus legacy
  per-visit files in ``old_ida_search/``;
* :func:`ida_load` - normalised table (one row per image series);
* :func:`ida_classify` - canonical type per row, computed once per distinct
  description;
* :func:`ida_build_acquisitions` - one row per subject, one token per SESSION
  (subject x type x study date) in every list;
* :func:`ida_combine_blocks` - subject x type coalescing of two sources (never
  mixes the lists of one source with another's);
* :func:`ida_match_index` - casts IDA subject IDs to the dtype of the clinical
  IDs before joining (text vs integer IDs silently duplicate subjects).
"""
from __future__ import annotations

import numpy as np
import os
import pandas as pd
import re

from ..imaging.scanner import norm_vendor

import logging as _logging
from ..logging_utils import make_print as _make_print

log = _logging.getLogger(__name__)
_print = _make_print(log)

# =============================================================================
# LONI IDA SEARCH HELPERS  (single-file export — added 2026-11-16)
# =============================================================================
# Since the 2026-11-16 download each cohort has ONE  idaSearch_<m>_<d>_<y>*.csv
# in Tables_<COHORT>/ ; the old per-visit files moved to
# Tables_<COHORT>/old_ida_search/ .  Columns of the new export:
#   Subject ID · Sex · [Research Group] · Visit · Study Date · Age · Modality ·
#   Description · Type · Imaging Protocol · Image ID  (+ Weight / Archive Date)
# The old export only had Subject ID / Sex / Age / Description [/ Type] and the
# visit was implied by the FILE NAME (one file per visit).
#
# What the helpers do (used by 4RTNI, NIFD, ABVIB, MCSA, A4 and AIBL):
#   ida_find()               newest-format file(s) in the cohort folder — no
#                            hard-coded download date in the pipelines any more
#   ida_load()               one lean read (only the needed columns, text
#                            dtype), 'Imaging Protocol' parsed once per distinct
#                            string -> vendor / model / field strength / tracer
#   ida_classify()           Description -> canonical type, computed once per
#                            DISTINCT description (not once per row)
#   ida_build_acquisitions() one row per subject; one token per SESSION
#                            (subject x type x study date) in every list:
#       all_{T}_date            study date of each session (ISO)
#       all_{T}_visit           IDA visit of each session (cohort token)
#       {T}_IDs                 LONI Image ID of each session (lowest ID,
#                               'Original' before 'Pre-processed')
#       all_{T}_image_ids       every Image ID of the session, joined by '/'
#       all_{T}_manufacturer / all_{T}_model / all_{T}_field_strength
#       {T}_manufacturer        vendor(s) seen for the subject
#       nb_total_{T}            number of sessions   (= tokens in every list)
#       nb_images_{T}           number of image series
#       1st_/Lst_{T}_visit, 1st_/Lst_{T}_ID
#     + MRI_type_disponible / PET_tracer_disponible
#   ida_combine_blocks()     subject x type coalescing of two sources (never
#                            mixes the lists of one source with another's)
#   ida_match_index()        casts the IDA subject IDs to the dtype of the
#                            clinical IDs (int vs str) before the outer join
#
# LEGACY FILL: when the new export lacks rows that the old per-visit files
# have (e.g. a search downloaded without the DTI / fMRI modalities), those
# rows are appended from old_ida_search/ — the session date and scanner are
# taken from the same subject+visit of the new export, Image ID stays 'NA'.
# Subjects touched are flagged in the column IDA_legacy_types. With a complete
# export the fill adds 0 rows.
# =============================================================================
IDA_FILE_RX     = re.compile(r"^idaSearch_(\d{1,2})_(\d{1,2})_(\d{4}).*\.csv$", re.I)
IDA_LEGACY_DIRS = ("old_ida_search",)
IDA_COLMAP = {
    "subject id": "subject", "sex": "sex", "research group": "group",
    "visit": "visit", "study date": "study_date", "age": "age",
    "modality": "modality", "description": "description", "type": "type",
    "imaging protocol": "protocol", "image id": "image_id",
}
IDA_PROTOCOL_KEYS = {
    "Manufacturer": "vendor_raw", "Mfg Model": "model",
    "Field Strength": "field_strength", "Weighting": "weighting",
    "Acquisition Type": "acq_type", "Radiopharmaceutical": "tracer_raw",
}
# Radiopharmaceutical (DICOM field inside 'Imaging Protocol') -> canonical tracer
IDA_TRACER_RULES = [
    (r"pib|pittsburgh",                 "ABETA_PIB"),
    (r"av-?45|florbetapir|\bfbp\b",     "ABETA_AV45"),
    (r"flutemetamol|vizamyl",           "ABETA_Flutemetamol"),
    (r"nav",                            "ABETA_NAV4694"),
    (r"florbetaben|\bfbb\b|neuraceq",   "ABETA_FBB"),
    (r"av-?1451|flortaucipir|t807",     "TAU_Flortaucipir"),
    (r"mk-?6240",                       "TAU_MK6240"),
    (r"pi-?2620",                       "TAU_PI2620"),
    (r"fdg|fluorodeoxyglucose",         "PET_FDG"),
]
_IDA_EMPTY = ("", "nan", "none", "nat", "na")


def ida_tracer(raw):
    """'11C-PIB' / 'Pittsburgh compound B C^11^' -> ABETA_PIB ; unknown -> None."""
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return None
    s = str(raw).strip().lower()
    if s in _IDA_EMPTY:
        return None
    for rx, canon in IDA_TRACER_RULES:
        if re.search(rx, s):
            return canon
    return None


def _ida_file_key(path):
    m = IDA_FILE_RX.match(os.path.basename(path))
    mo, d, y = (int(x) for x in m.groups()) if m else (0, 0, 0)
    return (y, mo, d, os.path.getmtime(path))


def _ida_header(path):
    for enc in ("utf-8", "latin-1", "cp1252"):
        try:
            return [c.strip().lower() for c in
                    pd.read_csv(path, nrows=0, encoding=enc).columns]
        except UnicodeDecodeError:
            continue
        except Exception:
            return []
    return []


def ida_find(folder: str) -> dict:
    """
    {'new': [...], 'legacy': [...]}  — IDA search CSVs of one cohort folder.
      new    : files of the single-file format (have 'Study Date' + 'Image ID'),
               looked up in the folder root, newest first
      legacy : old exports (no Image ID) from the root and from old_ida_search/
    """
    out = {"new": [], "legacy": []}
    if not os.path.isdir(folder):
        return out
    places = [folder] + [os.path.join(folder, d) for d in IDA_LEGACY_DIRS]
    for i, place in enumerate(places):
        if not os.path.isdir(place):
            continue
        for fn in sorted(os.listdir(place)):
            if not IDA_FILE_RX.match(fn):
                continue
            p = os.path.join(place, fn)
            hdr = _ida_header(p)
            if "subject id" not in hdr:
                continue
            is_new = "image id" in hdr and "study date" in hdr
            out["new" if (is_new and i == 0) else "legacy"].append(p)
    out["new"].sort(key=_ida_file_key, reverse=True)
    out["legacy"].sort()
    return out


def _ida_read_raw(path: str) -> pd.DataFrame:
    """Lean read: only the columns in IDA_COLMAP, everything as text."""
    for enc in ("utf-8", "latin-1", "cp1252"):
        try:
            df = pd.read_csv(path, dtype=str, keep_default_na=False, encoding=enc,
                             usecols=lambda c: c.strip().lower() in IDA_COLMAP,
                             on_bad_lines="warn")
            break
        except UnicodeDecodeError:
            continue
        except Exception as e:
            _print(f"    [ERROR] {os.path.basename(path)}: {e}")
            return pd.DataFrame()
    else:
        _print(f"    [ERROR] {os.path.basename(path)}: no encoding worked")
        return pd.DataFrame()
    df.columns = [IDA_COLMAP[c.strip().lower()] for c in df.columns]
    for c in df.columns:
        df[c] = df[c].str.strip()
    return df


def _ida_parse_protocol(df: pd.DataFrame) -> pd.DataFrame:
    """'Imaging Protocol' ('k=v;k=v;…') -> one column per key in
    IDA_PROTOCOL_KEYS. Each DISTINCT protocol string is parsed once."""
    prot = df["protocol"] if "protocol" in df.columns else pd.Series("", index=df.index)
    parsed = {}
    for p in pd.unique(prot):
        kv = {}
        for item in str(p).split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                kv[k.strip()] = v.strip()
        parsed[p] = kv
    for key, col in IDA_PROTOCOL_KEYS.items():
        lut = {p: (kv.get(key) or np.nan) for p, kv in parsed.items()}
        df[col] = prot.map(lut)
    df["vendor"] = df["vendor_raw"].map(
        {v: norm_vendor(v) for v in pd.unique(df["vendor_raw"].dropna())})
    fs = pd.to_numeric(df["field_strength"], errors="coerce")
    df["field_strength"] = fs.map(lambda x: f"{x:g}T" if pd.notna(x) and x > 0 else np.nan)
    df["tracer"] = df["tracer_raw"].map(
        {v: ida_tracer(v) for v in pd.unique(df["tracer_raw"].dropna())})
    return df.drop(columns=[c for c in ("protocol", "vendor_raw") if c in df.columns])


def ida_load(folder: str, legacy_visit=None, legacy_fill: bool = True,
             max_session_span_days: int = 31, label: str = "") -> pd.DataFrame:
    """
    Normalised IDA table of one cohort (one row per image series):
      subject, sex, group, visit, date (datetime64), age, modality, description,
      type, image_id, vendor, model, field_strength, weighting, acq_type,
      tracer_raw, tracer (canonical), source ('ida' | 'legacy')
    legacy_visit : fn(file_name) -> Visit label written like the new export
                   ('Month 6', 'Visit 2', '3' …) for the old per-visit files.
    Sentinels are cleaned: Age 0 / > 120 -> NaN, Sex 'X' -> NaN.
    """
    tag = label or os.path.basename(folder.rstrip("/\\"))
    found = ida_find(folder)
    cols = list(dict.fromkeys(IDA_COLMAP.values()))
    new = pd.DataFrame()
    if found["new"]:
        parts = [_ida_read_raw(p) for p in found["new"]]          # newest first
        parts = [p for p in parts if not p.empty]
        if parts:
            new = pd.concat(parts, ignore_index=True) if len(parts) > 1 else parts[0]
            if len(parts) > 1 and "image_id" in new.columns:      # union of exports
                new = new[~(new["image_id"].ne("") & new["image_id"].duplicated())]
        for p in found["new"]:
            _print(f"  [OK]   ida[{tag}]: {os.path.basename(p)}")
    else:
        _print(f"  [WARN] ida[{tag}]: no single-file export (idaSearch_*.csv with "
              f"'Study Date' + 'Image ID') in {folder}")
    new = new.reindex(columns=cols, fill_value="") if not new.empty else pd.DataFrame(columns=cols)
    new["source"] = "ida"

    # ── legacy per-visit files ───────────────────────────────────────────────
    leg = pd.DataFrame(columns=cols)
    if found["legacy"] and (legacy_fill or new.empty):
        parts = []
        for p in found["legacy"]:
            d = _ida_read_raw(p)
            if d.empty:
                continue
            if "visit" not in d.columns:
                if legacy_visit is None and not new.empty:
                    continue          # visit unknown -> rows cannot be matched to the new export
                d["visit"] = (legacy_visit(os.path.basename(p)) or "") if legacy_visit else ""
            parts.append(d)
        if parts:
            leg = pd.concat(parts, ignore_index=True).reindex(columns=cols, fill_value="")
    if not leg.empty:
        key = ["subject", "visit", "description"]
        leg["_occ"] = leg.groupby(key).cumcount()
        if not new.empty:
            have = set(zip(new["subject"], new["visit"], new["description"],
                           new.groupby(key).cumcount()))
            keep = [k not in have for k in zip(leg["subject"], leg["visit"],
                                               leg["description"], leg["_occ"])]
            leg = leg[keep]
        leg = leg.drop(columns="_occ")
        leg["source"] = "legacy"

    df = pd.concat([new, leg], ignore_index=True) if not leg.empty else new
    if df.empty:
        return df
    df["date"] = pd.to_datetime(df["study_date"], format="%m/%d/%Y", errors="coerce")
    df = df.drop(columns="study_date")
    age = pd.to_numeric(df["age"], errors="coerce")
    df["age"] = age.where((age > 0) & (age <= 120))
    df["sex"] = df["sex"].where(df["sex"].isin(["M", "F"]))
    df = _ida_parse_protocol(df)

    # legacy rows: date / scanner of the same subject+visit in the new export
    is_leg = df["source"].eq("legacy")
    n_leg = int(is_leg.sum())
    if n_leg and not new.empty:
        ref = df[~is_leg & df["date"].notna() & df["modality"].ne("PET")]
        g = ref.groupby(["subject", "visit"])
        span = (g["date"].max() - g["date"].min()).dt.days
        one = lambda s: s.dropna().iloc[0] if s.dropna().nunique() == 1 else np.nan
        lut = pd.DataFrame({
            "date": g["date"].min().where(span <= max_session_span_days),
            "vendor": g["vendor"].agg(one), "model": g["model"].agg(one),
            "field_strength": g["field_strength"].agg(one),
        })
        idx = pd.MultiIndex.from_arrays([df.loc[is_leg, "subject"], df.loc[is_leg, "visit"]])
        for c in lut.columns:
            df.loc[is_leg, c] = lut[c].reindex(idx).values
        df["date"] = pd.to_datetime(df["date"])
        dated = int(df.loc[is_leg, "date"].notna().sum())
        _print(f"  [WARN] ida[{tag}]: {n_leg:,} rows exist only in the OLD per-visit "
              f"files -> appended (legacy fill; {dated:,} got the session date of "
              f"the same subject+visit, no Image ID). Re-download the IDA search "
              f"with every modality to replace them.")
    elif n_leg:
        _print(f"  [WARN] ida[{tag}]: using {n_leg:,} rows of the OLD per-visit files "
              f"only (visit axis, no date / Image ID / scanner).")

    d = df["date"].dropna()
    rng = f"{d.min().date()} -> {d.max().date()}" if len(d) else "no dates"
    _print(f"         {len(df):,} rows · {df['subject'].nunique():,} subjects · "
          f"study dates {rng} · modalities "
          f"{df.loc[df['modality'].ne(''), 'modality'].value_counts().to_dict()}")
    if len(d) > 50 and ((d.dt.month == 1) & (d.dt.day == 1)).mean() > 0.9:
        _print(f"  [INFO] ida[{tag}]: study dates are year-only (01/01/YYYY, "
              f"de-identified) -> not usable as a time axis.")
    return df


def ida_classify(ida: pd.DataFrame, fn, col: str = "description") -> pd.Series:
    """fn(description) -> canonical type | None, evaluated once per DISTINCT value."""
    if ida.empty:
        return pd.Series(dtype=object)
    lut = {d: fn(d) for d in pd.unique(ida[col])}
    return ida[col].map(lut).where(lambda s: s.notna(), np.nan)


def ida_report_unclassified(ida: pd.DataFrame, type_col: str = "canon",
                            skip_fn=None, n: int = 8) -> None:
    """Print the most frequent descriptions that got NO canonical type and are
    not infrastructure (skip_fn) — the place to look when a sequence is missing."""
    if ida.empty:
        return
    un = ida[ida[type_col].isna()]
    if skip_fn is not None and not un.empty:
        lut = {d: bool(skip_fn(d)) for d in pd.unique(un["description"])}
        un = un[~un["description"].map(lut)]
    if un.empty:
        return
    top = (un.groupby(["modality", "description"]).size()
             .sort_values(ascending=False).head(n))
    ex = "; ".join(f"{m or '?'}:'{d}'({c})" for (m, d), c in top.items())
    _print(f"    [INFO] {len(un):,} rows without canonical type (not counted): {ex}")


def _ida_tok(v) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "NA"
    s = str(v).strip()
    return s if s and s.lower() not in ("nan", "none", "nat") else "NA"


def _ida_cell(tokens):
    """'|'-joined cell, NaN when no token carries information."""
    return "|".join(tokens) if any(t != "NA" for t in tokens) else np.nan


def ida_build_acquisitions(ida: pd.DataFrame, type_col: str = "canon",
                           use_dates: bool = True,
                           tok_col: str = "visit_tok", ord_col: str = "visit_ord",
                           first_age: bool = False,
                           value_cols: tuple = ("suvr", "centiloid")) -> pd.DataFrame:
    """
    Per-subject acquisition table (index 'ID') from a classified IDA table.
    A SESSION is subject x type x study date (x visit token when the date is
    missing or use_dates=False). Every all_{T}_* list and {T}_IDs carry exactly
    nb_total_{T} tokens, in chronological order — see the header of this block.
    value_cols: numeric per-row columns (when present in `ida`, e.g. PET 'suvr' /
    'centiloid' merged from a quantification table) written as all_{T}_{col},
    one value per session (mean of the session's rows, 'NA' when missing).
    IDA_sex: sex recorded in the IDA export (M / F), one value per subject.
    """
    need = {"subject", type_col}
    if ida is None or ida.empty or not need.issubset(ida.columns):
        return pd.DataFrame()
    d = ida[ida[type_col].notna() & ida["subject"].astype(str).ne("")].copy()
    if d.empty:
        return pd.DataFrame()
    for c in ("image_id", "type", "vendor", "model", "field_strength", "source"):
        if c not in d.columns:
            d[c] = np.nan
    d["_dstr"] = (d["date"].dt.strftime("%Y-%m-%d")
                  if use_dates and "date" in d.columns else np.nan)
    d["_dstr"] = d["_dstr"].astype(object).where(d["_dstr"].notna(), np.nan)
    d["_tok"] = d[tok_col] if tok_col in d.columns else np.nan
    d["_ord"] = pd.to_numeric(d[ord_col], errors="coerce") if ord_col in d.columns else np.nan
    d["_skey"] = np.where(d["_dstr"].notna(), "D:" + d["_dstr"].astype(str),
                          "V:" + d["_tok"].map(_ida_tok))
    d["_idnum"] = pd.to_numeric(d["image_id"], errors="coerce")
    d["_pre"] = d["type"].astype(str).str.lower().str.startswith("pre").astype(int)
    d["_iid"] = d["image_id"].where(d["image_id"].astype(str).str.strip().ne("")
                                    & d["image_id"].notna())
    d["_pair"] = d["vendor"].map(_ida_tok) + "\t" + d["model"].map(_ida_tok)
    d = d.sort_values(["subject", type_col, "_dstr", "_ord", "_pre", "_idnum"],
                      na_position="last", kind="mergesort")

    def _ids(s):
        v = list(dict.fromkeys(x for x in s if isinstance(x, str)))
        return "/".join(v) if v else "NA"

    def _pairs(s):
        p = sorted(set(s))
        if len(p) > 1:                       # drop the all-NA pair if a real one exists
            p = [x for x in p if x != "NA\tNA"] or p
        return "/".join(x.split("\t")[0] for x in p) + "\n" + "/".join(x.split("\t")[1] for x in p)

    def _uniq(s):
        v = sorted({x for x in s if isinstance(x, str) and x})
        return "/".join(v) if v else "NA"

    d = d.rename(columns={type_col: "_T"})
    g = d.groupby(["subject", "_T", "_skey"], sort=False)
    agg = dict(dstr=("_dstr", "first"), tok=("_tok", "first"),
               ids=("_iid", _ids), n=("_skey", "size"), pair=("_pair", _pairs),
               fs=("field_strength", _uniq),
               leg=("source", lambda s: bool((s == "legacy").any())))
    if first_age and "age" in d.columns:
        agg["age"] = ("age", "min")
    value_cols = [c for c in (value_cols or ()) if c in d.columns
                  and pd.to_numeric(d[c], errors="coerce").notna().any()]
    for c in value_cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
        agg[f"v_{c}"] = (c, "mean")
    sess = g.agg(**agg).reset_index().rename(columns={"_T": "ctype"})
    if "age" not in sess.columns:
        sess["age"] = np.nan

    acc = {}                                     # (subject, T) -> lists, in order
    for r in sess.itertuples(index=False):
        a = acc.get((r.subject, r.ctype))
        if a is None:
            a = acc[(r.subject, r.ctype)] = {
                "date": [], "tok": [], "raw": [], "rep": [], "ids": [],
                "man": [], "mod": [], "fs": [], "n": 0, "leg": False, "age": r.age,
                **{f"v_{c}": [] for c in value_cols}}
        man, mod = r.pair.split("\n")
        a["date"].append(_ida_tok(r.dstr)); a["tok"].append(_ida_tok(r.tok))
        a["raw"].append(r.tok)
        a["ids"].append(r.ids); a["rep"].append(r.ids.split("/")[0])
        a["man"].append(man); a["mod"].append(mod); a["fs"].append(r.fs)
        a["n"] += int(r.n); a["leg"] = a["leg"] or bool(r.leg)
        for c in value_cols:
            v = getattr(r, f"v_{c}")
            a[f"v_{c}"].append("NA" if pd.isna(v) else f"{round(float(v), 4):g}")

    rows = {}
    for (subj, T), a in acc.items():
        row = rows.setdefault(subj, {"ID": subj, "_mri": set(), "_pet": set(), "_leg": set()})
        (row["_mri"] if str(T).startswith("MRI_") else row["_pet"]).add(T)
        if a["leg"]:
            row["_leg"].add(T)
        real_ids = [x for x in a["rep"] if x != "NA"]
        real_tok = [x for x, t in zip(a["raw"], a["tok"]) if t != "NA"]
        vend = sorted({x for v in a["man"] for x in v.split("/") if x != "NA"})
        if use_dates:
            row[f"all_{T}_date"] = _ida_cell(a["date"])
        row[f"all_{T}_visit"]          = _ida_cell(a["tok"])
        row[f"nb_total_{T}"]           = len(a["tok"])
        row[f"nb_images_{T}"]          = a["n"]
        row[f"{T}_IDs"]                = _ida_cell(a["rep"])
        row[f"all_{T}_image_ids"]      = _ida_cell(a["ids"])
        row[f"1st_{T}_ID"]             = real_ids[0]  if real_ids else np.nan
        row[f"Lst_{T}_ID"]             = real_ids[-1] if real_ids else np.nan
        row[f"1st_{T}_visit"]          = real_tok[0]  if real_tok else np.nan
        row[f"Lst_{T}_visit"]          = real_tok[-1] if real_tok else np.nan
        row[f"all_{T}_manufacturer"]   = _ida_cell(a["man"])
        row[f"all_{T}_model"]          = _ida_cell(a["mod"])
        row[f"all_{T}_field_strength"] = _ida_cell(a["fs"])
        row[f"{T}_manufacturer"]       = "|".join(vend) if vend else np.nan
        for c in value_cols:
            row[f"all_{T}_{c}"] = _ida_cell(a[f"v_{c}"])
        if first_age:
            row[f"1st_{T}_age"] = a["age"]
    for row in rows.values():
        row["MRI_type_disponible"]   = "|".join(sorted(row.pop("_mri"))) or np.nan
        row["PET_tracer_disponible"] = "|".join(sorted(row.pop("_pet"))) or np.nan
        row["IDA_legacy_types"]      = "|".join(sorted(row.pop("_leg"))) or np.nan
    out = pd.DataFrame(list(rows.values())).set_index("ID")
    if "sex" in ida.columns:                     # sex as recorded in the IDA export
        sx = ida.loc[ida["sex"].isin(["M", "F"])].drop_duplicates("subject")
        out["IDA_sex"] = sx.set_index("subject")["sex"].reindex(out.index)
    lead = ["MRI_type_disponible", "PET_tracer_disponible", "IDA_legacy_types"]
    out = out[lead + [c for c in out.columns if c not in lead]]
    return out.dropna(axis=1, how="all")         # e.g. no scanner info in the export


_IDA_BLOCK_LISTS = ("date", "visit", "image_ids", "manufacturer", "model",
                    "field_strength", "suvr", "centiloid")
_IDA_ID_COLS = ("{T}_IDs", "all_{T}_image_ids", "1st_{T}_ID", "Lst_{T}_ID", "nb_images_{T}")


def ida_block_columns(columns, T: str) -> list:
    """Columns forming the per-acquisition block of type T (exact names, so
    MRI_T2w never captures MRI_T2w_FLAIR)."""
    names = {f"{T}_IDs", f"nb_total_{T}", f"nb_images_{T}", f"{T}_manufacturer",
             f"1st_{T}_ID", f"Lst_{T}_ID", f"1st_{T}_visit", f"Lst_{T}_visit"} \
        | {f"all_{T}_{s}" for s in _IDA_BLOCK_LISTS}
    return [c for c in columns if c in names]


def ida_combine_blocks(primary: pd.DataFrame, secondary: pd.DataFrame,
                       borrow_ids: bool = False) -> pd.DataFrame:
    """
    Coalesce two per-subject imaging tables at the SUBJECT x TYPE level:
    where `primary` has type T for a subject (nb_total_T filled) its whole
    block is kept, otherwise the block of `secondary`. Lists of the two sources
    are never mixed inside a block (cell-wise combine_first did that).
    Other columns: primary first, then secondary. *_disponible = union.
    borrow_ids: when both have T with the SAME number of sessions, the Image-ID
    columns missing in `primary` are taken from `secondary`.
    """
    ok = lambda x: isinstance(x, pd.DataFrame) and not x.empty
    if not ok(primary):
        return secondary if ok(secondary) else pd.DataFrame()
    if not ok(secondary):
        return primary
    idx = primary.index.union(secondary.index)
    P, S = primary.reindex(idx), secondary.reindex(idx)
    out = P.combine_first(S)
    out = out[list(dict.fromkeys(list(primary.columns) + list(secondary.columns)))]
    nan = pd.Series(np.nan, index=idx, dtype=object)
    types = sorted({c[len("nb_total_"):] for c in out.columns if c.startswith("nb_total_")})
    for T in types:
        nb = f"nb_total_{T}"
        has_p = P[nb].notna() if nb in P.columns else pd.Series(False, index=idx)
        has_s = S[nb].notna() if nb in S.columns else pd.Series(False, index=idx)
        same = has_p & has_s & (pd.to_numeric(P[nb], errors="coerce")
                                == pd.to_numeric(S[nb], errors="coerce")) \
            if (nb in P.columns and nb in S.columns) else pd.Series(False, index=idx)
        id_cols = {t.format(T=T) for t in _IDA_ID_COLS}
        for c in ida_block_columns(out.columns, T):
            p = P[c] if c in P.columns else nan
            s = S[c] if c in S.columns else nan
            col = p.where(has_p, s)
            if borrow_ids and c in id_cols:
                col = col.where(~(same & p.isna()), s)
            out[c] = col
    for c in ("MRI_type_disponible", "PET_tracer_disponible", "IDA_legacy_types"):
        if c in P.columns and c in S.columns:
            out[c] = [("|".join(sorted({t for v in (a, b) if isinstance(v, str)
                                        for t in v.split("|") if t})) or np.nan)
                      for a, b in zip(P[c], S[c])]
    out.index.name = "ID"
    return out


def ida_match_index(result: pd.DataFrame, ref_ids) -> pd.DataFrame:
    """Give the IDA table the ID dtype of the clinical tables. IDA subject IDs
    are text; a clinical table with integer IDs would otherwise not join and
    every subject would appear twice (once clinical, once imaging)."""
    if result is None or result.empty:
        return result
    ref = pd.Series(list(ref_ids)).dropna()
    if len(ref) and pd.api.types.is_numeric_dtype(ref):
        num = pd.to_numeric(pd.Series(result.index.astype(str)), errors="coerce")
        if num.notna().all() and (num % 1 == 0).all():
            result = result.copy()
            result.index = pd.Index(num.astype("int64").values, name="ID")
            return result
        _print(f"    [WARN] {int(num.isna().sum())} IDA subject IDs are not numeric "
              f"-> kept as text (they will not join integer clinical IDs)")
    result = result.copy()
    result.index = pd.Index(result.index.astype(str).str.strip(), name="ID")
    return result
