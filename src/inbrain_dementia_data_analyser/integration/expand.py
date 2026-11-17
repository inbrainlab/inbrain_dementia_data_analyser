"""Expanded master: every measure placed into time windows around the first diagnosis.
"""
from __future__ import annotations

import ast
import gc
import numpy as np
import os
import pandas as pd
import re
from datetime import date as _date

from ..config import ProjectConfig, load_config
from ..provenance import write_provenance

from ..config import ProjectConfig, load_config

import logging as _logging
from ..logging_utils import make_print as _make_print
from ..logging_utils import quiet_pandas_warnings

log = _logging.getLogger(__name__)
_print = _make_print(log)

@quiet_pandas_warnings
def expand_master(master_path: str | os.PathLike, out_path: str | os.PathLike, *,
                  config: ProjectConfig | None = None) -> str:
    """Place every measure of every subject into time windows.

    The anchor of a subject is its first-diagnosis date (or, when it has none,
    the measure's own first acquisition). Each acquisition goes to a 6-month
    window ``lt6 ... lt240`` / ``gt240`` after the anchor or ``pre6 ...`` /
    ``pregt240`` before it. Per measure and window the output holds
    ``<stem>_<w>`` (value; 1 for imaging presence), ``<stem>_<w>_acq`` (date or
    day offset), ``<stem>_<w>_id`` and the extra channels (image IDs, scanner,
    SUVR, centiloid; batch / run date for fluids). Several acquisitions in one
    window are stored as a Python-list string.

    The measures are the variables flagged ``expand: true`` in
    ``variables.yaml``, the imaging types of ``imaging.yaml`` and every family
    whose name starts with one of ``integration.yaml -> expand.auto_feature_prefixes``.

    Memory: only the needed columns are read (in row chunks), window columns
    are created only when they receive a value, and the output is written in
    blocks of rows.

    Parameters
    ----------
    master_path:
        ``MASTER_all_cohorts.csv`` written by :func:`~.merge.merge_cohorts`.
    out_path:
        Output CSV (``MASTER_expanded_240m.csv``).
    config:
        Project configuration.

    Returns
    -------
    str
        ``out_path``.
    """
    cfg = config or load_config()
    E = cfg.integration["expand"]
    MASTER_PATH, OUT_PATH = str(master_path), str(out_path)
    if not os.path.exists(MASTER_PATH):
        raise FileNotFoundError(f"Master table not found: {MASTER_PATH} (run merge_cohorts first)")
    FEATURES = cfg.expand_features()

    # Extra quantitative channels captured per window for imaging features
    # (only emitted when all_{feat}_{channel} exists, e.g. amyloid PET SUVR/Centiloid).
    # image_ids = every Image ID of the session ('/'-joined); {feat}_{w}_id keeps the
    # single representative ID that is aligned with _acq.
    IMAGING_VALUE_CHANNELS = list(E["imaging_value_channels"])

    # Every fluid biomarker / lab family present in the master is added
    # automatically as a 'scored' feature: all_{Plasma_*|CSF_*|Lab_*}_score
    # (ADNI plasma panels, extra ADNI CSF tables, A4 plasma, PREVENT-AD 4plex,
    # ADNI LABDATA, AIBL lab …). Set to () to keep only the FEATURES above.
    AUTO_FEATURE_PREFIXES = tuple(E["auto_feature_prefixes"])

    # Per-measurement metadata expanded per window next to the value, aligned
    # with _acq (only when all_{feat}_{channel} exists in the master):
    #   method/batch = assay batch / run, rundate = date the sample was assayed
    FLUID_META_CHANNELS = list(E["fluid_meta_channels"])

    # Per-subject metadata carried unchanged for every expanded feature
    # (laboratory, assay kit, number of censored <LoQ/>ULoQ values).
    FEATURE_META_SUFFIXES = list(E["feature_meta_suffixes"])

    # Features whose acquisition DATE lives under another feature (same draw/scan).
    # EMPTY on purpose: borrowing another marker's date list pairs dates and values
    # by position, which is wrong whenever the two markers were not measured at
    # exactly the same visits (e.g. ADNI Aβ40 exists only for part of the Aβ42
    # punctures). Every pipeline now writes its own aligned all_{feat}_date.
    DATE_ALIAS = dict(E.get("date_alias") or {})

    # PER-ROW date fallback for masters built with the OLD ADNI pipeline
    # (all_{m}_score present, all_{m}_date absent; date only in the shared
    # all_pct_date_collection). Applied ONLY when the subject has exactly as many
    # values as shared collection dates — values are a subset of the punctures,
    # so equal length guarantees the k-th value belongs to the k-th date.
    # Otherwise the subject's values are left untimed (counted and reported).
    SHARED_DATE_FALLBACK_COL = E["shared_date_fallback_column"]
    SHARED_DATE_FALLBACK_FEATURES = set(E["shared_date_fallback_features"])

    # Tokens meaning "no date" inside an aligned all_*_date list
    EMPTY_TOKENS = tuple(E["empty_tokens"])

    # If a subject has NO first-diagnosis date, anchor that feature to ITS OWN first
    # acquisition instead of dropping the subject. Keeps preclinical / CN-only
    # cohorts (e.g. A4, plasma-only subjects) from vanishing. The diagnosis spread
    # and conversion timing always stay anchored on the diagnosis.
    ANCHOR_FALLBACK_TO_FIRST_ACQ = bool(E["anchor_fallback_to_first_acquisition"])

    # identifiers / diagnosis / demographics carried into the expanded file
    CARRY_COLUMNS = list(E["carry_columns"])

    WINDOW_STEP = int(E["window_step_months"])
    WINDOW_MAX = int(E["window_max_months"])

    DAY_OFFSET_BANKS = set(E["day_offset_banks"])
    AGE_MONTH_BANKS = set(E["age_month_banks"])
    VISIT_NUM_BANKS = set(E["visit_number_banks"])
    DAYS_PER_MONTH = float(E["days_per_month"])

    # =============================================================================
    # WINDOW GRID
    # =============================================================================
    WIN_BOUNDS = list(range(WINDOW_STEP, WINDOW_MAX + 1, WINDOW_STEP))


    def assign_window(months):
        if months is None:
            return None
        if months < 0:                       # before first diagnosis -> binned pre-windows
            am = -months
            for m in WIN_BOUNDS:
                if am < m:
                    return f"pre{m}"         # pre6 = 0-6 mo before, pre12 = 6-12 mo before, ...
            return f"pregt{WINDOW_MAX}"
        for m in WIN_BOUNDS:
            if months < m:
                return f"lt{m}"
        return f"gt{WINDOW_MAX}"


    # =============================================================================
    # PARSING + TIME
    # =============================================================================
    def parse_list(val):
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return []
        s = str(val).strip()
        if s == "" or s.lower() in ("nan", "nat", "none"):
            return []
        if s.startswith("[") and s.endswith("]"):
            s2 = re.sub(r"\bnan\b", "None", s)
            s2 = re.sub(r"\bNaT\b", "None", s2)
            try:
                out = ast.literal_eval(s2)
                return [("" if x is None else str(x)) for x in out] if isinstance(out, list) else []
            except Exception:
                pass
        return [tok.strip() for tok in s.split("|")]


    def _is_num(x):
        try:
            float(str(x).replace(" ", ""))
            return ("-" not in str(x)) and ("/" not in str(x))
        except Exception:
            return False



    def parse_date_fast(s):
        """Fast ISO date parse ('YYYY-MM-DD'); falls back to pandas only if needed."""
        if s is None:
            return None
        t = str(s).strip()
        if t in ("", "nan", "NaT", "None"):
            return None
        t = t[:10]
        try:
            return _date.fromisoformat(t)
        except Exception:
            dt = pd.to_datetime(t, errors="coerce")
            return None if pd.isna(dt) else dt.date()


    def prep_anchor(anchor, bank):
        """Parse a baseline anchor ONCE. Type follows the VALUE, not the cohort:
        numeric -> ('num', x) (offset/age/visit, interpreted later via bank);
        ISO date -> ('date', date)."""
        s = str(anchor).strip()
        if s in ("", "nan", "NaT", "None"):
            return ("date", None)
        if _is_num(s):
            return ("num", float(s))
        return ("date", parse_date_fast(s))


    def months_since(acq_val, anchor_kind, anchor_pv, bank):
        """Months from the (pre-parsed) first-diagnosis anchor to an acquisition."""
        if acq_val is None or str(acq_val).strip() in ("", "nan", "NaT", "None"):
            return None
        s = str(acq_val).strip()
        if anchor_kind == "num":
            if anchor_pv is None or not _is_num(s):
                return None
            x = float(s)
            if bank in DAY_OFFSET_BANKS:
                return (x - anchor_pv) / DAYS_PER_MONTH
            if bank in AGE_MONTH_BANKS:
                return x - anchor_pv
            return (x - anchor_pv) * WINDOW_STEP            # visit-number axis
        # calendar
        if anchor_pv is None:
            return None
        d = parse_date_fast(s)
        if d is None:
            return (float(s) - 1.0) * WINDOW_STEP if _is_num(s) else None
        return (d - anchor_pv).days / DAYS_PER_MONTH


    def first_token(val):
        toks = parse_list(val)
        return toks[0] if toks else ""


    def first_or_list(items):
        items = [x for x in items if str(x).strip() not in ("", "nan", "None")]
        if not items:
            return np.nan
        return items[0] if len(items) == 1 else str(items)

    def aligned_or_first(items):
        """Like first_or_list, but KEEPS positions (empty -> "") so that the
        list of IDs / channels stays aligned index by index with the _acq list."""
        if not items:
            return np.nan
        cleaned = ["" if str(x).strip() in ("", "nan", "None") else str(x) for x in items]
        if len(cleaned) == 1:
            return np.nan if cleaned[0] == "" else cleaned[0]
        return str(cleaned)

    # =============================================================================
    # DECIDE WHICH INPUT COLUMNS TO READ  (avoid loading all 1692)
    # =============================================================================
    _print("Inspecting master header ...")
    header = pd.read_csv(MASTER_PATH, nrows=0).columns.tolist()
    hset = set(header)

    # auto-add fluid biomarker / lab families
    if AUTO_FEATURE_PREFIXES:
        _auto = sorted({m.group(1) for c in header
                        if (m := re.match(r"all_(.+)_score$", c))
                        and m.group(1).startswith(AUTO_FEATURE_PREFIXES)})
        _new = [f for f in _auto if f not in FEATURES]
        for f in _new:
            FEATURES[f] = "scored"
        _print(f"  auto-added {len(_new)} fluid/lab features: {_new[:12]}"
              + (" …" if len(_new) > 12 else ""))

    feature_raw = {}
    for feat in FEATURES:
        cols = [c for c in (f"all_{feat}_date", f"all_{feat}_visit", f"all_{feat}_score",
                            f"{feat}_IDs", f"all_{feat}_ID", f"nb_total_{feat}") if c in hset]
        if feat in DATE_ALIAS and f"all_{DATE_ALIAS[feat]}_date" in hset:
            cols.append(f"all_{DATE_ALIAS[feat]}_date")          # borrowed date source
        if FEATURES[feat] == "imaging":
            cols += [f"all_{feat}_{ch}" for ch in IMAGING_VALUE_CHANNELS
                     if f"all_{feat}_{ch}" in hset]
        else:
            cols += [f"all_{feat}_{ch}" for ch in FLUID_META_CHANNELS
                     if f"all_{feat}_{ch}" in hset]
        if cols:
            feature_raw[feat] = cols

    carry_present = [c for c in CARRY_COLUMNS if c in hset]
    # lab / kit / censored-count columns of every feature (and of the CSF block)
    meta_carry = [c for feat in FEATURES
                  for c in [f"{feat}_{sfx}" for sfx in FEATURE_META_SUFFIXES]
                           + [f"nb_of_{feat}_censored"]
                  if c in hset]
    carry_present = list(dict.fromkeys(carry_present + meta_carry))

    # --- CSF detection: shared collection date + per-marker *_values ----------
    CSF_DATE_COL = "all_pct_date_collection"
    csf_markers = sorted({m.group(1) for c in header
                          if (m := re.match(r"all_(.+)_values$", c))})
    csf_present = (CSF_DATE_COL in hset) and len(csf_markers) > 0
    csf_cols = []
    if csf_present:
        csf_cols = ([CSF_DATE_COL]
                    + [f"all_{m}_values" for m in csf_markers]
                    + [f"all_{m}_date" for m in csf_markers if f"all_{m}_date" in hset]
                    + [f"all_{m}_method" for m in csf_markers if f"all_{m}_method" in hset])
        meta_carry += [c for m in csf_markers
                       for c in (f"{m}_lab", f"{m}_kit", f"nb_of_{m}_censored") if c in hset]
        carry_present = list(dict.fromkeys(carry_present + meta_carry))

    needed = list(dict.fromkeys(
        carry_present
        + [c for cols in feature_raw.values() for c in cols]
        + csf_cols
        + ([SHARED_DATE_FALLBACK_COL] if SHARED_DATE_FALLBACK_COL in hset else [])))
    _print(f"  reading {len(needed)} of {len(header)} columns")

    # Read in row chunks, everything as TEXT (2026-11-16):
    #  * pd.read_csv(low_memory=False) tokenises the WHOLE file (all ~2000 columns)
    #    before applying usecols; in chunks the RAM peak follows only the columns
    #    actually used — this is what lets the NACC-sized master fit in Colab;
    #  * every list / value here is parsed from text anyway, and text keeps the
    #    image IDs intact (as numbers, 1321166 became '1321166.0').
    df = pd.concat(pd.read_csv(MASTER_PATH, usecols=needed, dtype=str, chunksize=4000),
                   ignore_index=True)
    n = len(df)
    COHORT_COL = "Databank" if "Databank" in df.columns else \
                 ("dataset" if "dataset" in df.columns else None)
    banks = df[COHORT_COL].astype(str).tolist() if COHORT_COL else [""] * n
    _print(f"  loaded {n} subjects x {df.shape[1]} columns | cohort col = '{COHORT_COL}'")

    # =============================================================================
    # SUBJECT-LEVEL SEIZURE SUMMARY  (added 2026-11-08)
    # =============================================================================
    # Seizure history is mostly recorded ONCE (NACC/OASIS UDS3 A5 only at the
    # initial visit; ADNI/ABVIB/PREVENT-AD only free text), so it rarely falls
    # inside a ±6-month window of an MRI. These per-subject columns are carried
    # unchanged into the expanded master so the cohort creator can keep/filter them.
    #   SEIZURES_ever          max coded SEIZURES (0 absent, 1 recent/active, 2 remote)
    #   EPILEP_ever            1 if any UDS D1 epilepsy diagnosis
    #   Antiepileptic_med_ever 1 if NACCEPMD reported at any visit
    #   SEIZURE_any            1 = any positive source, 0 = has >= 1 source and all
    #                          negative, NaN = no seizure information at all
    #   SEIZURE_source         which sources were positive (or 'none_positive')
    def _max_list(v):
        xs = []
        for t in parse_list(v):
            try:
                xs.append(float(t))
            except Exception:
                pass
        return max(xs) if xs else np.nan

    SEIZURE_SOURCES = {                      # column -> (summary name, positive test)
        "all_SEIZURES_score":          ("SEIZURES_ever",          lambda x: x in (1, 2)),
        "all_EPILEP_score":            ("EPILEP_ever",            lambda x: x == 1),
        "all_Antiepileptic_med_score": ("Antiepileptic_med_ever", lambda x: x == 1),
    }
    _seiz_pos = pd.DataFrame(index=df.index)
    for col, (nm, pos) in SEIZURE_SOURCES.items():
        if col in df.columns:
            df[nm] = df[col].map(_max_list)
            _seiz_pos[nm] = df[nm].map(lambda x: np.nan if pd.isna(x) else float(pos(x)))
    for col in ("SEIZURE_text_flag", "AED_text_flag"):
        if col in df.columns:
            _seiz_pos[col] = pd.to_numeric(df[col], errors="coerce")
    if _seiz_pos.shape[1]:
        df["SEIZURE_any"] = _seiz_pos.max(axis=1, skipna=True)
        df["SEIZURE_source"] = [
            ("|".join(c for c in _seiz_pos.columns if r[c] == 1) or "none_positive")
            if r.notna().any() else np.nan
            for _, r in _seiz_pos.iterrows()]
        _new = [c for c in ["SEIZURES_ever", "EPILEP_ever", "Antiepileptic_med_ever",
                            "SEIZURE_any", "SEIZURE_source"] if c in df.columns]
        carry_present = list(dict.fromkeys(carry_present + _new))
        if COHORT_COL:
            _print("  seizure summary (subjects with info / positive) by cohort:")
            _print(pd.DataFrame({
                "with_info": df["SEIZURE_any"].notna().groupby(df[COHORT_COL]).sum(),
                "positive":  df["SEIZURE_any"].eq(1).groupby(df[COHORT_COL]).sum(),
            }).astype(int).to_string())

    # anchor = first-diagnosis date/visit
    if "1st_diag_date" in df.columns:
        base_all = df["all_diag_date"] if "all_diag_date" in df.columns else pd.Series([""] * n)
        anchors = [s if str(s).strip() not in ("", "nan", "NaT", "None") else first_token(a)
                   for s, a in zip(df["1st_diag_date"], base_all)]
    elif "all_diag_date" in df.columns:
        anchors = [first_token(a) for a in df["all_diag_date"]]
    else:
        raise KeyError("Need '1st_diag_date' or 'all_diag_date' to anchor baseline.")

    # parse each subject's anchor ONCE (big speedup vs re-parsing per visit/feature)
    anchor_kind = [None] * n
    anchor_pv   = [None] * n
    for i in range(n):
        k, v = prep_anchor(anchors[i], banks[i])
        anchor_kind[i] = k
        anchor_pv[i] = v

    # flag how each subject's baseline is anchored (for stratifying later)
    baseline_anchor = ["diagnosis" if anchor_pv[i] is not None else "first_acq"
                       for i in range(n)]
    n_no_diag = baseline_anchor.count("first_acq")
    _print(f"  baseline anchor: {n - n_no_diag} on diagnosis, {n_no_diag} on first acquisition")

    # =============================================================================
    # EXPAND ONE FEATURE  (lazy: only allocate windows that receive data)
    # =============================================================================
    def expand_feature(feat, ftype, raw_cols):
        is_img = (ftype == "imaging")
        g = lambda name: df[name].tolist() if name in raw_cols and name in df.columns else [""] * n
        # date source: own date, else borrowed via DATE_ALIAS
        date_src = f"all_{feat}_date"
        if date_src not in raw_cols and feat in DATE_ALIAS:
            date_src = f"all_{DATE_ALIAS[feat]}_date"
        L_date  = df[date_src].tolist() if date_src in df.columns else [""] * n
        L_visit = g(f"all_{feat}_visit")
        L_score = g(f"all_{feat}_score")
        L_ids   = (df[f"{feat}_IDs"].tolist() if f"{feat}_IDs" in raw_cols else
                   df[f"all_{feat}_ID"].tolist() if f"all_{feat}_ID" in raw_cols else [""] * n)
        L_nbtot = df[f"nb_total_{feat}"].tolist() if f"nb_total_{feat}" in raw_cols else [""] * n
        # per-row shared-date fallback (e.g. ADNI UPENN: score present, own date empty)
        use_shared = (feat in SHARED_DATE_FALLBACK_FEATURES
                      and SHARED_DATE_FALLBACK_COL in df.columns)
        L_shared = df[SHARED_DATE_FALLBACK_COL].tolist() if use_shared else [""] * n
        channels = ([ch for ch in IMAGING_VALUE_CHANNELS if f"all_{feat}_{ch}" in df.columns]
                    if is_img else
                    [ch for ch in FLUID_META_CHANNELS if f"all_{feat}_{ch}" in df.columns])
        L_chan = {ch: df[f"all_{feat}_{ch}"].tolist() for ch in channels}

        data = {}                                       # colname -> list (lazy)
        n_fallback = 0
        n_shared_ok = n_shared_bad = 0

        def cell(col, i, value):
            arr = data.get(col)
            if arr is None:
                arr = [np.nan] * n
                data[col] = arr
            arr[i] = value

        for i in range(n):
            acq  = parse_list(L_date[i]) or parse_list(L_visit[i])      # unified acq (#4)
            vals = parse_list(L_score[i])
            if not acq and use_shared and vals:                         # old ADNI export
                shared = parse_list(L_shared[i])
                if shared and len(shared) == len(vals):
                    acq = shared
                    n_shared_ok += 1
                else:
                    n_shared_bad += 1                                   # cannot align -> untimed
            ids  = parse_list(L_ids[i])
            chans = {ch: parse_list(L_chan[ch][i]) for ch in channels}
            ak, av, bank = anchor_kind[i], anchor_pv[i], banks[i]

            # anchor fallback: no diagnosis date -> use this feature's first acquisition
            if av is None and ANCHOR_FALLBACK_TO_FIRST_ACQ and acq:
                base = next((a for a in acq if str(a).strip() not in EMPTY_TOKENS), None)
                if base is not None:
                    ak, av = prep_anchor(base, bank)
                    if av is not None:
                        n_fallback += 1

            bucket = {}
            untimed_ids = []
            m = max([len(acq), len(vals), len(ids)] + [len(c) for c in chans.values()])
            for j in range(m):
                a   = acq[j]  if j < len(acq)  else ""
                v   = vals[j] if j < len(vals) else ""
                sid = ids[j]  if j < len(ids)  else ""
                if str(a).strip() in EMPTY_TOKENS:
                    if sid and str(sid).strip() not in EMPTY_TOKENS:
                        untimed_ids.append(sid)
                    continue
                w = assign_window(months_since(a, ak, av, bank))
                if w is None:
                    continue
                b = bucket.setdefault(w, {"val": [], "acq": [], "id": [],
                                          **{ch: [] for ch in channels}})
                b["val"].append("1" if (is_img and str(v).strip() in ("", "nan", "None")) else v)
                b["acq"].append(a)
                # append the ID IN LOCKSTEP with acq (placeholder "" when absent)
                # ('NA' = session without Image ID in the IDA-based pipelines)
                b["id"].append(sid if str(sid).strip() not in EMPTY_TOKENS else "")
                for ch in channels:
                    cv = chans[ch][j] if j < len(chans[ch]) else ""
                    b[ch].append(cv if str(cv).strip() not in ("", "nan", "None") else "")

            for w, d in bucket.items():
                cell(f"{feat}_{w}",     i, first_or_list(d["val"]))
                cell(f"{feat}_{w}_acq", i, first_or_list(d["acq"]))
                # ID / channels keeping positions ("" keeps them aligned with _acq)
                if is_img and any(str(x).strip() not in EMPTY_TOKENS for x in d["id"]):
                    cell(f"{feat}_{w}_id", i, aligned_or_first(d["id"]))
                for ch in channels:
                    if any(str(x).strip() not in EMPTY_TOKENS for x in d[ch]):
                        cell(f"{feat}_{w}_{ch}", i, aligned_or_first(d[ch]))

            if is_img:
                nb = 0
                if str(L_nbtot[i]).strip() not in ("", "nan", "None"):
                    try:
                        nb = int(float(L_nbtot[i]))
                    except Exception:
                        nb = 0
                timed = sum(len(d["acq"]) for d in bucket.values())
                extra = max(len(untimed_ids), nb - timed, 0)
                if extra > 0:
                    cell(f"{feat}_untimed_n", i, extra)
                    if untimed_ids:
                        cell(f"{feat}_untimed_id", i, first_or_list(untimed_ids))

        if use_shared and (n_shared_ok or n_shared_bad):
            _print(f"    [{feat}] shared CSF date used for {n_shared_ok} subj; "
                  f"{n_shared_bad} subj left untimed (values ≠ collection dates -> "
                  f"re-run the ADNI pipeline to get all_{feat}_date)")
        return (pd.DataFrame(data, index=df.index) if data else None), n_fallback


    # =============================================================================
    # RUN  (one feature at a time; free raw columns immediately)
    # =============================================================================
    def diagnose_zero(feat, raw_cols):
        """Explain why a feature produced 0 columns."""
        sc = f"all_{feat}_score"
        has_score = sc in df.columns
        date_src = f"all_{feat}_date"
        if date_src not in df.columns and feat in DATE_ALIAS:
            date_src = f"all_{DATE_ALIAS[feat]}_date"
        has_date = date_src in df.columns
        has_visit = f"all_{feat}_visit" in df.columns
        score_nn = int(df[sc].notna().sum()) if has_score else 0
        date_nn  = int(df[date_src].notna().sum()) if has_date else 0
        visit_nn = int(df[f"all_{feat}_visit"].notna().sum()) if has_visit else 0
        if not has_score and not has_date and not has_visit:
            return "no source columns in master"
        if score_nn == 0 and date_nn == 0 and visit_nn == 0:
            return "columns exist but are entirely empty"
        if has_score and not has_date and not has_visit:
            return "has score but NO date/visit column -> add a DATE_ALIAS"
        if score_nn > 0 and date_nn == 0 and visit_nn == 0:
            return f"{score_nn} scores but date/visit empty -> add DATE_ALIAS"
        return "values present but none could be anchored (no 1st_diag_date and fallback failed)"


    frames = []
    expanded, skipped = [], []
    total_fallback = {}
    alias_src_cols = {f"all_{src}_date" for src in DATE_ALIAS.values()}   # keep until end
    for feat, ftype in FEATURES.items():
        raw_cols = feature_raw.get(feat)
        if not raw_cols:
            skipped.append(feat)
            _print(f"  skipped  {feat:<22} (no source columns)")
            continue
        fr, nfb = expand_feature(feat, ftype, raw_cols)
        if fr is not None and fr.shape[1] > 0:
            frames.append(fr)
            expanded.append(feat)
            extra = f" | {nfb} subj used first-acq fallback" if nfb else ""
            _print(f"  expanded {feat:<22} (+{fr.shape[1]} cols){extra}")
            if nfb:
                total_fallback[feat] = nfb
        else:
            reason = diagnose_zero(feat, raw_cols)
            skipped.append(feat)
            _print(f"  expanded {feat:<22} (+0 cols)  -> {reason}")
        # free RAM, but keep alias date sources needed by later features
        df.drop(columns=[c for c in raw_cols
                         if c in df.columns and c not in alias_src_cols],
                inplace=True)
        gc.collect()

    # alias date sources can be dropped now (not part of CARRY)
    df.drop(columns=[c for c in alias_src_cols if c in df.columns], inplace=True)

    # =============================================================================
    # CSF — shared collection date + per-marker *_values  (windowed like the rest)
    # =============================================================================
    def expand_csf(markers, date_col):
        """
        *_values markers. Date source per marker and subject:
          1. all_{m}_date  (aligned with all_{m}_values — current pipelines)
          2. shared all_pct_date_collection, ONLY if it has as many entries as the
             values (old exports); otherwise the subject is left untimed.
        all_{m}_method is aligned with the values and written per window.
        """
        L_date = df[date_col].tolist() if date_col in df.columns else [""] * n
        marker_vals = {m: df[f"all_{m}_values"].tolist() for m in markers
                       if f"all_{m}_values" in df.columns}
        marker_date = {m: df[f"all_{m}_date"].tolist() for m in markers
                       if f"all_{m}_date" in df.columns}
        marker_meth = {m: df[f"all_{m}_method"].tolist() for m in markers
                       if f"all_{m}_method" in df.columns}
        data = {}
        skipped = {}

        def cell(col, i, value):
            arr = data.get(col)
            if arr is None:
                arr = [np.nan] * n
                data[col] = arr
            arr[i] = value

        for i in range(n):
            shared = parse_list(L_date[i])
            bank = banks[i]
            for m, vlist in marker_vals.items():
                vals = parse_list(vlist[i])
                if not vals:
                    continue
                dates = parse_list(marker_date[m][i]) if m in marker_date else []
                if not dates:
                    if shared and len(shared) == len(vals):
                        dates = shared
                    else:
                        skipped[m] = skipped.get(m, 0) + 1
                        continue
                meth = parse_list(marker_meth[m][i]) if m in marker_meth else []
                ak, av = anchor_kind[i], anchor_pv[i]
                # fallback: no diagnosis date -> anchor to first collection
                if av is None and ANCHOR_FALLBACK_TO_FIRST_ACQ:
                    base = next((d for d in dates if str(d).strip() not in EMPTY_TOKENS), None)
                    if base is not None:
                        ak, av = prep_anchor(base, bank)
                bucket = {}
                for k, v in enumerate(vals):
                    if str(v).strip() in EMPTY_TOKENS:
                        continue
                    d = dates[k] if k < len(dates) else ""
                    if str(d).strip() in EMPTY_TOKENS:
                        continue
                    w = assign_window(months_since(d, ak, av, bank))
                    if w is None:
                        continue
                    b = bucket.setdefault(w, {"v": [], "a": [], "m": []})
                    b["v"].append(v)
                    b["a"].append(d)
                    b["m"].append(meth[k] if k < len(meth) else "")
                for w, b in bucket.items():
                    cell(f"{m}_{w}", i, first_or_list(b["v"]))
                    cell(f"{m}_{w}_acq", i, first_or_list(b["a"]))
                    if any(str(x).strip() not in EMPTY_TOKENS for x in b["m"]):
                        cell(f"{m}_{w}_method", i, aligned_or_first(b["m"]))
        if skipped:
            _print(f"  [CSF] subjects left untimed (no own date, values ≠ shared dates): {skipped}")
        return pd.DataFrame(data, index=df.index) if data else None


    if csf_present:
        # how many markers actually carry any value (diagnose silent +0)
        nonempty = [m for m in csf_markers
                    if f"all_{m}_values" in df.columns and df[f"all_{m}_values"].notna().any()]
        fr = expand_csf(csf_markers, CSF_DATE_COL)
        if fr is not None and fr.shape[1] > 0:
            frames.append(fr)
            expanded.append(f"CSF[{len(nonempty)} markers]")
            added = fr.shape[1]
        else:
            added = 0
        df.drop(columns=[c for c in csf_cols if c in df.columns], inplace=True)
        gc.collect()
        if added > 0:
            _print(f"  expanded CSF ({len(nonempty)}/{len(csf_markers)} markers have data) (+{added} cols)")
        elif not nonempty:
            _print(f"  expanded CSF (+0 cols)  -> all {len(csf_markers)} *_values columns "
                  f"are empty in the master (RBM panel not populated)")
        else:
            _print(f"  expanded CSF (+0 cols)  -> {nonempty} have values but none could be "
                  f"anchored/dated (check all_pct_date_collection)")

    # diagnosis spread (value + acq)
    if "all_diag" in df.columns:
        L_dval = df["all_diag"].tolist()
        L_ddate = (df["all_diag_date"].tolist() if "all_diag_date" in df.columns else [""] * n)
        ddata = {}

        def dcell(col, i, value):
            arr = ddata.get(col)
            if arr is None:
                arr = [np.nan] * n
                ddata[col] = arr
            arr[i] = value

        for i in range(n):
            dv = parse_list(L_dval[i]); da = parse_list(L_ddate[i])
            ak, av, bank = anchor_kind[i], anchor_pv[i], banks[i]
            bucket = {}
            for j in range(max(len(dv), len(da))):
                v = dv[j] if j < len(dv) else ""
                a = da[j] if j < len(da) else ""
                if str(v).strip() in ("", "nan", "None"):
                    continue
                w = (assign_window(months_since(a, ak, av, bank))
                     if str(a).strip() not in ("", "nan", "None") else "pre")
                b = bucket.setdefault(w, {"v": [], "a": []})
                b["v"].append(v); b["a"].append(a)
            for w, b in bucket.items():
                dcell(f"diag_{w}", i, first_or_list(b["v"]))
                dcell(f"diag_{w}_acq", i, first_or_list(b["a"]))
        if ddata:
            frames.append(pd.DataFrame(ddata, index=df.index))

    # conversion timing (months from first diagnosis to first label change)
    def conversion_info(diag_cell, date_cell, bank, ak, av):
        diags = parse_list(diag_cell); dates = parse_list(date_cell)
        if not diags:
            return (np.nan, np.nan, 0)
        base = diags[0].split("+")[0].upper()
        for j in range(1, len(diags)):
            lab = str(diags[j]).strip()
            if lab in ("", "nan", "None"):
                continue
            if lab.split("+")[0].upper() != base:
                a = dates[j] if j < len(dates) else ""
                mo = months_since(a, ak, av, bank)
                return (round(mo, 1) if mo is not None else np.nan, lab, 1)
        return (np.nan, np.nan, 0)


    if "all_diag" in df.columns:
        L_dval = df["all_diag"].tolist()
        L_ddate = (df["all_diag_date"].tolist() if "all_diag_date" in df.columns else [""] * n)
        info = [conversion_info(v, d, b, ak, av) for v, d, b, ak, av in zip(L_dval, L_ddate, banks, anchor_kind, anchor_pv)]
        frames.append(pd.DataFrame({
            "months_to_conversion": [x[0] for x in info],
            "conversion_to":        [x[1] for x in info],
            "converted":            [x[2] for x in info],
        }, index=df.index))

    # baseline-anchor flag (diagnosis vs first_acq) for downstream stratification
    frames.append(pd.DataFrame({"baseline_anchor": baseline_anchor}, index=df.index))

    # =============================================================================
    # ASSEMBLE + SAVE
    # =============================================================================
    # Written in row blocks straight from the pieces: the full (subjects x ~6500
    # object columns) table is never assembled in RAM. The old
    # pd.concat([df] + frames, axis=1) made a second complete copy of everything
    # just before saving — that copy was the memory peak of this cell.
    gc.collect()
    pieces = [df] + frames
    n_cols_out = sum(p.shape[1] for p in pieces)
    SAVE_CHUNK_ROWS = 2000
    with open(OUT_PATH, "w", newline="", encoding="utf-8") as _fh:
        for _i0 in range(0, n, SAVE_CHUNK_ROWS):
            _part = pd.concat([p.iloc[_i0:_i0 + SAVE_CHUNK_ROWS] for p in pieces], axis=1)
            _part.to_csv(_fh, index=False, header=(_i0 == 0))
            del _part
    del pieces, frames
    gc.collect()

    _print("\n" + "=" * 64)
    _print(f"EXPANDED: {n} subjects x {n_cols_out} columns")
    _print(f"  features expanded: {expanded}")
    if skipped:
        _print(f"  features absent:   {skipped}")
    _print(f"  window grid: pre, lt{WINDOW_STEP}..lt{WINDOW_MAX}, gt{WINDOW_MAX}")
    _print(f"Saved -> {OUT_PATH}")
    write_provenance(OUT_PATH, step="expand", config=cfg, inputs=[MASTER_PATH],
                     parameters={"n_subjects": int(n), "n_columns": int(n_cols_out),
                                 "features_expanded": expanded, "features_absent": skipped})
    return OUT_PATH
