"""Data behind the cohort figures.

:class:`FigureData` holds

* ``subjects`` - one row per subject of the cohort (index ``uid``): databank,
  diagnosis label of the cohort table and its diagnostic group;
* ``sets`` - subjects that have each modality (cross-sectional figures);
* ``history`` - the subjects' rows of the expanded master (longitudinal
  figures): every diagnosis with its date and the 6-month windows of each
  measure. Without it the longitudinal figures use the conversion columns of
  the cohort table, and the Sankey / follow-up figures are skipped.

Time is handled as in the cohort tools (:func:`integration.common.vrs_time`):
calendar dates, or day offsets for OASIS and A4, compared only within a subject.
"""
from __future__ import annotations

import logging
import os
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..config import ConfigError, ProjectConfig, load_config
from ..integration import common as vrs
from .settings import DxGrouper

log = logging.getLogger(__name__)

DAYS_PER_MONTH = 30.4375
NO_VISIT = "No visit"
LOST = "Lost to follow-up"
OTHER = "Other / unknown"
_EMPTY = {"", "nan", "none", "nat", "na", "<na>"}


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------
def read_cohort(cohort: Any) -> pd.DataFrame:
    """Cohort table from a CSV path, a DataFrame or a :class:`CohortResult`."""
    if hasattr(cohort, "table") and isinstance(getattr(cohort, "table"), pd.DataFrame):
        df = cohort.table.copy()
    elif isinstance(cohort, pd.DataFrame):
        df = cohort.copy()
    else:
        df = pd.read_csv(cohort, low_memory=False, dtype={"uid": str, "ID": str})
    bank = "Databank" if "Databank" in df.columns else "dataset" if "dataset" in df.columns else None
    if bank is None:
        raise ConfigError("cohort table: no 'Databank' column")
    if bank != "Databank":
        df = df.rename(columns={bank: "Databank"})
    if "uid" not in df.columns:
        if "ID" not in df.columns:
            raise ConfigError("cohort table: needs 'uid' or 'Databank' + 'ID'")
        df["uid"] = df["Databank"].astype(str) + "_" + df["ID"].astype(str)
    df["uid"] = df["uid"].astype(str)
    return df.drop_duplicates("uid").reset_index(drop=True)


def _uids_of_file(path: str | os.PathLike) -> set[str]:
    head = pd.read_csv(path, nrows=0).columns
    use = ["uid"] if "uid" in head else [c for c in ("Databank", "dataset", "ID") if c in head]
    d = pd.read_csv(path, usecols=use, dtype=str, low_memory=False)
    if "uid" in d.columns:
        return set(d["uid"].dropna())
    bank = "Databank" if "Databank" in d.columns else "dataset"
    return set(d[bank].astype(str) + "_" + d["ID"].astype(str))


def _present(s: pd.Series) -> pd.Series:
    if s.dtype == object or str(s.dtype).startswith("string"):
        t = s.astype("string").str.strip().str.lower()
        return s.notna() & ~t.isin(list(_EMPTY))
    return s.notna()


def cohort_columns_for(df: pd.DataFrame, measure: str) -> list[str]:
    """Columns of a cohort table that tell whether ``measure`` is present."""
    m = re.escape(measure)
    for cands in ([measure, f"{measure}_value"],
                  [f"{measure}_date", f"{measure}_acq", f"{measure}_loniID", f"{measure}_id"],
                  [f"match_{measure}_date", f"match_{measure}_id"]):
        found = [c for c in cands if c in df.columns]
        if found:
            return found
    rx = re.compile(rf"^{m}_V\d+_[^_]+_date$")              # longitudinal cohorts
    return [c for c in df.columns if rx.match(c)]


def read_history(expanded: Any, uids: set[str], measures: list[str],
                 config: ProjectConfig) -> pd.DataFrame:
    """Rows of ``uids`` in the expanded master, with the windows of ``measures``."""
    keep = ["1st_diag", "Lst_diag", "1st_diag_date", "Lst_diag_date", "all_diag", "all_diag_date"]
    if isinstance(expanded, pd.DataFrame):
        df = expanded.copy()
        if "uid" not in df.columns:
            df["uid"] = df["Databank"].astype(str) + "_" + df["ID"].astype(str)
    else:
        vrs.configure(config)
        df = vrs.vrs_read_expanded(str(expanded), measures, keep_cols=keep)
    if "Databank" not in df.columns and "dataset" in df.columns:
        df = df.rename(columns={"dataset": "Databank"})
    df["uid"] = df["uid"].astype(str)
    return df[df["uid"].isin(uids)].drop_duplicates("uid").reset_index(drop=True)


# ---------------------------------------------------------------------------
# time and diagnosis histories
# ---------------------------------------------------------------------------
def _split(cell) -> list:
    if cell is None or (isinstance(cell, float) and np.isnan(cell)):
        return []
    s = str(cell).strip()
    if not s or s.lower() in _EMPTY:
        return []
    if s.startswith("[") and s.endswith("]"):
        return [None if x is None else str(x) for x in vrs.vrs_split_cell(s)]
    return [t.strip() for t in s.split("|")]


def _to_days(tokens: pd.Series, banks: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Time tokens -> (days, kind) with the cohort-tool rules (same index for both)."""
    tok = tokens.astype("string")
    tok = tok.where(~tok.str.strip().str.lower().isin(list(_EMPTY))).fillna("")
    return vrs.vrs_time(tok, banks)


def visit_table(history: pd.DataFrame, grouper: DxGrouper) -> pd.DataFrame:
    """One row per diagnosis visit: uid, Databank, k, label, group, t (days), kind."""
    rows = []
    for uid, bank, dv, dd in zip(history["uid"], history["Databank"],
                                 history.get("all_diag", pd.Series([None] * len(history))),
                                 history.get("all_diag_date", pd.Series([None] * len(history)))):
        labs, dates = _split(dv), _split(dd)
        for k in range(max(len(labs), len(dates))):
            lab = labs[k] if k < len(labs) else None
            if lab is None or str(lab).strip().lower() in _EMPTY:
                continue
            rows.append((uid, bank, k, str(lab).strip(), dates[k] if k < len(dates) else None))
    v = pd.DataFrame(rows, columns=["uid", "Databank", "k", "label", "date"])
    if v.empty:
        return v.assign(group=pd.Series(dtype=object), t=pd.Series(dtype=float),
                        kind=pd.Series(dtype=object))
    v["group"] = v["label"].map(grouper)
    v["t"], v["kind"] = _to_days(v["date"], v["Databank"])
    return v


def trajectories(visits: pd.DataFrame) -> pd.DataFrame:
    """Per subject: baseline group, first conversion to another group, follow-up.

    Baseline = first diagnosis of the history. A conversion is the first later
    visit whose group differs from the baseline group (visits whose label has
    no group are skipped). Months are counted from the baseline visit, within
    the same time axis.
    """
    out = []
    for uid, g in visits.sort_values(["uid", "k"]).groupby("uid", sort=False):
        b = g.iloc[0]
        base_g, t0, kind0 = b["group"], b["t"], b["kind"]
        same_axis = g[(g["kind"] == kind0) & g["t"].notna()] if pd.notna(t0) else g.iloc[0:0]
        fu = (same_axis["t"].max() - t0) / DAYS_PER_MONTH if len(same_axis) else np.nan
        conv = g.iloc[1:]
        conv = conv[conv["group"].notna() & (conv["group"] != base_g)] if base_g is not None else conv.iloc[0:0]
        if len(conv):
            c = conv.iloc[0]
            mo = ((c["t"] - t0) / DAYS_PER_MONTH
                  if pd.notna(c["t"]) and pd.notna(t0) and c["kind"] == kind0 else np.nan)
            out.append((uid, b["Databank"], b["label"], base_g, 1, c["label"], c["group"], mo, fu, len(g)))
        else:
            out.append((uid, b["Databank"], b["label"], base_g, 0, None, None, np.nan, fu, len(g)))
    return pd.DataFrame(out, columns=["uid", "Databank", "baseline_label", "baseline_group", "converted",
                                      "conversion_to_label", "conversion_to_group",
                                      "months_to_conversion", "followup_months", "n_visits"]
                        ).set_index("uid")


def trajectories_from_cohort(df: pd.DataFrame, grouper: DxGrouper) -> pd.DataFrame:
    """Same table from the conversion columns of a cohort table (no expanded master).

    Uses 1st_diag, conversion_to, converted, months_to_conversion and the
    first / last diagnosis dates. A conversion inside one group (e.g. LBD ->
    LBD+PD) is not counted.
    """
    need = {"1st_diag", "converted"}
    if not need.issubset(df.columns):
        return pd.DataFrame()
    d = df.set_index("uid")
    base = d["1st_diag"].map(grouper)
    to = d.get("conversion_to", pd.Series(index=d.index, dtype=object)).map(grouper)
    conv = pd.to_numeric(d["converted"], errors="coerce").fillna(0).astype(int)
    conv = conv.where(~(to.isna() | (to == base)), 0)
    fu = pd.Series(np.nan, index=d.index)
    if {"1st_diag_date", "Lst_diag_date"}.issubset(d.columns):
        t0, k0 = _to_days(d["1st_diag_date"], d["Databank"])
        t1, k1 = _to_days(d["Lst_diag_date"], d["Databank"])
        fu = ((t1 - t0) / DAYS_PER_MONTH).where(k0 == k1)
    mo = pd.to_numeric(d.get("months_to_conversion", pd.Series(index=d.index, dtype=float)),
                       errors="coerce").where(conv == 1)
    return pd.DataFrame({"Databank": d["Databank"], "baseline_label": d["1st_diag"],
                         "baseline_group": base, "converted": conv,
                         "conversion_to_label": d.get("conversion_to").where(conv == 1)
                         if "conversion_to" in d.columns else None,
                         "conversion_to_group": to.where(conv == 1),
                         "months_to_conversion": mo, "followup_months": fu,
                         "n_visits": pd.to_numeric(d.get("nb_of_visits", np.nan), errors="coerce")})


def window_start(window: str, step: int) -> float | None:
    """'lt6' -> 0, 'lt12' -> 6 ... (months since baseline); pre / gt -> None."""
    m = re.fullmatch(r"lt(\d+)", str(window))
    return float(int(m.group(1)) - step) if m else None


# ---------------------------------------------------------------------------
# FigureData
# ---------------------------------------------------------------------------
@dataclass
class FigureData:
    settings: dict
    grouper: DxGrouper
    cohort: pd.DataFrame
    subjects: pd.DataFrame                 # index uid: Databank, dx_label, dx
    reference: str | None
    sets: dict[str, set] = field(default_factory=dict)
    modality_source: dict[str, str] = field(default_factory=dict)
    history: pd.DataFrame | None = None
    notes: list[str] = field(default_factory=list)
    n_removed: int = 0
    removed_labels: Counter = field(default_factory=Counter)
    _visits: pd.DataFrame | None = None
    _traj: pd.DataFrame | None = None

    # -- basic views ------------------------------------------------------
    @property
    def order(self) -> list[str]:
        return self.grouper.order

    @property
    def colors(self) -> dict[str, str]:
        return self.grouper.colors

    @property
    def U(self) -> set:
        return set(self.subjects.index)

    @property
    def N(self) -> int:
        return len(self.subjects)

    @property
    def ref_label(self) -> str:
        return self.reference or "cohort"

    def modalities(self) -> list[str]:
        """Modalities found in this cohort (reference first)."""
        mods = [m for m in self.settings["modalities"] if m in self.sets]
        if self.reference in mods:
            mods.remove(self.reference)
            mods.insert(0, self.reference)
        return mods

    def comp_pct(self, uids) -> pd.Series:
        d = self.subjects.loc[self.subjects.index.intersection(list(uids)), "dx"]
        return d.value_counts(normalize=True).reindex(self.order).fillna(0) * 100

    def steps(self) -> list[tuple[str, set]]:
        """Composite-cohort steps: modalities (not the reference), bundles merged."""
        mods = [m for m in self.modalities() if m != self.reference]
        bundles = {b: [x for x in mem] for b, mem in (self.settings.get("bundles") or {}).items()
                   if all(x in self.sets for x in mem)}
        used, out = set(), []
        for m in mods:
            if m in used:
                continue
            b = next((b for b, mem in bundles.items() if m in mem), None)
            if b:
                s = set.intersection(*(self.sets[x] for x in bundles[b]))
                out.append((b, s)); used.update(bundles[b])
            else:
                out.append((m, self.sets[m])); used.add(m)
        return out

    # -- longitudinal views -------------------------------------------------
    @property
    def visits(self) -> pd.DataFrame | None:
        if self._visits is None and self.history is not None:
            self._visits = visit_table(self.history, self.grouper)
        return self._visits

    @property
    def traj(self) -> pd.DataFrame:
        if self._traj is None:
            coh = self.cohort[self.cohort["uid"].isin(self.subjects.index)]
            if self.visits is not None and len(self.visits):
                th = trajectories(self.visits)
                rest = coh[~coh["uid"].isin(th.index)]
                tc = trajectories_from_cohort(rest, self.grouper) if len(rest) else pd.DataFrame()
                self._traj = pd.concat([th, tc]) if len(tc) else th
                self._traj_source = "expanded master (all_diag / all_diag_date)"
                if len(rest):
                    self._traj_source += (f"; {len(rest):,} subjects not in the expanded master: "
                                          + ("cohort columns" if len(tc) else "left out"))
            else:
                self._traj = trajectories_from_cohort(coh, self.grouper)
                self._traj_source = "cohort table (1st_diag / conversion_to / months_to_conversion)"
        return self._traj

    def states_by_bin(self) -> pd.DataFrame:
        """Diagnostic state of each subject in each time bin (rows uid, columns bin start)."""
        L = self.settings["longitudinal"]
        W, MAX = int(L["window_months"]), int(L["max_months"])
        bins = list(range(0, MAX + 1, W))
        v = self.visits
        tr = self.traj
        tr = tr[tr["baseline_group"].notna()]
        if v is None or v.empty or tr.empty:
            return pd.DataFrame(columns=bins)
        tr = tr[tr.index.isin(v["uid"])]               # subjects with a diagnosis history only
        v = v[v["uid"].isin(tr.index)].copy()
        base = v.sort_values("k").groupby("uid").first()
        v["t0"] = v["uid"].map(base["t"])
        v["kind0"] = v["uid"].map(base["kind"])
        v = v[(v["kind"] == v["kind0"]) & v["t"].notna() & v["t0"].notna()]
        v["mo"] = (v["t"] - v["t0"]) / DAYS_PER_MONTH
        v = v[v["mo"] >= -0.5]
        v["bin"] = (np.floor(v["mo"].clip(lower=0) / W) * W).astype(int)
        v["state"] = v["group"].fillna(OTHER)
        last = v.sort_values(["uid", "k"]).groupby(["uid", "bin"])["state"].last().unstack()
        last = last.reindex(index=tr.index, columns=bins)
        last[0] = tr["baseline_group"]
        last_seen = v.groupby("uid")["bin"].max().reindex(tr.index).fillna(0)
        out = last.copy().astype(object)
        for b in bins[1:]:
            empty = out[b].isna()
            out.loc[empty & (last_seen > b), b] = NO_VISIT
            out.loc[empty & (last_seen <= b), b] = LOST
        return out

    def availability(self, label: str) -> pd.DataFrame | None:
        """Subjects with ``label`` in each time bin, by baseline group (rows bin, cols group)."""
        if self.history is None:
            return None
        spec = self.settings["modalities"][label]
        measure = spec.get("measure")
        if not measure:
            return None
        L = self.settings["longitudinal"]
        W, MAX = int(L["window_months"]), int(L["max_months"])
        step = int(self.settings.get("integration_window", 6))
        ev = vrs.vrs_events(self.history, measure, label=label)
        if ev.empty:
            return None
        ev["start"] = ev["window"].map(lambda w: window_start(w, step))
        ev = ev[ev["start"].notna()]
        ev["bin"] = (np.floor(ev["start"] / W) * W).astype(int)
        ev = ev[ev["bin"] <= MAX]
        base = self.traj["baseline_group"]
        ev["group"] = ev["uid"].map(base)
        ev = ev[ev["group"].notna()]
        tab = ev.drop_duplicates(["uid", "bin"]).groupby(["bin", "group"]).size().unstack(fill_value=0)
        return tab.reindex(index=range(0, MAX + 1, W), columns=self.order).fillna(0).astype(int)


def build_figure_data(cohort: Any, *, settings: dict, expanded: Any = None,
                      dx_column: str | None = None, config: ProjectConfig | None = None,
                      need_history: bool = True) -> FigureData:
    """Collect everything the figures need from a cohort table (and the expanded master)."""
    cfg = config or load_config()
    vrs.configure(cfg)
    grouper = DxGrouper(settings)
    df = read_cohort(cohort)
    notes: list[str] = []
    dx_col = dx_column or next((c for c in ("diag_at_anchor", "1st_diag", "Lst_diag") if c in df.columns), None)
    if dx_col is None:
        raise ConfigError("cohort table: no diagnosis column (diag_at_anchor / 1st_diag); give dx_column")
    notes.append(f"Diagnosis = {dx_col}.")
    subj = pd.DataFrame({"Databank": df["Databank"].values, "dx_label": df[dx_col].values},
                        index=df["uid"].values)
    subj.index.name = "uid"
    subj["dx"] = subj["dx_label"].map(grouper)

    mods: dict = settings["modalities"]
    history = None
    need_exp = [lab for lab, s in mods.items()
                if s.get("source") in ("ever", "baseline")
                or (not s.get("file") and not s.get("column") and s.get("measure")
                    and not cohort_columns_for(df, s["measure"]))]
    measures = sorted({s["measure"] for s in mods.values() if s.get("measure")})
    if expanded is not None and (need_history or need_exp):
        history = read_history(expanded, set(subj.index), measures + ["diag"], cfg)
    elif expanded is None:
        notes.append("No expanded master given: Sankey and follow-up figures are skipped; "
                     "conversions come from the cohort table.")

    sets: dict[str, set] = {}
    source: dict[str, str] = {}
    for lab, s in mods.items():
        if s.get("file"):
            sets[lab] = _uids_of_file(s["file"]) & set(subj.index)
            source[lab] = f"file {Path(s['file']).name}"
            continue
        cols = ([s["column"]] if s.get("column") else
                cohort_columns_for(df, s["measure"]) if s.get("measure") else [])
        cols = [c for c in cols if c in df.columns]
        if cols and s.get("source", "cohort") == "cohort":
            pres = np.zeros(len(df), dtype=bool)
            for c in cols:
                pres |= _present(df[c]).to_numpy()
            sets[lab] = set(df.loc[pres, "uid"])
            source[lab] = "cohort column " + " / ".join(cols[:3]) + (" ..." if len(cols) > 3 else "")
            continue
        if history is not None and s.get("measure"):
            ev = vrs.vrs_events(history, s["measure"], label=lab)
            if s.get("source") == "baseline":
                ev = ev[ev["window"] == f"lt{int(settings.get('integration_window', 6))}"]
            if len(ev):
                sets[lab] = set(ev["uid"])
                source[lab] = (f"expanded master, {s['measure']} "
                               + ("first window" if s.get("source") == "baseline" else "any time"))
                continue
        notes.append(f"Modality {lab!r} not found in this cohort - left out.")

    ref = settings.get("reference")
    if ref is not None and ref not in sets:
        notes.append(f"Reference {ref!r} not found: percentages are of all cohort subjects.")
        ref = None
    keep = subj.index
    if ref is not None:
        keep = keep[keep.isin(list(sets[ref]))]
    subj = subj.loc[keep]
    removed = subj["dx"].isna()
    removed_labels = Counter(subj.loc[removed, "dx_label"].fillna("missing").astype(str))
    n_removed = int(removed.sum())
    subj = subj[~removed]
    U = set(subj.index)
    sets = {k: v & U for k, v in sets.items()}
    if ref is not None:
        sets[ref] = U
    fd = FigureData(settings=settings, grouper=grouper, cohort=df, subjects=subj, reference=ref,
                    sets=sets, modality_source=source, history=history, notes=notes,
                    n_removed=n_removed, removed_labels=removed_labels)
    log.info("figures: %d subjects (%d without a diagnostic group removed); modalities: %s",
             fd.N, n_removed, ", ".join(f"{m}={len(sets[m])}" for m in fd.modalities()))
    return fd
