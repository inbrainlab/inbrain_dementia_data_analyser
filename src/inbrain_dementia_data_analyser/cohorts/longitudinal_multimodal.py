"""Longitudinal multimodal cohorts: a visit schedule where every visit carries a set of measures.

Definition used (to be confirmed by the lab):

1. **Events.** Every acquisition of the schedule modality (e.g. ``MRI_T1w``)
   is an anchor event. An event is *complete* when every ``required`` measure
   (clinical scales, PET, plasma, CSF, other MRI ...) has an acquisition within
   +/- ``margin_months`` of it (per-measure margins allowed). This is the
   anchor pairing of :func:`~.anchor.create_anchor_cohort` with
   ``anchor_events="all"``.
2. **Sessions.** Complete events of a subject closer than
   ``session_window_days`` are one visit; the event with most optional
   matches represents it (ties -> the earliest).
3. **Schedule.** Visit targets are ``visit_targets_months`` after the subject's
   first complete visit; each target takes the closest visit within
   +/- ``tolerance_months``, with the same globally optimal assignment and the
   same non-overlapping tolerances as :func:`~.longitudinal.create_longitudinal_cohort`.
4. Subjects with at least ``min_visits`` matched targets are kept.

Because a visit only counts when it is complete, a subject whose second MRI
lacks a CDR within the margin loses that visit even if the MRI itself is on
schedule.

Outputs (one row per subject x matched visit, i.e. long format):
``<prefix>_longitudinal_multimodal.csv`` with ``visit``, ``target_months``,
``offset_days`` and every column of the matched anchor event (values, dates,
image IDs, deltas of each measure, diagnosis and age at the visit), plus
``<prefix>_longitudinal_multimodal_loniIDs.csv/.txt``. The complete events
themselves are written as ``<prefix>_events_pairs.csv`` (anchor pairing output).
"""
from __future__ import annotations

import logging
import os
from typing import Sequence

import numpy as np
import pandas as pd

from ..config import ProjectConfig, load_config
from ..integration.common import configure as configure_common, vrs_time
from ..provenance import write_provenance
from .anchor import DEFAULT_KEEP_COLUMNS, create_anchor_cohort
from .longitudinal import _HAS_SCIPY, _target_offsets_days, _tolerances_days
from .outputs import CohortResult
from ..logging_utils import quiet_pandas_warnings

log = logging.getLogger(__name__)

try:
    from scipy.optimize import linear_sum_assignment
except ImportError:  # pragma: no cover - scipy is a declared dependency
    linear_sum_assignment = None

_BIG = 10 ** 9


def _sessions(ev: pd.DataFrame, window_days: float) -> list[pd.Series]:
    """Group the events of one subject x time axis into sessions; pick one event each."""
    ev = ev.sort_values("t")
    groups, cur, start = [], [], None
    for _, row in ev.iterrows():
        if cur and row["t"] - start > window_days:
            groups.append(cur)
            cur = []
        if not cur:
            start = row["t"]
        cur.append(row)
    if cur:
        groups.append(cur)
    out = []
    for grp in groups:
        g = pd.DataFrame(grp)
        g = g.sort_values(["n_optional_matched", "t"], ascending=[False, True])
        out.append(g.iloc[0])
    return out


def _assign(offsets: list[float], target_days: list[int], tols: list[int]) -> dict[int, int]:
    """Target index -> session index (globally optimal; greedy without scipy)."""
    n_t, n_s = len(target_days), len(offsets)
    cost = np.full((n_t, n_s), float(_BIG))
    for i in range(n_t):
        for j in range(n_s):
            d = abs(offsets[j] - target_days[i])
            if d <= tols[i]:
                cost[i, j] = d
    assign: dict[int, int] = {}
    if linear_sum_assignment is not None:
        rows, cols = linear_sum_assignment(cost)
        for i, j in zip(rows, cols):
            if cost[i, j] < _BIG:
                assign[int(i)] = int(j)
        return assign
    pairs = sorted((cost[i, j], i, j) for i in range(n_t) for j in range(n_s) if cost[i, j] < _BIG)
    used_t, used_s = set(), set()
    for _, i, j in pairs:
        if i in used_t or j in used_s:
            continue
        assign[i] = j
        used_t.add(i)
        used_s.add(j)
    return assign


@quiet_pandas_warnings
def create_longitudinal_multimodal_cohort(
        input_path: str | os.PathLike, *,
        schedule_modality: str = "MRI_T1w",
        required: Sequence[str] = ("CDR",),
        optional: Sequence[str] = (),
        margin_months: float = 6,
        margin_months_by_feature: dict | None = None,
        visit_targets_months: Sequence[int] = (0, 12, 24),
        tolerance_months: float = 6,
        min_visits: int = 2,
        session_window_days: float = 31,
        clamp_windows: bool = True,
        keep_columns: Sequence[str] = DEFAULT_KEEP_COLUMNS,
        select_databanks: str | Sequence[str] = "ALL",
        cohort_filter_path: str | os.PathLike | None = None,
        output_prefix: str = "longitudinal_multimodal",
        out_dir: str | os.PathLike | None = None,
        config: ProjectConfig | None = None) -> CohortResult:
    """Subjects followed over time with a complete set of measures at every visit.

    Parameters
    ----------
    input_path:
        ``MASTER_expanded_240m.csv``.
    schedule_modality:
        Measure (stem or group) whose acquisitions define the visits.
    required:
        Measures that must be present within the margin at a visit (e.g.
        ``["CDR", "PET_ABETA", "PLASMA_PTAU217"]``).
    optional:
        Measures reported when present within the margin.
    margin_months, margin_months_by_feature:
        +/- months between the visit and each measure (per-measure override).
    visit_targets_months, tolerance_months, clamp_windows:
        Visit schedule, counted from the first complete visit.
    min_visits:
        Minimum number of matched targets to keep a subject.
    session_window_days:
        Complete events closer than this are one visit.
    keep_columns, select_databanks, cohort_filter_path:
        As in :func:`~.anchor.create_anchor_cohort`.
    output_prefix, out_dir:
        Output names / folder (default: the folder of ``input_path``).
    config:
        Project configuration.

    Returns
    -------
    CohortResult
        ``table``: one row per subject x matched visit; ``images``: matched
        image IDs; ``summary``: one row per subject (number of matched visits,
        follow-up span).
    """
    cfg = config or load_config()
    configure_common(cfg)
    out_dir = str(out_dir or os.path.dirname(str(input_path)))
    targets = [int(t) for t in visit_targets_months]
    tol_days = int(round(tolerance_months * cfg.integration["time_axis"].get("days_per_month", 30.4375)))

    events = create_anchor_cohort(
        input_path, anchor=schedule_modality, required=required, optional=optional,
        margin_months=margin_months, margin_months_by_feature=margin_months_by_feature,
        anchor_events="all", keep_columns=keep_columns, select_databanks=select_databanks,
        cohort_filter_path=cohort_filter_path, output_prefix=f"{output_prefix}_events",
        out_dir=out_dir, config=cfg)
    ev = events.table.copy()
    if ev.empty:
        raise ValueError("No complete visit found: relax 'required' or the margins.")
    bank_col = "Databank" if "Databank" in ev.columns else None
    bank = ev[bank_col] if bank_col else pd.Series("", index=ev.index)
    # day codes such as OASIS 'd0129' -> 129 (days since entry)
    acq = ev[f"{schedule_modality}_acq"].astype("string").str.strip() \
        .str.replace(r"^d(?:ay)?[_ ]?(\d+)$", r"\1", regex=True)
    ev["t"], ev["_kind"] = vrs_time(acq, bank)
    ev = ev[ev["t"].notna()]

    rows, subj_rows = [], []
    for (uid, kind), grp in ev.groupby(["uid", "_kind"], sort=False):
        sess = _sessions(grp, session_window_days)
        base = sess[0]
        target_days = _target_offsets_days({"t": base["t"]}, targets, kind)
        tols = _tolerances_days(target_days, tol_days, clamp_windows)
        offsets = [s["t"] - base["t"] for s in sess]
        assign = _assign(offsets, target_days, tols)
        if len(assign) < min_visits:
            continue
        prev = next((r for r in subj_rows if r["uid"] == uid), None)
        if prev is not None and prev["n_visits_matched"] >= len(assign):
            continue                      # another time axis already gave more visits
        if prev is not None:
            subj_rows.remove(prev)
            rows = [r for r in rows if r["uid"] != uid]
        for i, j in sorted(assign.items()):
            r = sess[j].to_dict()
            r.update({"visit": f"V{i + 1}_{targets[i]}m", "target_months": targets[i],
                      "offset_days": int(round(offsets[j]))})
            rows.append(r)
        subj_rows.append({"uid": uid, "time_axis": kind, "n_visits_matched": len(assign),
                          "n_complete_sessions": len(sess),
                          "followup_days": int(round(max(offsets[j] for j in assign.values()))),
                          "visits": "|".join(f"V{i + 1}_{targets[i]}m" for i in sorted(assign))})

    table = pd.DataFrame(rows)
    summary = pd.DataFrame(subj_rows)
    if not table.empty:
        table = table.drop(columns=["t", "_kind"], errors="ignore")
        lead = [c for c in ("uid", "Databank", "ID", "visit", "target_months", "offset_days")
                if c in table.columns]
        table = table[lead + [c for c in table.columns if c not in lead]]
        table = table.sort_values(["uid", "target_months"]).reset_index(drop=True)
    imgs = events.images
    if not imgs.empty and not table.empty:
        key = table[["uid", "anchor_event_n", "visit"]]
        imgs = imgs.merge(key, on=["uid", "anchor_event_n"], how="inner")
    else:
        imgs = pd.DataFrame()

    os.makedirs(out_dir, exist_ok=True)
    files = {"table": os.path.join(out_dir, f"{output_prefix}_longitudinal_multimodal.csv"),
             "summary": os.path.join(out_dir, f"{output_prefix}_longitudinal_multimodal_subjects.csv"),
             "events": events.files.get("pairs", "")}
    table.to_csv(files["table"], index=False)
    summary.to_csv(files["summary"], index=False)
    if not imgs.empty:
        files["images"] = os.path.join(out_dir, f"{output_prefix}_longitudinal_multimodal_loniIDs.csv")
        files["image_ids"] = os.path.join(out_dir, f"{output_prefix}_longitudinal_multimodal_loniIDs.txt")
        imgs.to_csv(files["images"], index=False)
        ids = imgs["loniID"].astype(str).str.strip().drop_duplicates().tolist()
        with open(files["image_ids"], "w") as fh:
            fh.write(",".join(ids))
    n_subj = int(summary["uid"].nunique()) if not summary.empty else 0
    log.info("LONGITUDINAL MULTIMODAL: %d subjects, %d visits (>= %d of %s) -> %s",
             n_subj, len(table), min_visits, targets, files["table"])
    if not summary.empty and "Databank" in table.columns:
        log.info("%s", table.drop_duplicates("uid")["Databank"].value_counts().to_string())
    write_provenance(files["table"], step="longitudinal_multimodal_cohort", config=cfg,
                     inputs=[str(input_path)], table=table,
                     parameters={"schedule_modality": schedule_modality, "required": list(required),
                                 "optional": list(optional), "margin_months": margin_months,
                                 "margin_months_by_feature": margin_months_by_feature or {},
                                 "visit_targets_months": targets, "tolerance_months": tolerance_months,
                                 "min_visits": min_visits, "session_window_days": session_window_days,
                                 "scipy": _HAS_SCIPY})
    return CohortResult(table=table, images=imgs, summary=summary, files=files)
