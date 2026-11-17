"""Draw every cohort figure and write the images, source tables and a summary."""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from ..config import ConfigError, ProjectConfig, load_config
from .cross_sectional import CROSS_SECTIONAL
from .data import FigureData, build_figure_data
from .longitudinal import LONGITUDINAL
from .settings import load_figure_settings
from .style import Panel, apply_style, plt

log = logging.getLogger(__name__)

#: figure key -> (type of figure, what it shows, purpose); used for the summary file
CATALOGUE = {
    "modality_connectome_graph": (
        "Network graph with a diagnosis ring around each node",
        "Modalities as nodes (size = subjects), edges = shared subjects, ring = diagnostic mix",
        "Which modalities coexist in the same subjects and how the diagnostic profile changes between them"),
    "modality_overlap_matrix": (
        "Annotated lower-triangle heatmap", "n and % of shared subjects for each pair of modalities",
        "Exact overlap values behind the connectome"),
    "composite_cohort_lattice": (
        "Hasse diagram with a diagnosis strip in each box",
        "Composite cohorts adding one modality at a time, with n, % and diagnostic mix",
        "Sample cost of requiring each extra modality"),
    "composite_cohort_dx_stack": (
        "100 % stacked horizontal bars", "Diagnosis distribution in each composite cohort",
        "Selection bias introduced by each combination of modalities"),
    "dx_group_counts_log": (
        "Vertical bars, log scale", "Subjects per diagnostic group", "Group imbalance, including rare groups"),
    "dx_subgroup_labels": (
        "Horizontal bars coloured by group", "Original labels inside groups that pool several labels",
        "Transparency of the diagnostic grouping"),
    "dx_by_databank_heatmap": (
        "Contingency heatmap (log shading, totals row)", "Diagnostic group x source databank",
        "Which databanks contribute each diagnosis"),
    "modality_venn": (
        "Two-set Venn diagram (area-proportional)", "Subjects with only A, only B or both",
        "Overlap of two measures (default MMSE / MoCA) to justify harmonisation or a choice"),
    "venn_region_by_databank": (
        "100 % stacked bars (grey scale)", "Source databank in each Venn region",
        "Whether the use of A or B depends on the protocol of each databank"),
    "venn_region_by_dx": (
        "100 % stacked bars", "Diagnosis distribution in each Venn region",
        "Whether the A and B groups differ clinically"),
    "converters_transition_upset": (
        "UpSet plot with arrows", "Diagnostic transitions of converters by type, median (IQR) months to conversion",
        "Conversion trajectories and time to event"),
    "non_converters_stability_upset": (
        "UpSet plot with hatching", "Non-converters per diagnosis, median (IQR) follow-up, single assessments",
        "Diagnostic stability and real follow-up per group"),
    "dx_transition_sankey": (
        "Longitudinal Sankey / alluvial", "Diagnosis in each time window since baseline, including no visit / lost",
        "Stability, conversions, reversions and attrition over time"),
    "followup_availability_area_by_dx": (
        "Stacked area per measure", "Subjects with each measure in each time window, by baseline diagnosis",
        "Attrition and feasible time horizon of longitudinal analyses"),
}
ALL_FIGURES = {**CROSS_SECTIONAL, **LONGITUDINAL}
CROSS_KEYS = list(CROSS_SECTIONAL)
LONG_KEYS = list(LONGITUDINAL)


@dataclass
class FigureReport:
    """Result of :func:`make_cohort_figures`."""
    out_dir: Path
    files: dict[str, list[Path]] = field(default_factory=dict)
    skipped: dict[str, str] = field(default_factory=dict)
    tables_path: Path | None = None
    summary_path: Path | None = None
    data: FigureData | None = None

    def __repr__(self) -> str:
        return (f"FigureReport({len(self.files)} figures in {self.out_dir}, "
                f"{len(self.skipped)} skipped)")


def _select(only) -> list[str]:
    if only in (None, "all"):
        return list(ALL_FIGURES)
    if only == "cross_sectional":
        return CROSS_KEYS
    if only == "longitudinal":
        return LONG_KEYS
    keys = [only] if isinstance(only, str) else list(only)
    alias = {"mmse_moca_venn": "modality_venn", "dx_transition_sankey_6m": "dx_transition_sankey"}
    groups = {"all": list(ALL_FIGURES), "cross_sectional": CROSS_KEYS, "longitudinal": LONG_KEYS}
    keys = list(dict.fromkeys(x for k in keys for x in groups.get(k, [alias.get(k, k)])))
    bad = [k for k in keys if k not in ALL_FIGURES]
    if bad:
        raise ConfigError(f"unknown figure(s) {bad}; choose from {list(ALL_FIGURES)} "
                          "or 'cross_sectional' / 'longitudinal' / 'all'")
    return keys


def make_cohort_figures(cohort: Any, out_dir: str | os.PathLike, *, expanded: Any = None,
                        spec: str | os.PathLike | dict | None = None, figures="all",
                        dx_column: str | None = None, config: ProjectConfig | None = None,
                        show: bool = False) -> FigureReport:
    """Draw the cohort figures.

    Parameters
    ----------
    cohort:
        Output of a cohort tool (CSV path, DataFrame or ``CohortResult``):
        cross-sectional (anchor / multimodal) or longitudinal. One row per
        subject; needs ``uid`` (or ``Databank`` + ``ID``) and a diagnosis column.
    out_dir:
        Folder for the images, ``figure_source_tables.xlsx`` and ``FIGURES.md``.
    expanded:
        Expanded master (``MASTER_expanded_240m.csv``, path or DataFrame).
        Needed by the Sankey and follow-up figures and for modalities that are
        not columns of the cohort table; the conversion figures use it when given.
    spec:
        YAML file or dict overriding ``configs/figures.yaml`` (modalities,
        reference, bundles, venn, diagnostic groups, time windows, style).
    figures:
        ``"all"``, ``"cross_sectional"``, ``"longitudinal"`` or a list of figure keys.
    dx_column:
        Diagnosis column of the cohort table (default ``diag_at_anchor``, else ``1st_diag``).
    show:
        Also display the figures (notebooks).

    Returns
    -------
    FigureReport
        Paths of the files written, figures skipped (with the reason) and the data used.
    """
    cfg = config or load_config()
    settings = load_figure_settings(spec, cfg)
    keys = _select(figures)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    apply_style(settings)
    need_hist = any(k in LONG_KEYS for k in keys)
    fd = build_figure_data(cohort, settings=settings, expanded=expanded, dx_column=dx_column,
                           config=cfg, need_history=need_hist)
    rep = FigureReport(out_dir=out, data=fd)
    formats = settings["style"].get("formats") or ["png"]
    tables: dict[str, pd.DataFrame] = {}
    panels: list[Panel] = []
    for k in keys:
        try:
            p = ALL_FIGURES[k](fd)
        except Exception as exc:                       # one broken figure must not stop the others
            log.exception("figure %s failed", k)
            rep.skipped[k] = f"error: {exc}"
            continue
        panels.append(p)
        if p.fig is None:
            rep.skipped[p.name] = p.note
            log.info("  [skip] %s: %s", p.name, p.note)
            continue
        paths = []
        for fmt in formats:
            path = out / f"{p.name}.{fmt}"
            p.fig.savefig(path, bbox_inches="tight", dpi=int(settings["style"].get("dpi", 300)))
            paths.append(path)
        rep.files[p.name] = paths
        log.info("  [ok]   %s", paths[0].name)
        if p.table is not None:
            tables[p.name] = p.table
        for tk, tv in p.extra_tables.items():
            tables[tk] = tv
        if show:
            plt.show()
        plt.close(p.fig)
    if tables:
        rep.tables_path = out / "figure_source_tables.xlsx"
        with pd.ExcelWriter(rep.tables_path) as w:
            for name, t in tables.items():
                t.to_excel(w, sheet_name=name[:31])
    rep.summary_path = out / "FIGURES.md"
    rep.summary_path.write_text(_summary(fd, rep, panels), encoding="utf-8")
    return rep


def _summary(fd: FigureData, rep: FigureReport, panels: list[Panel]) -> str:
    lines = ["# Cohort figures", "",
             f"Subjects: {fd.N:,} ({fd.ref_label}); {fd.n_removed:,} without a diagnostic group left out.", ""]
    lines += ["| figure | type | shows | purpose | file |", "|---|---|---|---|---|"]
    for p in panels:
        key = p.name if p.name in CATALOGUE else ("modality_venn" if p.name.endswith("_venn")
                                                  else "dx_transition_sankey" if "sankey" in p.name else p.name)
        typ, shows, purpose = CATALOGUE.get(key, ("", "", ""))
        f = rep.files.get(p.name)
        lines.append(f"| `{p.name}` | {typ} | {shows} | {purpose} | "
                     + (f"`{f[0].name}`" if f else f"not drawn: {p.note}") + " |")
    lines += ["", "## Modalities", "", "| modality | subjects | found in |", "|---|---|---|"]
    for m in fd.modalities():
        lines.append(f"| {m} | {len(fd.sets[m]):,} | {fd.modality_source.get(m, 'reference')} |")
    if fd.notes:
        lines += ["", "## Notes", ""] + [f"- {n}" for n in fd.notes]
    if fd.removed_labels:
        lines += ["", "Labels without a group: "
                  + ", ".join(f"{k} ({v})" for k, v in fd.removed_labels.most_common())]
    return "\n".join(lines) + "\n"
