"""Figures to describe a cohort built with the cohort tools.

Typical use::

    from inbrain_dementia_data_analyser.figures import make_cohort_figures
    rep = make_cohort_figures("my_cohort_summary.csv", "figures_out",
                              expanded=".../Integrated_ALL/MASTER_expanded_240m.csv")
    print(rep.files, rep.skipped)

Cross-sectional figures (modalities, composite cohorts, diagnostic mix, Venn)
use the cohort table; longitudinal figures (conversions, Sankey, follow-up)
use the subjects' diagnosis history in the expanded master. Defaults live in
``configs/figures.yaml``; a spec (YAML or dict) overrides them.
Requires matplotlib.
"""
from .report import ALL_FIGURES, CATALOGUE, FigureReport, make_cohort_figures
from .settings import DxGrouper, load_figure_settings

__all__ = ["ALL_FIGURES", "CATALOGUE", "DxGrouper", "FigureReport", "load_figure_settings",
           "make_cohort_figures"]
