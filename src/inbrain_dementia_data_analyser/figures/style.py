"""Shared look of the cohort figures (matplotlib)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - depends on the environment
    raise ImportError("The figures need matplotlib: pip install matplotlib "
                      "(or pip install 'inbrain_dementia_data_analyser[figures]')") from exc

INK, INK2, GRID, EDGE, OFF = "#0b0b0b", "#52514e", "#e6e5e0", "#8a8984", "#dcdbd5"
NO_VISIT_COLOR, LOST_COLOR, OTHER_COLOR = "#c9c8c1", "#ecebe6", "#a9a8a1"
BLUES = ["#c9dcf5", "#2a78d6", "#0d3a75"]
TRANSITION_COLORS = {"Progression": "#c2410c", "Reversion": "#2a78d6", "Change of aetiology": "#6250d6"}


def apply_style(settings: dict) -> None:
    st = settings.get("style", {})
    plt.rcParams.update({
        "font.family": st.get("font", "DejaVu Sans"), "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": EDGE,
        "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "savefig.dpi": int(st.get("dpi", 300)), "figure.facecolor": "white",
        "savefig.facecolor": "white"})


def blue_cmap():
    return matplotlib.colors.LinearSegmentedColormap.from_list("inbrain_blue", BLUES)


def greys(n: int) -> list[str]:
    """``n`` grey shades from dark to light."""
    if n <= 0:
        return []
    vals = np.linspace(0.17, 0.88, n) if n > 1 else [0.4]
    return [matplotlib.colors.to_hex((v, v * 0.995, v * 0.97)) for v in vals]


@dataclass
class Panel:
    """One figure: matplotlib figure (None when it cannot be drawn), source table, note."""
    name: str
    fig: object | None
    table: pd.DataFrame | None = None
    note: str = ""
    extra_tables: dict = field(default_factory=dict)


def title(ax, text: str) -> None:
    ax.set_title(text, loc="left", fontweight="bold", color=INK)


def footnote(fig, text: str, y: float = -0.02) -> None:
    if text:
        fig.text(0.01, y, text, fontsize=7.8, color=INK2, va="top", wrap=True)


def stacked_pct(ax, fd, labels, comps, totals=None, min_lbl=5.0) -> None:
    """100 % stacked horizontal bars of the diagnostic groups (thick bars)."""
    y = np.arange(len(labels))[::-1]
    left = np.zeros(len(labels))
    for g in fd.order:
        v = np.array([c[g] for c in comps], dtype=float)
        ax.barh(y, v, left=left, height=0.92, color=fd.colors[g], edgecolor="white", linewidth=1.0,
                label=g)
        for yi, l, w in zip(y, left, v):
            if w >= min_lbl:
                ax.text(l + w / 2, yi, f"{w:.0f}%", ha="center", va="center", fontsize=8,
                        color="white" if g in fd.grouper.light_text else INK)
        left += v
    ax.set_yticks(y, labels, fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_ylim(y.min() - 0.5, y.max() + 0.5)
    ax.set_xlabel("% of subjects")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    if totals is not None:
        for yi, t in zip(y, totals):
            ax.text(101, yi, f"n={t:,}", va="center", fontsize=8.5, color=INK2)


def bottom_legend(fig, handles, ncol=None):
    """Legend centred just below the figure (kept by bbox_inches='tight')."""
    fig.legend(handles=handles, ncol=ncol or min(len(handles), 7), loc="upper center",
               bbox_to_anchor=(0.5, 0.0), frameon=False, fontsize=8.5)


def dx_legend(ax, fd, ncol=None, y=None, **kw):
    from matplotlib.patches import Patch
    bottom_legend(ax.figure, [Patch(color=fd.colors[g], label=g) for g in fd.order], ncol)
