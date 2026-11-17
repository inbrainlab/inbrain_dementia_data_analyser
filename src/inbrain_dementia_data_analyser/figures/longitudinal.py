"""Longitudinal figures: diagnostic trajectories, conversions and follow-up
(diagnosis history of each cohort subject in the expanded master)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Patch

from .data import LOST, NO_VISIT, OTHER, FigureData
from .style import (GRID, INK, INK2, LOST_COLOR, NO_VISIT_COLOR, OFF, OTHER_COLOR, TRANSITION_COLORS,
                    Panel, footnote, plt, title)


def _iqr(s: pd.Series) -> tuple[float, float, float]:
    s = s.dropna()
    if s.empty:
        return (np.nan, np.nan, np.nan)
    return (s.median(), s.quantile(.25), s.quantile(.75))


def _upset(fd: FigureData, bars, bar_colors, top_text, set_sizes, set_title, ylab, matrix_fn, legend,
           heading, note, hatch_vals=None):
    rows = fd.order
    nC = len(bars)
    fig = plt.figure(figsize=(max(9, 0.62 * nC + 4.5), 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], width_ratios=[1, max(3.2, 0.62 * nC / 1.6)],
                          hspace=0.04, wspace=0.42)
    axT = fig.add_subplot(gs[0, 1])
    axM = fig.add_subplot(gs[1, 1], sharex=axT)
    axS = fig.add_subplot(gs[1, 0], sharey=axM)
    x = np.arange(nC)
    axT.bar(x, bars, width=0.72, color=bar_colors)
    if hatch_vals is not None:
        axT.bar(x, hatch_vals, width=0.72, color="none", edgecolor="white", hatch="////", linewidth=0)
    ymax = max(max(bars), 1)
    for xi, b, t in zip(x, bars, top_text):
        axT.text(xi, b + ymax * 0.02, f"{int(b):,}", ha="center", va="bottom", fontsize=8.5, color=INK)
        axT.text(xi, ymax * 1.42, t, ha="center", va="top", fontsize=8.3, color=INK, rotation=90)
    axT.axhline(ymax * 1.17, color=GRID, lw=0.8)
    axT.text(-0.6, ymax * 1.45, ylab, ha="left", va="bottom", fontsize=9, color=INK2, fontweight="bold")
    axT.set_ylim(0, ymax * 1.45); axT.set_ylabel("Subjects")
    axT.set_yticks([t for t in axT.get_yticks() if 0 <= t <= ymax * 1.05])
    axT.tick_params(axis="x", bottom=False, labelbottom=False); axT.spines["bottom"].set_visible(False)
    axT.yaxis.grid(True, color=GRID, lw=0.8); axT.set_axisbelow(True)
    axT.set_title(heading, loc="left", fontweight="bold", color=INK, pad=30)
    axL = fig.add_subplot(gs[0, 0]); axL.axis("off")
    if legend:
        axL.legend(handles=legend, loc="lower left", bbox_to_anchor=(0.0, 0.02), frameon=False, fontsize=8.8,
                   title="Legend", title_fontproperties={"weight": "bold", "size": 9})
    yr = {r: len(rows) - 1 - i for i, r in enumerate(rows)}
    for i, r in enumerate(rows):
        if i % 2 == 0:
            axM.axhspan(yr[r] - 0.5, yr[r] + 0.5, color="#f4f4f1", zorder=0)
    for xi in x:
        axM.scatter([xi] * len(rows), list(yr.values()), s=60, color=OFF, zorder=1)
        matrix_fn(axM, xi, yr)
    axM.set_ylim(-0.5, len(rows) - 0.5); axM.set_xlim(-0.6, nC - 0.4)
    axM.tick_params(axis="x", bottom=False, labelbottom=False); axM.tick_params(axis="y", left=False, labelleft=False)
    for s in axM.spines.values():
        s.set_visible(False)
    yy = [yr[r] for r in rows]
    axS.barh(yy, [set_sizes.get(r, 0) for r in rows], height=0.6, color=[fd.colors[r] for r in rows])
    smax = max(max(set_sizes.values(), default=0), 1)
    for r in rows:
        v = set_sizes.get(r, 0)
        axS.text(v + smax * 0.03, yr[r], f"{v:,}", va="center", ha="right", fontsize=8, color=INK)
    axS.set_xlim(smax * 1.45, 0)
    axS.yaxis.tick_right()
    axS.tick_params(axis="y", length=0, labelright=True, labelleft=False, labelsize=9.5, pad=6)
    axS.set_yticks(yy, rows)
    for sp in ("left", "right", "top"):
        axS.spines[sp].set_visible(False)
    axS.set_xlabel(set_title, fontsize=8.5); axS.tick_params(axis="x", labelsize=7.5)
    footnote(fig, note)
    return fig


def converters_transition_upset(fd: FigureData) -> Panel:
    tr = fd.traj
    name = "converters_transition_upset"
    if tr.empty:
        return Panel(name, None, note="no conversion information (expanded master or 1st_diag/converted)")
    min_n = int(fd.settings["longitudinal"]["upset_min_n"])
    conv = tr[(tr["converted"] == 1) & tr["baseline_group"].notna() & tr["conversion_to_group"].notna()].copy()
    n_unknown = int(((tr["converted"] == 1) & tr["baseline_group"].isna()).sum())
    if conv.empty:
        return Panel(name, None, note="no converters in this cohort")
    conv["type"] = [fd.grouper.transition(a, b) for a, b in zip(conv["baseline_group"], conv["conversion_to_group"])]
    cg = (conv.groupby(["baseline_group", "conversion_to_group", "type"])["months_to_conversion"]
          .agg(n="size", median="median", q1=lambda s: s.quantile(.25), q3=lambda s: s.quantile(.75))
          .reset_index().rename(columns={"baseline_group": "from", "conversion_to_group": "to"})
          .sort_values("n", ascending=False))
    plot = cg[cg["n"] >= min_n].reset_index(drop=True)
    n_small = int(cg.loc[cg["n"] < min_n, "n"].sum())
    if plot.empty:
        return Panel(name, None, table=cg, note=f"no transition with n >= {min_n}")

    def mfn(ax, xi, yr):
        r = plot.iloc[xi]
        y0, y1 = yr[r["from"]], yr[r["to"]]
        ax.add_patch(FancyArrowPatch((xi, y0), (xi, y1), arrowstyle="-|>", mutation_scale=13, lw=2.2,
                                     color=TRANSITION_COLORS[r["type"]], shrinkA=5, shrinkB=6, zorder=2))
        ax.scatter([xi], [y0], s=90, color=INK, zorder=3)
        ax.scatter([xi], [y1], s=90, facecolor="white", edgecolor=INK, lw=1.8, zorder=3)

    sets = {g: int(((conv["baseline_group"] == g) | (conv["conversion_to_group"] == g)).sum()) for g in fd.order}
    txt = [(f"{r['median']:.1f} ({r.q1:.0f}–{r.q3:.0f})" if pd.notna(r["median"]) else "time n/a")
           for _, r in plot.iterrows()]
    legend = ([Patch(color=v, label=k) for k, v in TRANSITION_COLORS.items()] +
              [plt.Line2D([], [], marker="o", ls="", color=INK, label="baseline diagnosis"),
               plt.Line2D([], [], marker="o", ls="", mfc="white", mec=INK, label="converted-to diagnosis")])
    note = (f"Baseline = first diagnosis; destination = first later diagnosis in another group; time = months "
            f"from the first diagnosis (source: {fd._traj_source}). Columns with n < {min_n} not shown "
            f"({n_small} subjects). {n_unknown} converters without a group at baseline excluded.\n"
            "Direction: CN < MCI < dementia (rank of diagnosis_groups.stage_rank); dementia -> dementia = change "
            "of aetiology.")
    fig = _upset(fd, plot["n"].values, [TRANSITION_COLORS[t] for t in plot["type"]], txt, sets,
                 "Converters involving\nthis diagnosis", "Median months to conversion (IQR)", mfn, legend,
                 f"Converters in the {fd.ref_label} cohort (n = {len(conv):,}): diagnostic transitions", note)
    return Panel(name, fig, cg.round(1))


def non_converters_stability_upset(fd: FigureData) -> Panel:
    tr = fd.traj
    name = "non_converters_stability_upset"
    if tr.empty:
        return Panel(name, None, note="no conversion information")
    nc = tr[(tr["converted"] == 0) & tr["baseline_group"].notna()].copy()
    n_unknown = int(((tr["converted"] == 0) & tr["baseline_group"].isna()).sum())
    if nc.empty:
        return Panel(name, None, note="no non-converters")
    ng = []
    for g in fd.order:
        d = nc[nc["baseline_group"] == g]
        if not len(d):
            continue
        f = d["followup_months"][d["followup_months"] > 0]
        med, q1, q3 = _iqr(f)
        ng.append(dict(dx=g, n=len(d), n_single=int((d["followup_months"].fillna(0) <= 0).sum()),
                       median=med, q1=q1, q3=q3, n_fu=len(f)))
    ng = pd.DataFrame(ng).sort_values("n", ascending=False).reset_index(drop=True)

    def mfn(ax, xi, yr):
        ax.scatter([xi], [yr[ng.dx[xi]]], s=90, color=INK, zorder=3)

    txt = [(f"{r['median']:.1f} ({r.q1:.0f}–{r.q3:.0f})" if r.n_fu else "no follow-up") for _, r in ng.iterrows()]
    note = ("Diagnosis = first diagnosis (no change of group recorded). Follow-up = last - first diagnosis "
            "(months); median / IQR among subjects with > 0 months. Hatched = single assessment.\n"
            f"Source: {fd._traj_source}. "
            f"{n_unknown:,} non-converters without a group at baseline excluded.")
    fig = _upset(fd, ng["n"].values, [fd.colors[g] for g in ng["dx"]], txt, dict(zip(ng["dx"], ng["n"])),
                 "Non-converters", "Median months of follow-up without conversion (IQR)", mfn,
                 [Patch(facecolor="#8a8984", hatch="////", edgecolor="white", label="single assessment (no follow-up)")],
                 f"Non-converters in the {fd.ref_label} cohort (n = {len(nc):,}): stable diagnosis", note,
                 hatch_vals=ng["n_single"].values)
    return Panel(name, fig, ng.round(1))


# ---------------------------------------------------------------------------
# Sankey of diagnostic states per time bin
# ---------------------------------------------------------------------------
def _sigmoid_band(ax, x0, x1, a0, a1, b0, b1, color, alpha=0.45):
    """Band from [a0, a1] at x0 to [b0, b1] at x1 (smooth)."""
    t = np.linspace(0, 1, 40)
    s = t * t * (3 - 2 * t)
    xs = x0 + (x1 - x0) * t
    lo = a0 + (b0 - a0) * s
    hi = a1 + (b1 - a1) * s
    ax.fill_between(xs, lo, hi, color=color, alpha=alpha, lw=0, zorder=1)


def dx_transition_sankey(fd: FigureData) -> Panel:
    W = int(fd.settings["longitudinal"]["window_months"])
    name = f"dx_transition_sankey_{W}m"
    if fd.history is None:
        return Panel(name, None, note="needs the expanded master (diagnosis per visit)")
    st = fd.states_by_bin()
    if st.empty:
        return Panel(name, None, note="no diagnosis history")
    min_flow = int(fd.settings["longitudinal"]["sankey_min_flow"])
    active = [b for b in st.columns if (st[b] != LOST).any()]
    st = st[[b for b in st.columns if b <= max(active)]] if active else st
    bins = list(st.columns)
    states = fd.order + [OTHER, NO_VISIT, LOST]
    color = dict(fd.colors, **{OTHER: OTHER_COLOR, NO_VISIT: NO_VISIT_COLOR, LOST: LOST_COLOR})
    N = len(st)
    gap = N * 0.012
    node_w = 0.18
    counts = {b: st[b].value_counts().reindex(states).fillna(0).astype(int) for b in bins}
    # node positions (top -> bottom in the order of `states`)
    ypos = {}
    top = N + gap * len(states)
    for b in bins:
        n_nodes = int((counts[b] > 0).sum())
        y = top - (top - (N + gap * max(n_nodes - 1, 0))) / 2          # centre each column
        for s in states:
            h = counts[b][s]
            ypos[(b, s)] = (y - h, y)
            y -= h + (gap if h else 0)
    fig, ax = plt.subplots(figsize=(1.35 * len(bins) + 3, 7.5))
    flows = []
    for i in range(len(bins) - 1):
        b0, b1 = bins[i], bins[i + 1]
        ft = pd.crosstab(st[b0], st[b1]).reindex(index=states, columns=states).fillna(0).astype(int)
        out_off = {s: ypos[(b0, s)][1] for s in states}
        in_off = {s: ypos[(b1, s)][1] for s in states}
        for s0 in states:
            for s1 in states:
                n = ft.loc[s0, s1]
                if n <= 0:
                    continue
                a1, a0 = out_off[s0], out_off[s0] - n
                c1, c0 = in_off[s1], in_off[s1] - n
                out_off[s0] -= n
                in_off[s1] -= n
                flows.append((b0, b1, s0, s1, int(n)))
                if n >= min_flow and not (s0 == LOST and s1 == LOST):
                    _sigmoid_band(ax, i + node_w / 2, i + 1 - node_w / 2, a0, a1, c0, c1,
                                  color[s0] if s0 not in (NO_VISIT, LOST) else color[s1],
                                  alpha=0.25 if s0 in (NO_VISIT, LOST) or s1 in (NO_VISIT, LOST) else 0.5)
    for i, b in enumerate(bins):
        for s in states:
            lo, hi = ypos[(b, s)]
            if hi - lo <= 0:
                continue
            ax.add_patch(plt.Rectangle((i - node_w / 2, lo), node_w, hi - lo, fc=color[s], ec="white",
                                       lw=0.6, zorder=3))
            if hi - lo >= N * 0.035:
                ax.text(i, (lo + hi) / 2, f"{int(round(hi - lo)):,}", ha="center", va="center", fontsize=7, zorder=4,
                        color="white" if s in fd.grouper.light_text else INK, rotation=90)
        ax.text(i, N + gap * len(states) + N * 0.02, f"{b} m", ha="center", va="bottom", fontsize=9, color=INK2)
    ax.set_xlim(-0.6, len(bins) - 0.4)
    ax.set_ylim(-N * 0.02, N * 1.08 + gap * len(states))
    ax.axis("off")
    used = [s for s in states if any(counts[b][s] for b in bins)]
    ax.legend(handles=[Patch(color=color[s], label=s) for s in used], ncol=min(len(used), 5),
              loc="upper center", bbox_to_anchor=(0.5, 0.0), frameon=False, fontsize=8.5)
    title(ax, f"Diagnostic trajectory in {W}-month windows since the first diagnosis (n = {N:,})")
    n_out = fd.N - N
    footnote(fig, (f"{n_out:,} cohort subjects have no diagnosis history in the expanded master and are not "
                   "shown. " if n_out > 0 else "") + f"State in each window = last diagnosis of the window; '{NO_VISIT}' = no diagnosis in the "
                  f"window but one later; '{LOST}' = no later diagnosis. '{OTHER}' = label without a group "
                  f"(e.g. Dementia-NOS). Flows with fewer than {min_flow} subjects are not drawn.", y=-0.05)
    tab = pd.DataFrame(flows, columns=["from_month", "to_month", "from_state", "to_state", "N"])
    return Panel(name, fig, tab, extra_tables={"sankey_states": st})


def followup_availability_area_by_dx(fd: FigureData) -> Panel:
    name = "followup_availability_area_by_dx"
    if fd.history is None:
        return Panel(name, None, note="needs the expanded master")
    labels = [m for m in fd.settings["longitudinal"]["availability"] if m in fd.settings["modalities"]]
    tabs = {m: t for m in labels if (t := fd.availability(m)) is not None and t.values.sum() > 0}
    if not tabs:
        return Panel(name, None, note="none of the availability modalities has acquisitions")
    W = int(fd.settings["longitudinal"]["window_months"])
    n = len(tabs)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.8 * nrow + 0.8), squeeze=False)
    st = fd.states_by_bin()
    followed = None
    if not st.empty:
        followed = (st != LOST).sum()
    long = []
    for ax, (m, t) in zip(axes.ravel(), tabs.items()):
        x = t.index.values
        ax.stackplot(x, [t[g].values for g in fd.order], colors=[fd.colors[g] for g in fd.order],
                     labels=fd.order, alpha=0.9, edgecolor="white", linewidth=0.4)
        if followed is not None:
            f = followed.reindex(x)
            ax.plot(x, f.values, color=INK2, lw=1, ls="--", label="followed (any diagnosis)")
        ax.set_xticks(x[::max(1, len(x) // 6)])
        ax.set_xlabel("months since first diagnosis (window start)")
        ax.set_ylabel("subjects")
        ax.yaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
        title(ax, f"{m}  (baseline n with {m} = {int(t.iloc[0].sum()):,})")
        long.append(t.assign(modality=m).reset_index().rename(columns={"index": "month"}))
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    h, l = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=min(len(l), 8), loc="lower center", bbox_to_anchor=(0.5, 0.025), frameon=False,
               fontsize=8.5)
    fig.suptitle(f"Subjects with each measure in {W}-month windows, by baseline diagnosis",
                 x=0.01, ha="left", fontweight="bold", color=INK)
    fig.tight_layout(rect=(0, 0.08, 1, 0.97))
    footnote(fig, "Windows of the expanded master (time since the first diagnosis, or since the first "
                  "acquisition when the subject has no diagnosis date). Dashed = subjects with a diagnosis "
                  "in that window or later.", y=0.015)
    return Panel(name, fig, pd.concat(long, ignore_index=True))


LONGITUDINAL = {
    "converters_transition_upset": converters_transition_upset,
    "non_converters_stability_upset": non_converters_stability_upset,
    "dx_transition_sankey": dx_transition_sankey,
    "followup_availability_area_by_dx": followup_availability_area_by_dx,
}
