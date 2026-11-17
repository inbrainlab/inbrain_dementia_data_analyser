"""Cross-sectional figures: which modalities each subject has, composite
cohorts and the diagnostic mix (one row per subject of the cohort)."""
from __future__ import annotations

import itertools
import math

import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyBboxPatch, Patch, PathPatch, Rectangle, Wedge
from matplotlib.path import Path

from .data import FigureData
from .style import (EDGE, GRID, INK, INK2, Panel, blue_cmap, bottom_legend, dx_legend, footnote, greys, plt,
                    stacked_pct, title)


def _removed_note(fd: FigureData) -> str:
    if not fd.n_removed:
        return ""
    top = ", ".join(f"{k} {v}" for k, v in fd.removed_labels.most_common(6))
    return (f"{fd.n_removed:,} subjects whose diagnosis has no group were left out ({top}"
            f"{', ...' if len(fd.removed_labels) > 6 else ''}).")


def _ref_text(fd: FigureData) -> str:
    return f"% of {fd.ref_label}" if fd.reference else "% of cohort"


# ---------------------------------------------------------------------------
# 1. modality connectome
# ---------------------------------------------------------------------------
def modality_connectome_graph(fd: FigureData) -> Panel:
    mods = fd.modalities()
    if len(mods) < 2:
        return Panel("modality_connectome_graph", None, note="needs at least two modalities")
    N = fd.N
    S = fd.sets
    center = fd.reference if fd.reference in mods else None
    ring = [m for m in mods if m != center]
    ang = np.deg2rad(90 - np.arange(len(ring)) * 360 / len(ring))
    pos = {m: (math.cos(a), math.sin(a)) for m, a in zip(ring, ang)}
    if center:
        pos[center] = (0.0, 0.0)
    cmap = blue_cmap()
    fig, ax = plt.subplots(figsize=(8.6, 8.4))
    edges = [(a, b, len(S[a] & S[b])) for a, b in itertools.combinations(mods, 2)]
    for a, b, w in sorted(edges, key=lambda e: e[2]):
        (x1, y1), (x2, y2) = pos[a], pos[b]
        rad = 0 if center in (a, b) else 0.28
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        cx, cy = mx - rad * (y2 - y1), my + rad * (x2 - x1)
        p = Path([(x1, y1), (cx, cy), (x2, y2)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
        ax.add_patch(PathPatch(p, fc="none", ec=cmap(w / N), lw=0.8 + 13 * w / N, alpha=0.85,
                               capstyle="round", zorder=1))
    gap, rw = 0.015, 0.06
    for k, (x, y) in pos.items():
        n = len(S[k])
        r = 0.08 + 0.17 * math.sqrt(n / N)
        ax.add_patch(Circle((x, y), r, fc="#fcfcfb", ec=INK if k == center else EDGE,
                            lw=1.8 if k == center else 1.2, zorder=3))
        comp, th = fd.comp_pct(S[k]), 90.0
        for g in fd.order:
            sweep = 3.6 * comp[g]
            if sweep > 0:
                ax.add_patch(Wedge((x, y), r + gap + rw, th - sweep, th, width=rw, fc=fd.colors[g],
                                   ec="white", lw=0.8, zorder=4))
                th -= sweep
        ax.text(x, y + (0.07 if k == center else 0.012), k, ha="center", va="center", fontsize=10.5,
                fontweight="bold", color=INK, zorder=5)
        if k == center:
            ax.text(x, y - 0.1, f"n={n:,}", ha="center", va="center", fontsize=8.5, color=INK2, zorder=5)
        else:
            off = r + gap + rw + 0.07
            ux, uy = x / (math.hypot(x, y) or 1), y / (math.hypot(x, y) or 1)
            ax.text(x + ux * (off + 0.13), y + uy * (off + 0.11), f"n={n:,}\n{100 * n / N:.1f}%",
                    ha="center", va="center", fontsize=8.5, color=INK2, zorder=5)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, N))
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.035, pad=0.0, shrink=0.55)
    cb.set_label("shared subjects (edge colour & width)", fontsize=8.5, color=INK2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8)
    ax.legend(handles=[Patch(color=fd.colors[g], label=g) for g in fd.order], title="Ring = diagnosis mix",
              loc="upper left", bbox_to_anchor=(-0.02, 1.0), frameon=False, fontsize=8,
              title_fontproperties={"weight": "bold", "size": 8.5})
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.68); ax.set_aspect("equal"); ax.axis("off")
    title(ax, "Modality connectome\nnode size = subjects with the modality · edge = shared subjects")
    footnote(fig, f"n = {N:,} subjects ({_ref_text(fd)} next to each node). " + _removed_note(fd))
    tab = pd.DataFrame(edges, columns=["modality_1", "modality_2", "shared_N"])
    comp = pd.DataFrame({m: fd.comp_pct(S[m]).round(2) for m in mods}).T
    comp.insert(0, "N", [len(S[m]) for m in mods])
    return Panel("modality_connectome_graph", fig, tab, extra_tables={"modality_dx_pct": comp})


# ---------------------------------------------------------------------------
# 2. overlap matrix
# ---------------------------------------------------------------------------
def modality_overlap_matrix(fd: FigureData) -> Panel:
    mods = fd.modalities()
    if len(mods) < 2:
        return Panel("modality_overlap_matrix", None, note="needs at least two modalities")
    N, S = fd.N, fd.sets
    M = np.array([[len(S[a] & S[b]) for b in mods] for a in mods])
    fig, ax = plt.subplots(figsize=(1.25 * len(mods) + 2.5, 1.15 * len(mods) + 1.8))
    ax.imshow(np.ma.masked_array(M, np.triu(np.ones_like(M, bool), 1)), cmap=blue_cmap(),
              vmin=0, vmax=N, aspect="equal")
    for i in range(len(mods)):
        for j in range(i + 1):
            ax.text(j, i, f"{M[i, j]:,}\n{100 * M[i, j] / N:.1f}%", ha="center", va="center", fontsize=8.5,
                    color="white" if M[i, j] / N > 0.45 else INK, fontweight="bold" if i == j else None)
    ax.set_xticks(range(len(mods)), mods); ax.set_yticks(range(len(mods)), mods)
    ax.set_xticks(np.arange(-.5, len(mods)), minor=True); ax.set_yticks(np.arange(-.5, len(mods)), minor=True)
    ax.grid(which="minor", color="white", lw=2.5); ax.tick_params(which="both", length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    title(ax, f"Shared subjects between modalities\n(diagonal = modality n; {_ref_text(fd)})")
    return Panel("modality_overlap_matrix", fig, pd.DataFrame(M, index=mods, columns=mods))


# ---------------------------------------------------------------------------
# composite cohorts
# ---------------------------------------------------------------------------
def composites(fd: FigureData) -> list[tuple[tuple[str, ...], set]]:
    """(step names, subjects) for each composite cohort, smallest first."""
    steps = dict(fd.steps())
    names = list(steps)
    spec = fd.settings.get("composites", "auto")
    if isinstance(spec, list):
        combos = [()]
        for c in spec:
            c = tuple(x for x in names if x in c)          # keep step order
            if c and c not in combos:
                combos.append(c)
    elif len(names) <= 4:
        combos = [c for k in range(len(names) + 1) for c in itertools.combinations(names, k)]
    else:
        # greedy lattice: at each level add every remaining step to the largest
        # composite of the previous level (k + (k-1) + ... + 1 boxes)
        combos, cur = [()], ()
        while len(cur) < len(names):
            nxt = [tuple(x for x in names if x in set(cur) | {y}) for y in names if y not in cur]
            combos += nxt
            size = lambda c: len(set(fd.U).intersection(*(steps[x] for x in c)))
            cur = max(nxt, key=size)
    combos = list(dict.fromkeys(combos))
    out = []
    for c in combos:
        s = set(fd.U)
        for x in c:
            s &= steps[x]
        out.append((c, s))
    return out


def _cname(fd: FigureData, c: tuple, sep=" + ") -> str:
    return sep.join([fd.ref_label] + list(c))


def composite_cohort_lattice(fd: FigureData) -> Panel:
    C = composites(fd)
    if len(C) < 2:
        return Panel("composite_cohort_lattice", None, note="needs at least one modality besides the reference")
    N = fd.N
    cols: dict[int, list] = {}
    for c, s in C:
        cols.setdefault(len(c), []).append((c, s))
    BW, BH, gap_y = 0.86, 1.75, 0.5
    pos = {}
    for k, items in cols.items():
        items.sort(key=lambda t: -len(t[1]))
        h = len(items) * BH + (len(items) - 1) * gap_y
        for i, (c, _) in enumerate(items):
            pos[c] = (k, h / 2 - BH / 2 - i * (BH + gap_y))
    rows = max(len(v) for v in cols.values())
    fig, ax = plt.subplots(figsize=(3.3 * len(cols) + 1.5, 1.05 * rows * (BH + gap_y) + 1.5))
    sets = dict(C)
    for a in sets:
        for b in sets:
            if len(b) == len(a) + 1 and set(a) <= set(b):
                (x1, y1), (x2, y2) = pos[a], pos[b]
                xs, xe = x1 + BW / 2, x2 - BW / 2
                xm = (xs + xe) / 2
                p = Path([(xs, y1), (xm, y1), (xm, y2), (xe, y2)],
                         [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                ax.add_patch(PathPatch(p, fc="none", ec="#2a78d6", lw=1 + 16 * len(sets[b]) / N,
                                       alpha=0.35, capstyle="butt", zorder=1))
    for c, (x, y) in pos.items():
        s = sets[c]
        ax.add_patch(FancyBboxPatch((x - BW / 2, y - BH / 2), BW, BH, boxstyle="round,pad=0,rounding_size=0.06",
                                    fc="#fcfcfb", ec=INK if not c else EDGE, lw=1.4 if not c else 1, zorder=2))
        ax.text(x, y + BH / 2 - 0.09, _cname(fd, c, "\n+ "), ha="center", va="top", fontsize=8.2, color=INK,
                zorder=3, linespacing=1.15)
        ax.text(x, y - BH / 2 + 0.56, f"n={len(s):,}", ha="center", va="center", fontsize=9.5,
                fontweight="bold", color=INK, zorder=3)
        ax.text(x, y - BH / 2 + 0.33, f"{100 * len(s) / N:.1f}% of {fd.ref_label}", ha="center",
                va="center", fontsize=7.8, color=INK2, zorder=3)
        comp, x0, wt = fd.comp_pct(s), x - BW / 2 + 0.05, BW - 0.1
        for g in fd.order:
            w = wt * comp[g] / 100
            if w > 0:
                ax.add_patch(Rectangle((x0, y - BH / 2 + 0.1), w, 0.16, fc=fd.colors[g], ec="none", zorder=3))
            x0 += w
    ys = [y for _, y in pos.values()]
    ax.set_xlim(-0.5, max(cols) + 0.5)
    ax.set_ylim(min(ys) - BH / 2 - 0.3, max(ys) + BH / 2 + 0.3)
    ax.axis("off")
    title(ax, "Composite cohorts - each column adds a modality\n(link width ∝ n; coloured strip = diagnosis mix)")
    dx_legend(ax, fd, y=0.0)
    tab = pd.DataFrame({"composite": [_cname(fd, c) for c, _ in C], "n_modalities_added": [len(c) for c, _ in C],
                        "N": [len(s) for _, s in C], "pct_of_reference": [round(100 * len(s) / N, 2) for _, s in C]})
    return Panel("composite_cohort_lattice", fig, tab)


def composite_cohort_dx_stack(fd: FigureData) -> Panel:
    C = composites(fd)
    labels = [_cname(fd, c) for c, _ in C]
    comps = [fd.comp_pct(s) for _, s in C]
    fig, ax = plt.subplots(figsize=(11, 0.55 * len(C) + 1.6))
    stacked_pct(ax, fd, labels, comps, totals=[len(s) for _, s in C])
    title(ax, "Diagnosis in each composite cohort")
    dx_legend(ax, fd, y=-0.6 / max(1, len(C) ** 0.5))
    tab = pd.DataFrame([c.round(2) for c in comps], index=labels)
    tab.insert(0, "N", [len(s) for _, s in C])
    return Panel("composite_cohort_dx_stack", fig, tab)


# ---------------------------------------------------------------------------
# diagnosis
# ---------------------------------------------------------------------------
def dx_group_counts_log(fd: FigureData) -> Panel:
    N = fd.N
    cnt = fd.subjects["dx"].value_counts().reindex(fd.order).fillna(0).astype(int)
    x = np.arange(len(fd.order))
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.bar(x, np.maximum(cnt.values, 0.8), color=[fd.colors[g] for g in fd.order], width=0.92,
           edgecolor="white", linewidth=1)
    for xi, n in zip(x, cnt.values):
        ax.text(xi, max(n, 1) * 1.12, f"{n:,}\n{100 * n / N:.1f}%", ha="center", va="bottom", fontsize=8.5,
                color=INK)
    ax.set_yscale("log")
    ax.set_ylim(0.8, max(cnt.max(), 1) * 6)
    ax.set_xticks(x, [g.replace(" / ", "/\n") for g in fd.order], fontsize=9)
    ax.set_ylabel("Subjects (log scale)")
    ax.yaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    title(ax, f"Diagnostic groups ({fd.ref_label} cohort, n = {N:,})")
    footnote(fig, _removed_note(fd))
    return Panel("dx_group_counts_log", fig,
                 pd.DataFrame({"N": cnt, "pct": (100 * cnt / max(N, 1)).round(2)}))


def dx_subgroup_labels(fd: FigureData, top_mixed: int = 6) -> Panel:
    sub = fd.subjects
    multi = [g for g in fd.order if sub.loc[sub["dx"] == g, "dx_label"].nunique() > 1]
    if not multi:
        return Panel("dx_subgroup_labels", None, note="no group holds more than one label")
    rows = []
    for g in multi:
        d = sub[sub["dx"] == g]["dx_label"].value_counts()
        if g == fd.grouper.compound and len(d) > top_mixed + 1:
            rows += [(g, k, int(v)) for k, v in d.head(top_mixed).items()]
            rows.append((g, f"other {g.lower()}", int(d.iloc[top_mixed:].sum())))
        else:
            rows += [(g, k, int(v)) for k, v in d.items()]
    yy = np.arange(len(rows))[::-1]
    fig, ax = plt.subplots(figsize=(7.5, 0.32 * len(rows) + 1.4))
    ax.barh(yy, [r[2] for r in rows], color=[fd.colors[r[0]] for r in rows], height=0.92, edgecolor="white")
    mx = max(r[2] for r in rows)
    for yi, r in zip(yy, rows):
        ax.text(r[2] + mx * 0.01, yi, f"{r[2]:,}", va="center", fontsize=8, color=INK)
    ax.set_yticks(yy, [r[1] for r in rows], fontsize=8.5)
    ax.set_xlabel("Subjects")
    ax.set_ylim(-0.5, len(rows) - 0.5); ax.set_xlim(0, mx * 1.15)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    title(ax, "Labels within " + ", ".join(multi))
    ax.legend(handles=[Patch(color=fd.colors[g], label=g) for g in multi], loc="lower right",
              frameon=False, fontsize=8.5)
    return Panel("dx_subgroup_labels", fig, pd.DataFrame(rows, columns=["group", "label", "N"]))


def dx_by_databank_heatmap(fd: FigureData) -> Panel:
    ct = pd.crosstab(fd.subjects["Databank"], fd.subjects["dx"]).reindex(columns=fd.order).fillna(0).astype(int)
    if ct.empty:
        return Panel("dx_by_databank_heatmap", None, note="no subjects")
    ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
    ct.loc["Total"] = ct.sum()
    vals = np.log10(ct.values + 1)
    fig, ax = plt.subplots(figsize=(0.9 * len(fd.order) + 3, 0.62 * len(ct) + 2))
    ax.imshow(vals, cmap=plt.matplotlib.colors.LinearSegmentedColormap.from_list(
        "b", ["#f4f7fc", "#2a78d6", "#0d3a75"]), aspect="equal")
    hi = vals.max() * 0.6 if vals.max() else 1
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            v = ct.values[i, j]
            ax.text(j, i, f"{v:,}" if v else "–", ha="center", va="center", fontsize=8.5,
                    color="white" if vals[i, j] > hi else INK,
                    fontweight="bold" if i == ct.shape[0] - 1 else None)
    ax.set_xticks(range(len(fd.order)), fd.order, fontsize=9, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(ct.shape[0]), ct.index, fontsize=9)
    ax.set_xticks(np.arange(-.5, ct.shape[1]), minor=True); ax.set_yticks(np.arange(-.5, ct.shape[0]), minor=True)
    ax.grid(which="minor", color="white", lw=2); ax.tick_params(which="both", length=0)
    ax.axhline(ct.shape[0] - 1.5, color="white", lw=5)
    for s in ax.spines.values():
        s.set_visible(False)
    title(ax, "Diagnostic group × databank\n(subjects; log shading)")
    return Panel("dx_by_databank_heatmap", fig, ct)


# ---------------------------------------------------------------------------
# Venn of two modalities
# ---------------------------------------------------------------------------
def _lens(r1, r2, d):
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    a1 = r1 ** 2 * math.acos((d ** 2 + r1 ** 2 - r2 ** 2) / (2 * d * r1))
    a2 = r2 ** 2 * math.acos((d ** 2 + r2 ** 2 - r1 ** 2) / (2 * d * r2))
    a3 = 0.5 * math.sqrt(max(0.0, (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2)))
    return a1 + a2 - a3


def _venn_distance(r1, r2, target):
    lo, hi = abs(r1 - r2), r1 + r2
    if target <= 0:
        return hi * 1.05
    if target >= math.pi * min(r1, r2) ** 2:
        return lo
    for _ in range(80):
        mid = (lo + hi) / 2
        if _lens(r1, r2, mid) > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def venn_regions(fd: FigureData):
    a, b = fd.settings["venn"]
    if a not in fd.sets or b not in fd.sets:
        return None
    A, B = fd.sets[a], fd.sets[b]
    return a, b, {f"{a} only": A - B, "Both": A & B, f"{b} only": B - A}


def _venn_name(fd: FigureData) -> str:
    a, b = fd.settings["venn"]
    slug = lambda s: "".join(ch for ch in s.lower() if ch.isalnum()) or "x"
    return f"{slug(a)}_{slug(b)}_venn"


def modality_venn(fd: FigureData) -> Panel:
    name = _venn_name(fd)
    vr = venn_regions(fd)
    if vr is None:
        return Panel(name, None, note="one of the Venn modalities is not in this cohort")
    a, b, reg = vr
    nA, nB = len(fd.sets[a]), len(fd.sets[b])
    union = len(fd.sets[a] | fd.sets[b])
    if union == 0:
        return Panel(name, None, note="no subject has either modality")
    r1, r2 = math.sqrt(max(nA, 1e-9) / math.pi), math.sqrt(max(nB, 1e-9) / math.pi)
    d = _venn_distance(r1, r2, len(reg["Both"]))
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    c1, c2 = (0.0, 0.0), (d, 0.0)
    for (cx, cy), r, col in ((c1, r1, "#2a78d6"), (c2, r2, "#eb6834")):
        ax.add_patch(Circle((cx, cy), r, fc=col, alpha=0.5, ec="white", lw=2))
    pct = lambda n: f"{n:,}\n({100 * n / union:.1f}%)"
    xa = (-r1 + min(d - r2, r1)) / 2 if nA else 0
    xb = (max(d - r2, r1) + d + r2) / 2 if nB else d
    xm = (max(d - r2, -r1) + min(r1, d + r2)) / 2
    ax.text(xa, 0, pct(len(reg[f"{a} only"])), ha="center", va="center", fontsize=10.5)
    ax.text(xb, 0, pct(len(reg[f"{b} only"])), ha="center", va="center", fontsize=10.5)
    if reg["Both"]:
        ax.text(xm, 0, pct(len(reg["Both"])), ha="center", va="center", fontsize=10.5)
    ax.text(-r1 * 0.55, r1 * 0.95, a, ha="right", va="bottom", fontsize=12, fontweight="bold")
    ax.text(d + r2 * 0.55, r2 * 0.95, b, ha="left", va="bottom", fontsize=12, fontweight="bold")
    ax.set_xlim(-r1 * 1.45, d + r2 * 1.45)
    ax.set_ylim(-max(r1, r2) * 1.25, max(r1, r2) * 1.4)
    ax.set_aspect("equal"); ax.axis("off")
    title(ax, f"{a} and/or {b} (union = {union:,})")
    tab = pd.DataFrame({"region": list(reg), "N": [len(v) for v in reg.values()],
                        "pct_of_union": [round(100 * len(v) / union, 2) for v in reg.values()]})
    return Panel(name, fig, tab)


def venn_region_by_databank(fd: FigureData) -> Panel:
    vr = venn_regions(fd)
    if vr is None:
        return Panel("venn_region_by_databank", None, note="Venn modalities not in this cohort")
    _, _, reg = vr
    rk = list(reg)
    region = pd.Series({u: k for k, s in reg.items() for u in s}, name="region")
    if region.empty:
        return Panel("venn_region_by_databank", None, note="no subject has either modality")
    bank = fd.subjects.loc[region.index, "Databank"]
    dbs = bank.value_counts().index.tolist()
    col = dict(zip(dbs, greys(len(dbs))))
    ctab = pd.crosstab(region, bank).reindex(index=rk, columns=dbs).fillna(0)
    pct = ctab.div(ctab.sum(axis=1).replace(0, np.nan), axis=0).fillna(0) * 100
    fig, ax = plt.subplots(figsize=(10, 0.7 * len(rk) + 1.8))
    yy, left = np.arange(len(rk))[::-1], np.zeros(len(rk))
    for i, db in enumerate(dbs):
        w = pct[db].values
        ax.barh(yy, w, left=left, color=col[db], height=0.92, edgecolor="white", linewidth=1, label=db)
        for yi, l, wi in zip(yy, left, w):
            if wi >= 5:
                ax.text(l + wi / 2, yi, f"{wi:.0f}%", ha="center", va="center", fontsize=8,
                        color="white" if i < len(dbs) / 2 else INK)
        left += w
    for yi, t in zip(yy, ctab.sum(axis=1).values):
        ax.text(101, yi, f"n={int(t):,}", va="center", fontsize=8.5, color=INK2)
    ax.set_yticks(yy, rk); ax.set_xlim(0, 100); ax.set_ylim(-0.5, len(rk) - 0.5)
    ax.set_xlabel("% of subjects")
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    bottom_legend(fig, [Patch(color=col[db], label=db) for db in dbs], min(len(dbs), 6))
    title(ax, "Source databank per Venn region")
    return Panel("venn_region_by_databank", fig, ctab.astype(int))


def venn_region_by_dx(fd: FigureData) -> Panel:
    vr = venn_regions(fd)
    if vr is None:
        return Panel("venn_region_by_dx", None, note="Venn modalities not in this cohort")
    _, _, reg = vr
    rk = list(reg)
    fig, ax = plt.subplots(figsize=(10, 0.7 * len(rk) + 1.8))
    stacked_pct(ax, fd, rk, [fd.comp_pct(reg[k]) for k in rk], totals=[len(reg[k]) for k in rk])
    title(ax, "Diagnosis per Venn region")
    dx_legend(ax, fd, y=-0.32)
    tab = pd.DataFrame({k: fd.subjects.loc[list(reg[k]), "dx"].value_counts().reindex(fd.order).fillna(0)
                        for k in rk}).T.astype(int)
    return Panel("venn_region_by_dx", fig, tab)


CROSS_SECTIONAL = {
    "modality_connectome_graph": modality_connectome_graph,
    "modality_overlap_matrix": modality_overlap_matrix,
    "composite_cohort_lattice": composite_cohort_lattice,
    "composite_cohort_dx_stack": composite_cohort_dx_stack,
    "dx_group_counts_log": dx_group_counts_log,
    "dx_subgroup_labels": dx_subgroup_labels,
    "dx_by_databank_heatmap": dx_by_databank_heatmap,
    "modality_venn": modality_venn,
    "venn_region_by_databank": venn_region_by_databank,
    "venn_region_by_dx": venn_region_by_dx,
}
