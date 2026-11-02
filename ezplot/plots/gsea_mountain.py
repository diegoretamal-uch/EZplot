"""GSEA enrichment ("mountain") plot: running ES curve, member hit barcode, and ranked metric.

Supports both single-pathway detailed view (with NES, P, FDR q, leading-edge cutoff and optional
top leading-edge gene callouts) and multi-pathway overlay (stacked hit barcodes per pathway).
"""

from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

from .. import style
from ..gsea.engine import enrichment_curve

PALETTE = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#332288", "#882255"]


@dataclass
class GseaMountainParams:
    pathways: list = field(default_factory=list)   # pathway IDs; empty = top 1 from results
    title: str = ""                                # empty = auto (pathway name or collection title)
    subtitle: str = "Positive NES: {treatment} · Negative NES: {control}"
    show_ranking: bool = True                      # bottom panel with ranked metric distribution
    show_stats: bool = True                        # NES / p / FDR annotation box or legend
    label_genes: int = 0                           # top N leading-edge genes to label on curve
    line_width: float = 2.4
    width: float = 8.5
    height: float = 0.0                            # 0 = automatic


def select_pathways(res, gene_sets, requested=None, max_sets=6):
    """Pick pathway rows from `res` matching `requested` (or the top pathway if empty)."""
    if res is None or res.empty:
        return []
    by_name = {r["pathway"]: r for _, r in res.iterrows()}
    if requested:
        chosen = []
        for p in requested:
            if p in by_name and p in gene_sets:
                chosen.append(by_name[p])
            else:
                # fallback: match by pretty_name or case-insensitive substring
                q = str(p).strip().lower()
                match = next(
                    (
                        r
                        for name, r in by_name.items()
                        if name in gene_sets
                        and (name.lower() == q or style.pretty_name(name, None).lower() == q)
                    ),
                    None,
                )
                if match is not None:
                    chosen.append(match)
        if chosen:
            return chosen[:max_sets]
    for _, r in res.iterrows():
        if r["pathway"] in gene_sets:
            return [r]
    return []


def _place_gene_callouts(ax, ranks, es_vals, genes, n_genes, es_sign, color):
    """Annotate the top leading-edge genes along the ES curve without overlaps or crossed lines."""
    k = len(genes)
    if k == 0:
        return
    ymin, ymax = ax.get_ylim()
    yr = max(1e-6, ymax - ymin)

    # Sort hits by rank ascending; on both the left climb (ES > 0) and right climb (ES < 0),
    # es_vals increases monotonically with rank among the top leading-edge hits.
    order = np.argsort(ranks)
    ranks = np.asarray(ranks, dtype=float)[order]
    es_vals = np.asarray(es_vals, dtype=float)[order]
    genes = np.asarray(genes)[order]

    ax.scatter(ranks, es_vals, s=28, color=color, edgecolors="white", linewidths=0.85, zorder=6)

    # Monotonic vertical slots and half-plane separation guarantee zero line-crossings and zero text overlaps
    min_dy = 0.088 * yr
    gap_x = n_genes * 0.008
    if es_sign >= 0:
        y_low = max(0.07 * yr, float(es_vals.min()) - 0.01 * yr)
        y_high = min(ymax - 0.07 * yr, max(float(es_vals.max()) + 0.14 * yr, y_low + (k - 1) * min_dy))
        if y_high - y_low < (k - 1) * min_dy:
            y_low = max(ymin + 0.07 * yr, y_high - (k - 1) * min_dy)
        tys = np.linspace(y_low, y_high, k) if k > 1 else np.array([float(es_vals[0]) + 0.10 * yr])
        lx = min(n_genes * 0.58, float(ranks.max()) + n_genes * 0.068)
        line_end_x = lx - gap_x
        ha = "left"
    else:
        y_high = min(-0.07 * yr, float(es_vals.max()) + 0.01 * yr)
        y_low = max(ymin + 0.07 * yr, min(float(es_vals.min()) - 0.14 * yr, y_high - (k - 1) * min_dy))
        if y_high - y_low < (k - 1) * min_dy:
            y_high = min(ymax - 0.07 * yr, y_low + (k - 1) * min_dy)
        tys = np.linspace(y_low, y_high, k) if k > 1 else np.array([float(es_vals[0]) - 0.10 * yr])
        lx = max(n_genes * 0.42, float(ranks.min()) - n_genes * 0.068)
        line_end_x = lx + gap_x
        ha = "right"

    for rx, ey, ty, gene in zip(ranks, es_vals, tys, genes):
        ax.plot([rx, line_end_x], [ey, ty], color="#7a7f87", lw=0.85, zorder=6)
        ax.text(
            lx, ty, gene,
            fontsize=10.5,
            fontstyle="italic",
            color=style.TEXT,
            ha=ha,
            va="center",
            zorder=7,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.88),
        )



def plot_gsea_mountain(res, stats, gene_sets, collection_title, rank_label, control, treatment, p=None):
    """Render a single- or multi-pathway GSEA enrichment (mountain) plot."""
    p = p or GseaMountainParams()
    style.setup()

    rows = select_pathways(res, gene_sets, p.pathways)
    m = len(rows)
    curves = []
    for r in rows:
        c = enrichment_curve(stats, gene_sets[r["pathway"]])
        if c is not None:
            curves.append((r, c))
    m = len(curves)

    multi = m > 1
    leg_extra = 0.28 * m if (multi and p.show_stats) else 0.0
    bar_extra = 0.18 * max(0, m - 1)
    auto_h = (6.1 if p.show_ranking else 4.6) + leg_extra + bar_extra
    height = p.height or auto_h
    fig = plt.figure(figsize=(p.width, height))

    left, w_ax = 0.13, 0.82
    top_header = 1.0 - 0.24 / height
    sub_y = 1.0 - 0.58 / height

    if m == 1:
        default_title = style.pretty_name(curves[0][0]["pathway"], width=None)
    else:
        default_title = f"{collection_title} — Enrichment curves"
    title_str = (
        p.title.format(
            collection=collection_title,
            pathway= style.pretty_name(curves[0][0]["pathway"], width=None) if m else "",
        )
        if p.title
        else default_title
    )
    fig.text(left, top_header, title_str, ha="left", va="top", fontsize=17.5, fontweight="bold",
             color=style.TEXT)
    fig.text(left, sub_y, p.subtitle.format(treatment=treatment, control=control),
             ha="left", va="top", fontsize=12.0, color="#444444")

    if m == 0:
        ax = fig.add_axes([left, 0.15, w_ax, 0.68])
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_color("#cccccc")
        ax.text(0.5, 0.5, "No gene sets available to plot", transform=ax.transAxes,
                ha="center", va="center", fontsize=13, color=style.MUTED)
        return fig

    n_genes = curves[0][1]["n_genes"]
    ranked_vals = curves[0][1]["ranked_values"]
    zero_cross = curves[0][1]["zero_cross"]

    # Vertical layout in physical inches converted to figure fractions
    bot_margin = 0.68 / height
    es_top = 1.0 - (0.92 + leg_extra) / height
    bar_h = (0.44 + 0.18 * max(0, m - 1)) / height
    gap = 0.09 / height

    if p.show_ranking:
        rank_h = 1.25 / height
        rank_bot = bot_margin
        bar_bot = rank_bot + rank_h + gap
        es_bot = bar_bot + bar_h + gap
    else:
        rank_bot = rank_h = 0.0
        bar_bot = bot_margin
        es_bot = bar_bot + bar_h + gap

    ax_es = fig.add_axes([left, es_bot, w_ax, max(0.18, es_top - es_bot)])
    ax_bar = fig.add_axes([left, bar_bot, w_ax, bar_h])
    ax_rank = fig.add_axes([left, rank_bot, w_ax, rank_h]) if p.show_ranking else None

    # ── Panel 1: Running Enrichment Score (ES) ──
    ax_es.axhline(0, color="#666666", ls="--", lw=1.1, zorder=1)
    ax_es.grid(True, ls=":", color="#dcdcdc", lw=0.8, zorder=0)
    ax_es.set_axisbelow(True)

    all_y = np.concatenate([c["ys"] for _, c in curves])
    ymin, ymax = float(all_y.min()), float(all_y.max())
    span = max(0.25, ymax - ymin)
    pad_top = span * (0.24 if p.label_genes > 0 else 0.14)
    pad_bot = span * (0.24 if p.label_genes > 0 else 0.14)
    ax_es.set_ylim(min(-0.06, ymin - pad_bot), max(0.06, ymax + pad_top))
    ax_es.set_xlim(0, n_genes)

    legend_handles = []
    for idx, (row, c) in enumerate(curves):
        if m == 1:
            color = style.POS if row["NES"] >= 0 else style.NEG
        else:
            color = PALETTE[idx % len(PALETTE)]

        ax_es.plot(c["xs"], c["ys"], color=color, lw=p.line_width, zorder=4)
        if m == 1:
            ax_es.fill_between(c["xs"], 0, c["ys"], color=color, alpha=0.10, zorder=2)

        ax_es.plot([c["peak_rank"], c["peak_rank"]], [0, c["es"]], color=color, ls="--",
                   lw=1.25, alpha=0.85, zorder=3)
        ax_es.scatter([c["peak_rank"]], [c["es"]], s=38, color=color, edgecolors="white",
                      linewidths=0.9, zorder=5)

        pname = style.pretty_name(row["pathway"], width=None)
        legend_handles.append(
            Line2D(
                [], [], color=color, lw=p.line_width,
                label=f"{pname}  (NES = {row['NES']:+.2f}, FDR = {row['padj']:.2e})",
            )
        )

    if p.label_genes and p.label_genes > 0:
        row0, c0 = curves[0]
        col0 = (style.POS if row0["NES"] >= 0 else style.NEG) if m == 1 else PALETTE[0]
        k_lab = min(int(p.label_genes), len(c0["leading_genes"]))
        _place_gene_callouts(
            ax_es,
            c0["leading_ranks"][:k_lab],
            c0["leading_es"][:k_lab],
            c0["leading_genes"][:k_lab],
            n_genes,
            row0["NES"],
            col0,
        )

    if p.show_stats:
        if m == 1:
            row0, c0 = curves[0]
            n_le = len(c0["leading_genes"])
            stats_txt = (
                f"NES = {row0['NES']:+.2f}\n"
                f"$P$ = {row0['pval']:.2e}\n"
                f"FDR $q$ = {row0['padj']:.2e}\n"
                f"Leading edge: {n_le}/{int(row0['size'])}"
            )
            if row0["NES"] >= 0:
                bx, by, ha, va = 0.975, 0.94, "right", "top"
            else:
                bx, by, ha, va = 0.025, 0.06, "left", "bottom"
            ax_es.text(
                bx, by, stats_txt, transform=ax_es.transAxes, ha=ha, va=va,
                fontsize=11.5, color=style.TEXT, linespacing=1.35, zorder=8,
                bbox=dict(boxstyle="round,pad=0.45", facecolor="white", edgecolor="#cccccc",
                          linewidth=0.9, alpha=0.92),
            )
        else:
            fig.legend(
                handles=legend_handles,
                loc="upper left",
                bbox_to_anchor=(left, es_top + (leg_extra - 0.04) / height),
                frameon=False,
                fontsize=11.0,
                handlelength=1.8,
                labelspacing=0.35,
            )

    for s in ("top", "right"):
        ax_es.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax_es.spines[s].set_linewidth(1.2)
    ax_es.tick_params(labelsize=11.5, labelbottom=False)
    ax_es.set_ylabel("Enrichment score (ES)", fontsize=13, labelpad=8)

    # ── Panel 2: Hit Barcode ──
    ax_bar.set_xlim(0, n_genes)
    ax_bar.set_ylim(0, 1)
    ax_bar.set_yticks([])

    if m == 1:
        _, c0 = curves[0]
        ax_bar.vlines(c0["hits"], 0.24, 0.96, color="#222222", lw=0.75, alpha=0.75, zorder=3)
        # Diverging rank colour strip at the bottom of the barcode
        grad = np.linspace(1.0, -1.0, 256)[None, :]
        zero_frac = np.clip(zero_cross / max(1, n_genes), 0.05, 0.95)
        cmap_bar = LinearSegmentedColormap.from_list(
            "rank_bar",
            [(0.0, style.NEG), (1.0 - zero_frac, "#f7f7f7"), (1.0, style.POS)],
        )
        ax_bar.imshow(grad, aspect="auto", cmap=cmap_bar, extent=(0, n_genes, 0.0, 0.22), zorder=2)
    else:
        for idx, (_, c) in enumerate(curves):
            color = PALETTE[idx % len(PALETTE)]
            y_top = 1.0 - idx / m
            y_bot = 1.0 - (idx + 1) / m
            pad = 0.08 / m
            ax_bar.vlines(c["hits"], y_bot + pad, y_top - pad, color=color, lw=0.8, alpha=0.85, zorder=3)
            if idx > 0:
                ax_bar.axhline(y_top, color="#e0e0e0", lw=0.7, zorder=1)

    for s in ax_bar.spines.values():
        s.set_color("#888888")
        s.set_linewidth(1.0)

    if not p.show_ranking:
        ax_bar.tick_params(axis="x", labelsize=11.5)
        ax_bar.set_xlabel(f"Gene rank in ordered dataset (n = {n_genes:,})", fontsize=13, labelpad=8)
    else:
        ax_bar.tick_params(labelbottom=False, bottom=False)

    # ── Panel 3: Ranked Metric Distribution ──
    if ax_rank is not None:
        ax_rank.set_xlim(0, n_genes)
        step = max(1, n_genes // 1200)
        rx = np.arange(1, n_genes + 1, step)
        if rx[-1] != n_genes:
            rx = np.append(rx, n_genes)
        ry = ranked_vals[rx - 1]

        pos_mask = ry >= 0
        neg_mask = ry <= 0
        ax_rank.fill_between(rx, 0, np.where(pos_mask, ry, 0), color=style.POS, alpha=0.30, zorder=2)
        ax_rank.fill_between(rx, 0, np.where(neg_mask, ry, 0), color=style.NEG, alpha=0.30, zorder=2)
        ax_rank.plot(rx, ry, color="#444444", lw=0.9, zorder=3)
        ax_rank.axhline(0, color="#666666", lw=0.9, zorder=1)
        ax_rank.grid(True, ls=":", color="#e2e2e2", lw=0.7, zorder=0)

        if 1 < zero_cross < n_genes:
            ax_rank.axvline(zero_cross, color="#666666", ls=":", lw=1.1, zorder=4)
            ax_rank.text(
                zero_cross, 0.06, f" Zero cross ({zero_cross:,})",
                transform=ax_rank.get_xaxis_transform(),
                ha="left" if zero_cross < n_genes * 0.7 else "right",
                va="bottom", fontsize=9.5, color="#555555",
            )

        ax_rank.text(0.012, 0.88, f"Higher in {treatment}", transform=ax_rank.transAxes,
                     ha="left", va="top", fontsize=10.5, fontweight="bold", color=style.POS)
        ax_rank.text(0.988, 0.10, f"Higher in {control}", transform=ax_rank.transAxes,
                     ha="right", va="bottom", fontsize=10.5, fontweight="bold", color=style.NEG)

        for s in ("top", "right"):
            ax_rank.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax_rank.spines[s].set_linewidth(1.2)
        ax_rank.tick_params(labelsize=11.0)
        ax_rank.set_ylabel(rank_label, fontsize=12, labelpad=8)
        ax_rank.set_xlabel(f"Gene rank in ordered dataset (n = {n_genes:,})", fontsize=13, labelpad=8)

    return fig
