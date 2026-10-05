"""GSEA dot plot in the format of the template deck (slides 2-3).

Top-N gene sets with FDR q < cutoff (smallest q first), x = NES, colour = -log10(q) on
viridis clipped at qmax, area = gene-set size.
"""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

from .. import style

@dataclass
class GseaDotParams:
    title: str = "{collection} — FDR q < {q:g}"
    subtitle: str = "Positive NES: {treatment} · Negative NES: {control}"
    q: float = 0.05
    top: int = 10
    qmax: float = 4.0
    size_k: float = 1.4
    wrap: int = 50
    cmap: str = "viridis"
    width: float = 10.0
    height: float = 0.0          # 0 = automatic from the number of sets
    direction: str = "both"      # both / positive / negative


def size_legend_values(max_size):
    for vals in ([5, 15, 25], [10, 50, 100], [25, 100, 200], [25, 100, 300],
                 [25, 100, 400], [50, 200, 700]):
        if max_size <= vals[-1] * 1.2:
            return vals
    return [100, 400, 800]


def select(res, p):
    sig = res[res["padj"] < p.q]
    if p.direction == "positive":
        sig = sig[sig["NES"] > 0]
    elif p.direction == "negative":
        sig = sig[sig["NES"] < 0]
    return sig.head(p.top).sort_values("NES")   # res is already sorted by padj, |NES|


def plot_gsea_dot(res, collection_title, ylabel, control, treatment, p=None):
    """res: engine output (ezplot.gsea.engine.prerank)."""
    p = p or GseaDotParams()
    style.setup()
    sig = select(res, p)
    n = len(sig)
    height = p.height or max(5.93, 2.4 + 0.42 * n)
    fig = plt.figure(figsize=(p.width, height))
    # fixed physical margins (inches) so tall figures keep the template proportions
    ax_b, ax_t = 1.19 / height, 1 - 1.07 / height
    ax = fig.add_axes([0.50, ax_b, 0.365, ax_t - ax_b])
    cb_top = ax_t - 0.18 / height
    cb_h = min(1.96 / height, cb_top - ax_b)
    cax = fig.add_axes([0.89, cb_top - cb_h, 0.017, cb_h])

    fig.text(0.5, 1 - 0.27 / height, p.title.format(collection=collection_title, q=p.q),
             ha="center", va="top", fontsize=17, color=style.TEXT)
    fig.text(0.5, 1 - 0.62 / height, p.subtitle.format(treatment=treatment, control=control),
             ha="center", va="top", fontsize=12.5, color="#333333")

    lim = max(2.6, np.ceil((sig["NES"].abs().max() if n else 0) * 5 + 1.5) / 5)
    ax.set_xlim(-lim, lim)
    ax.axvline(0, color=style.ZERO_LINE, ls="--", lw=1.8, zorder=1)
    ax.grid(True, ls="--", color="#cccccc", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("Normalized enrichment score (NES)", fontsize=13.5, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=13.5, labelpad=12)
    ax.tick_params(labelsize=12.5)

    norm = Normalize(-np.log10(p.q), p.qmax)
    if n:
        y = np.arange(n)[::-1]
        logq = np.minimum(-np.log10(sig["padj"].clip(lower=1e-300)), p.qmax)
        ax.scatter(sig["NES"], y, s=p.size_k * sig["size"], c=logq, cmap=p.cmap,
                   norm=norm, edgecolors="none", zorder=3)
        labels = [style.pretty_name(name, p.wrap) for name in sig["pathway"]]
        ax.set_yticks(y, labels)
        if any("\n" in l for l in labels):
            ax.tick_params(axis="y", labelsize=11)
        ax.set_ylim(-0.6, n - 0.4)

        vals = size_legend_values(sig["size"].max())
        handles = [ax.scatter([], [], s=p.size_k * v, color="#808080", edgecolors="none")
                   for v in vals]
        ax.legend(handles, vals, title="Gene-set size", loc="upper left",
                  bbox_to_anchor=(0.875, cb_top - cb_h - 0.35 / height), bbox_transform=fig.transFigure,
                  fontsize=10.5, title_fontsize=10.5, labelspacing=1.0, borderpad=0.7,
                  handletextpad=1.0, framealpha=1, edgecolor="#cccccc")
    else:
        ax.set_yticks([])
        ax.set_ylim(0, 1)
        ax.text(0.5, 0.5, f"No gene sets at FDR q < {p.q:g}", transform=ax.transAxes,
                ha="center", va="center", fontsize=13, color=style.MUTED,
                bbox=dict(fc="white", ec="none"))

    cb = fig.colorbar(plt.cm.ScalarMappable(norm, p.cmap), cax=cax)
    ticks = [-np.log10(p.q)] + [t for t in range(2, int(p.qmax) + 1) if t > -np.log10(p.q) + 0.3]
    cb.set_ticks(ticks, labels=[f"{ticks[0]:.2f}"] + [str(t) for t in ticks[1:]])
    cb.ax.tick_params(labelsize=11.5)
    cb.set_label("$-$log$_{10}$(FDR $q$$-$value)\n"
                 f"($q$ < 10$^{{-{p.qmax:g}}}$ shown as 10$^{{-{p.qmax:g}}}$)",
                 fontsize=12.5, labelpad=8)

    return fig
