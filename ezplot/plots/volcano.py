"""Volcano plot in the format of the template deck (slide 3, "Differential gene expression").

x = log2FC (oriented: positive = higher in treatment), y = -log10(adjusted P). Genes with
FDR < cutoff and |log2FC| >= threshold are coloured; the rest are grey. Significant genes
are labelled (italic) with leader lines, placed by adjustText.
"""

from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from .. import style


@dataclass
class VolcanoParams:
    title: str = "Differential gene expression"
    fdr: float = 0.05
    lfc: float = 1.0
    label_mode: str = "significant"     # significant / top / custom / none
    label_top: int = 20
    label_rank: str = "fdr"             # top N by: fdr / lfc / distance
    label_genes: list = field(default_factory=list)
    neg_color: str = style.NEG
    pos_color: str = style.POS
    ns_color: str = style.NS
    neg_label: str = "Negative log$_2$FC"
    pos_label: str = "Positive log$_2$FC"
    show_thresholds: bool = True
    xlim: float = 0.0                   # 0 = automatic, symmetric
    ymax: float = 0.0                   # 0 = automatic
    point_size: float = 70
    label_fontsize: float = 13
    legend_fontsize: float = 12.5
    width: float = 8.1
    height: float = 7.0


def classify(df, p):
    sig = (df["padj"] < p.fdr) & (df["lfc"].abs() >= p.lfc)
    return np.where(sig & (df["lfc"] > 0), "pos", np.where(sig & (df["lfc"] < 0), "neg", "ns"))


def plot_volcano(genes, p=None):
    """genes: DataFrame with columns symbol, lfc (oriented), padj."""
    from adjustText import adjust_text

    p = p or VolcanoParams()
    style.setup()
    df = genes[["symbol", "lfc", "padj"]].dropna().copy()
    floor = df.loc[df["padj"] > 0, "padj"].min() if (df["padj"] > 0).any() else 1e-300
    df["y"] = -np.log10(df["padj"].clip(lower=floor / 10))
    df["cls"] = classify(df, p)
    counts = df["cls"].value_counts()

    fig = plt.figure(figsize=(p.width, p.height))
    ax = fig.add_axes([0.09, 0.10, 0.88, 0.66])

    ns = df[df["cls"] == "ns"]
    ax.scatter(ns["lfc"], ns["y"], s=p.point_size * 0.35, color=p.ns_color, alpha=0.8,
               edgecolors="none", zorder=2, rasterized=len(ns) > 3000)
    for cls, color in (("neg", p.neg_color), ("pos", p.pos_color)):
        d = df[df["cls"] == cls]
        ax.scatter(d["lfc"], d["y"], s=p.point_size, color=color, edgecolors="white",
                   linewidths=0.8, zorder=3)

    xl = p.xlim or np.ceil(df["lfc"].abs().max() * 1.05 * 2) / 2
    ax.set_xlim(-xl, xl)
    ytop = p.ymax or df["y"].max() * 1.08
    ax.set_ylim(0, max(ytop, -np.log10(p.fdr) * 1.5))
    if p.show_thresholds:
        kw = dict(color="#888888", ls=(0, (3, 2)), lw=1.3, zorder=1)
        ax.axhline(-np.log10(p.fdr), **kw)
        if p.lfc > 0:
            ax.axvline(-p.lfc, **kw)
            ax.axvline(p.lfc, **kw)
        ax.text(0.995, 0.012, f"FDR < {p.fdr:g}; |log$_2$FC| ≥ {p.lfc:g}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=13,
                color="#6b7078", bbox=dict(fc="white", ec="none", pad=1, alpha=0.85))

    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.4)
    ax.tick_params(labelsize=14, width=1.4, length=6)
    ax.set_xlabel("log$_2$ fold change", fontsize=18, labelpad=8)
    ax.set_ylabel("$-$log$_{10}$ (adjusted $P$ value)", fontsize=18, labelpad=8)

    fig.text(0.09, 0.975, p.title, ha="left", va="top", fontsize=22, fontweight="bold",
             color=style.TEXT)
    handles = [
        Line2D([], [], ls="", marker="o", ms=9, color=p.neg_color,
               label=f"{p.neg_label} (n={counts.get('neg', 0):,})"),
        Line2D([], [], ls="", marker="o", ms=9, color=p.ns_color,
               label=f"Not significant (n={counts.get('ns', 0):,})"),
        Line2D([], [], ls="", marker="o", ms=9, color=p.pos_color,
               label=f"{p.pos_label} (n={counts.get('pos', 0):,})"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.02, 0.84, 0.97, 0.05),
               mode="expand", ncol=3, frameon=False, fontsize=p.legend_fontsize,
               handletextpad=0.3, columnspacing=0.8, borderaxespad=0)

    lab = _labels(df, p)
    if len(lab):
        texts = []
        for _, r in lab.iterrows():
            dx = 0.04 * xl * (1 if r["lfc"] > 0 else -1)
            texts.append(ax.text(r["lfc"] + dx, r["y"], r["symbol"], fontsize=p.label_fontsize,
                                 fontstyle="italic", color=style.TEXT, zorder=5,
                                 ha="left" if r["lfc"] > 0 else "right", va="center"))
        sig = df[df["cls"] != "ns"]
        adjust_text(texts, target_x=lab["lfc"].to_numpy(), target_y=lab["y"].to_numpy(),
                    objects=None, ax=ax, x=sig["lfc"].to_numpy(), y=sig["y"].to_numpy(),
                    expand=(1.3, 1.6), force_text=(0.4, 0.6), force_static=(0.3, 0.5),
                    only_move={"text": "xy", "static": "xy", "explode": "xy", "pull": "xy"},
                    arrowprops=dict(arrowstyle="-", color="#8a8f97", lw=1.0,
                                    shrinkA=2, shrinkB=4),
                    min_arrow_len=4)
    return fig


def _labels(df, p):
    if p.label_mode == "none":
        return df.iloc[:0]
    if p.label_mode == "custom":
        want = {g.strip().lower() for g in p.label_genes if g.strip()}
        return df[df["symbol"].str.lower().isin(want)]
    sig = df[df["cls"] != "ns"]
    if p.label_mode == "top":
        return rank_genes(sig, p.label_rank).head(p.label_top)
    return sig.sort_values("padj").head(60)   # "significant": all, capped for legibility


def rank_genes(d, by="fdr"):
    """Order genes for labelling: smallest FDR, largest |log2FC|, or farthest from the
    origin as drawn (both axes scaled to their range, so neither dominates)."""
    if by == "lfc":
        key = -d["lfc"].abs()
    elif by == "distance":
        x = d["lfc"].abs() / (d["lfc"].abs().max() or 1)
        y = d["y"] / (d["y"].max() or 1)
        key = -(x ** 2 + y ** 2)
    else:
        key = d["padj"]
    return d.assign(_k=key).sort_values(["_k", "padj"], kind="mergesort").drop(columns="_k")
