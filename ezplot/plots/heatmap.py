"""Clustered DEG heatmap in the format of the template deck (slide 3).

iDEP's DEG_Heatmap_Data.csv is row-centred but not scaled; rows are divided by their
standard deviation (ddof=1, like R's scale()) to give row z-scores. Rows and columns are
clustered hierarchically (average linkage, Euclidean), colour scale RdYlBu_r from -vlim
to vlim, with a condition bar above the columns.
"""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial import distance

from .. import style

LINKAGE_METHODS = {
    "average": "Average (UPGMA)",
    "ward": "Ward (minimum variance)",
    "complete": "Complete (maximum)",
    "single": "Single (minimum)",
    "weighted": "Weighted (WPGMA)",
}
DISTANCE_METRICS = {
    "euclidean": "Euclidean",
    "correlation": "Pearson correlation",
    "cityblock": "Manhattan (L1)",
    "cosine": "Cosine",
}
TOP_RANK_MODES = {
    "fdr": "FDR (smallest)",
    "lfc": "|log₂FC| (largest)",
    "pi": "|log₂FC| × −log₁₀(FDR)",
}


@dataclass
class HeatmapParams:
    title: str = "Differentially Expressed Genes (n = {n})"
    zscore: bool = True
    vlim: float = 2.0
    cmap: str = "RdYlBu_r"
    cluster_rows: bool = True
    cluster_cols: bool = True
    method: str = "average"
    metric: str = "euclidean"
    top_n: int = 0                 # 0 = all genes in file; > 0 = top N genes by `top_by`
    top_by: str = "fdr"            # fdr / lfc / pi
    show_genes: bool = True
    italic_genes: bool = True
    cell_borders: bool = True      # grey lines between cells (and between samples in the bar)
    gene_fontsize: float = 15
    sample_fontsize: float = 14
    title_fontsize: float = 20
    cbar_fontsize: float = 13      # colour-scale label; ticks one point smaller
    legend_fontsize: float = 13    # legend entries; title one point larger
    width: float = 10.5
    height: float = 0.0            # 0 = automatic from the number of genes
    legend_title: str = "Condition"


def row_zscore(x):
    sd = x.std(axis=1, ddof=1).replace(0, np.nan)
    return x.sub(x.mean(axis=1), axis=0).div(sd, axis=0).fillna(0.0)


def _empirical_stats(x, sample_roles):
    """Fallback |log2FC| and p-value from the heatmap matrix when DE table is not loaded."""
    ctrl = [s for s in x.columns if sample_roles.get(s) == "control"]
    trt = [s for s in x.columns if sample_roles.get(s) == "treatment"]
    sd = x.std(axis=1, ddof=1).fillna(0.0)
    if not ctrl or not trt:
        return sd, pd.Series(1.0, index=x.index)
    lfc_abs = (x[trt].mean(axis=1) - x[ctrl].mean(axis=1)).abs()
    if len(ctrl) >= 2 and len(trt) >= 2:
        _, pval = stats.ttest_ind(x[trt].to_numpy(), x[ctrl].to_numpy(), axis=1, equal_var=False)
        pval = pd.Series(np.nan_to_num(pval, nan=1.0), index=x.index).clip(lower=1e-300, upper=1.0)
    else:
        pval = pd.Series(1.0, index=x.index)
    return lfc_abs, pval


def select_top_genes(x, sample_roles, p, de_genes=None):
    """Filter heatmap matrix `x` to the top `p.top_n` genes by `fdr`, `lfc`, or `pi`."""
    if not p.top_n or p.top_n <= 0 or len(x) <= p.top_n:
        return x
    k = int(p.top_n)
    by = p.top_by or "fdr"
    emp_lfc, emp_p = _empirical_stats(x, sample_roles)

    if de_genes is not None and not de_genes.empty:
        de = (
            de_genes.dropna(subset=["symbol"])
            .assign(_abs_lfc=de_genes["lfc"].abs())
            .sort_values(["padj", "_abs_lfc"], ascending=[True, False], kind="mergesort")
            .drop_duplicates("symbol", keep="first")
            .set_index("symbol")
        )
        lfc_abs = de["lfc"].abs().reindex(x.index).fillna(emp_lfc)
        padj = de["padj"].reindex(x.index).fillna(emp_p)
    else:
        lfc_abs, padj = emp_lfc, emp_p

    pos_p = padj[padj > 0]
    floor = pos_p.min() / 10.0 if len(pos_p) else 1e-300
    neg_log10_fdr = -np.log10(padj.clip(lower=floor, upper=1.0))

    if by == "lfc":
        rank_df = pd.DataFrame({"_k1": -lfc_abs, "_k2": padj}, index=x.index)
    elif by == "pi":
        pi_score = lfc_abs * neg_log10_fdr
        rank_df = pd.DataFrame({"_k1": -pi_score, "_k2": padj}, index=x.index)
    else:  # "fdr"
        rank_df = pd.DataFrame({"_k1": padj, "_k2": -lfc_abs}, index=x.index)

    chosen = rank_df.sort_values(["_k1", "_k2"], kind="mergesort").head(k).index
    return x.loc[x.index.isin(chosen)]


def _order(data, method, metric):
    if len(data) < 2:
        return np.arange(len(data)), None
    d = distance.pdist(data, metric=metric)
    d = np.nan_to_num(d, nan=0.0, posinf=0.0, neginf=0.0)
    link = hierarchy.linkage(d, method=method)
    return np.array(hierarchy.leaves_list(link)), link


def plot_heatmap(x, sample_roles, role_labels, role_colors, p=None, de_genes=None):
    """x: genes x samples. sample_roles: sample -> 'control'/'treatment'/other key.

    role_labels/role_colors: role -> legend label / colour.
    """
    p = p or HeatmapParams()
    style.setup()
    x_sub = select_top_genes(x, sample_roles, p, de_genes=de_genes)
    z = row_zscore(x_sub) if p.zscore else x_sub.copy()
    n, m = z.shape

    ri, rlink = _order(z.to_numpy(), p.method, p.metric) if p.cluster_rows else (np.arange(n), None)
    ci, clink = _order(z.to_numpy().T, p.method, p.metric) if p.cluster_cols else (np.arange(m), None)
    z = z.iloc[ri, ci]

    height = p.height or 4.0 + 0.26 * n
    fig = plt.figure(figsize=(p.width, height))
    left, width, bottom, top = 0.11, 0.62, 0.13, 0.70
    ax_h = fig.add_axes([left, bottom, width, top - bottom])
    ax_c = fig.add_axes([left, top + 0.006, width, 0.025])
    ax_cd = fig.add_axes([left, top + 0.036, width, 0.13])
    ax_rd = fig.add_axes([0.01, bottom, left - 0.015, top - bottom])
    ax_cb = fig.add_axes([0.93, bottom + 0.02, 0.018, 0.50])

    # without borders, edges take the cell colour so PDF/SVG viewers show no hairline seams
    mesh = ax_h.pcolormesh(z.to_numpy(), cmap=p.cmap, vmin=-p.vlim, vmax=p.vlim,
                           edgecolors="#aaaaaa" if p.cell_borders else "face", linewidth=0.5)
    ax_h.set_xlim(0, m)
    ax_h.set_ylim(n, 0)
    ax_h.set_xticks(np.arange(m) + 0.5, z.columns, rotation=90, fontsize=p.sample_fontsize)
    if p.show_genes:
        ax_h.yaxis.tick_right()
        ax_h.set_yticks(np.arange(n) + 0.5, z.index, fontsize=p.gene_fontsize,
                        fontstyle="italic" if p.italic_genes else "normal")
    else:
        ax_h.set_yticks([])
    ax_h.tick_params(length=0, pad=6)
    ax_h.tick_params(axis="x", pad=8)
    for s in ax_h.spines.values():
        s.set_color("#777777")

    roles = [sample_roles.get(s) for s in z.columns]
    colors = [role_colors.get(r, "#bbbbbb") for r in roles]
    ax_c.imshow([[plt.matplotlib.colors.to_rgb(c) for c in colors]], aspect="auto",
                extent=(0, m, 0, 1))
    if p.cell_borders:
        ax_c.vlines(np.arange(1, m), 0, 1, color="#aaaaaa", lw=0.5)
    ax_c.set_xticks([])
    ax_c.set_yticks([])
    for s in ax_c.spines.values():
        s.set_color("#999999")

    tree = dict(no_labels=True, color_threshold=0, above_threshold_color="#333333")
    for ax, link, orient in ((ax_cd, clink, "top"), (ax_rd, rlink, "left")):
        if link is not None:
            with plt.rc_context({"lines.linewidth": 1.2}):
                hierarchy.dendrogram(link, ax=ax, orientation=orient, **tree)
            if orient == "left":
                ax.set_ylim(10 * n, 0)
        ax.set_axis_off()

    cb = fig.colorbar(mesh, cax=ax_cb)
    cb.set_ticks(np.linspace(-p.vlim, p.vlim, 5) if p.vlim == 2 else
                 [-p.vlim, 0, p.vlim])
    cb.ax.tick_params(labelsize=p.cbar_fontsize - 1)
    cb.set_label("Row z-score" if p.zscore else "Centred log$_2$ expression", fontsize=p.cbar_fontsize,
                 labelpad=10)
    cb.outline.set_visible(False)

    fig.suptitle(p.title.format(n=n), fontsize=p.title_fontsize, fontweight="bold",
                 x=left + width / 2, y=0.955)
    present = [r for r in role_labels if r in roles]
    fig.legend(handles=[Patch(color=role_colors[r], label=role_labels[r]) for r in present],
               title=p.legend_title,
               title_fontproperties={"weight": "bold", "size": p.legend_fontsize + 1},
               fontsize=p.legend_fontsize, loc="upper left", bbox_to_anchor=(0.76, 0.93), frameon=False,
               handlelength=1.4, handleheight=0.9)
    return fig
