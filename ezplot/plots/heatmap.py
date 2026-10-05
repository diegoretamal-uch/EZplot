"""Clustered DEG heatmap in the format of the template deck (slide 3).

iDEP's DEG_Heatmap_Data.csv is row-centred but not scaled; rows are divided by their
standard deviation (ddof=1, like R's scale()) to give row z-scores. Rows and columns are
clustered hierarchically (average linkage, Euclidean), colour scale RdYlBu_r from -vlim
to vlim, with a condition bar above the columns.
"""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from scipy.cluster import hierarchy

from .. import style


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
    show_genes: bool = True
    italic_genes: bool = True
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


def _order(data, method, metric):
    if len(data) < 2:
        return np.arange(len(data)), None
    link = hierarchy.linkage(data, method=method, metric=metric)
    return np.array(hierarchy.leaves_list(link)), link


def plot_heatmap(x, sample_roles, role_labels, role_colors, p=None):
    """x: genes x samples. sample_roles: sample -> 'control'/'treatment'/other key.

    role_labels/role_colors: role -> legend label / colour.
    """
    p = p or HeatmapParams()
    style.setup()
    z = row_zscore(x) if p.zscore else x.copy()
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

    mesh = ax_h.pcolormesh(z.to_numpy(), cmap=p.cmap, vmin=-p.vlim, vmax=p.vlim,
                           edgecolors="#aaaaaa", linewidth=0.5)
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
