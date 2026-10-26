import io as _io
import zipfile

import matplotlib.pyplot as plt
from PIL import Image

from ezplot import export
from ezplot.plots.gsea_dot import GseaDotParams
from ezplot.plots.heatmap import HeatmapParams
from ezplot.plots.volcano import VolcanoParams
from ezplot.project import Project

from .conftest import EXAMPLE


def project():
    pr = Project()
    pr.load(EXAMPLE / "example_DE_matrix.csv", "de.csv")
    pr.load(EXAMPLE / "example_DEG_Heatmap_Data.csv", "hm.csv")
    return pr


def test_all_figures_render_and_export():
    pr = project()
    for make in (pr.volcano, pr.heatmap_fig, lambda: pr.gsea_fig("hallmarks")):
        fig = make()
        for fmt in export.FORMATS:
            assert len(export.figure_bytes(fig, fmt, 150)) > 1000
        plt.close(fig)


def test_raster_width_and_dpi():
    fig = project().volcano()
    img = Image.open(_io.BytesIO(export.figure_bytes(fig, "tiff", 300, width_mm=89)))
    assert img.width == round(89 / 25.4 * 300)
    assert round(img.info["dpi"][0]) == 300
    plt.close(fig)


def test_heatmap_without_cell_borders():
    pr = project()
    fig = pr.heatmap_fig(HeatmapParams(cell_borders=False))
    mesh = fig.axes[0].collections[0]
    assert (mesh.get_edgecolor() == mesh.get_facecolor()).all()   # edges blend into cells
    assert not fig.axes[1].collections                              # no separators in the bar
    plt.close(fig)


def test_volcano_label_modes_and_empty_gsea():
    pr = project()
    for mode in ("significant", "top", "custom", "none"):
        plt.close(pr.volcano(VolcanoParams(label_mode=mode, label_genes=["Emc8"])))
    plt.close(pr.gsea_fig("hallmarks", GseaDotParams(q=1e-30)))   # "no sets" branch


def test_bundle():
    pr = project()
    data = export.bundle({"volcano": pr.volcano}, {"t": pr.run_gsea("hallmarks")},
                         {"x": 1}, ("png", "svg"), 100)
    names = zipfile.ZipFile(_io.BytesIO(data)).namelist()
    assert {"figures/volcano.png", "figures/volcano.svg", "tables/t.tsv",
            "ezplot_parameters.json"} <= set(names)


def test_heatmap_font_sizes():
    fig = project().heatmap_fig(HeatmapParams(title_fontsize=31, cbar_fontsize=17,
                                              legend_fontsize=11))
    assert fig._suptitle.get_fontsize() == 31
    cbar = fig.axes[-1]
    assert cbar.yaxis.label.get_fontsize() == 17
    assert cbar.get_yticklabels()[0].get_fontsize() == 16
    leg = fig.legends[0]
    assert leg.get_texts()[0].get_fontsize() == 11
    assert leg.get_title().get_fontsize() == 12
    plt.close(fig)


def test_volcano_label_ranking():
    import pandas as pd

    from ezplot.plots.volcano import rank_genes
    d = pd.DataFrame({"symbol": ["a", "b", "c"], "lfc": [-1.2, 6.0, 3.0],
                      "padj": [1e-10, 1e-3, 1e-6], "y": [10.0, 3.0, 6.0]})
    assert list(rank_genes(d, "fdr").symbol) == ["a", "c", "b"]
    assert list(rank_genes(d, "lfc").symbol) == ["b", "c", "a"]
    # scaled (x/6, y/10): b=(1.0,0.3) 1.09 > a=(0.2,1.0) 1.04 > c=(0.5,0.6) 0.61
    assert list(rank_genes(d, "distance").symbol) == ["b", "a", "c"]
    pr = project()
    for by in ("fdr", "lfc", "distance"):
        plt.close(pr.volcano(VolcanoParams(label_mode="top", label_top=5, label_rank=by)))


def test_volcano_labels_no_overlap_high_n():
    import numpy as np
    import pandas as pd

    from ezplot.plots.volcano import plot_volcano

    rng = np.random.default_rng(42)
    n = 1500
    lfc = rng.normal(0, 1.2, n)
    padj = np.clip(10 ** (-np.abs(lfc) * rng.uniform(0.5, 2.5, n) + rng.normal(0, 0.5, n)), 1e-15, 1.0)
    genes = pd.DataFrame({"symbol": [f"Gene{i:03d}" for i in range(n)], "lfc": lfc, "padj": padj})

    fig = plot_volcano(genes, VolcanoParams(label_mode="top", label_top=60))
    ax = fig.axes[0]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    gene_texts = [t for t in ax.texts if t.get_text().startswith("Gene")]
    assert len(gene_texts) == 60

    bboxes = [t.get_window_extent(renderer) for t in gene_texts]
    overlaps = 0
    for i in range(len(bboxes)):
        for j in range(i + 1, len(bboxes)):
            a, b = bboxes[i], bboxes[j]
            if (min(a.xmax, b.xmax) - max(a.xmin, b.xmin) > 1.0) and (
                min(a.ymax, b.ymax) - max(a.ymin, b.ymin) > 1.0
            ):
                overlaps += 1
    assert overlaps == 0
    plt.close(fig)

