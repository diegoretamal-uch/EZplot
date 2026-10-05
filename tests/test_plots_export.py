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
