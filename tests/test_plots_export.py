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


def test_gsea_size_legend_constant_gap_and_default_nperm():
    import numpy as np
    import pandas as pd

    from ezplot.plots.gsea_dot import plot_gsea_dot
    from ezplot.project import GseaSettings

    assert GseaSettings().nperm == 10000

    res = pd.DataFrame({
        "pathway": [f"HALLMARK_SET_{i}" for i in range(6)],
        "ES": [0.5, 0.4, 0.3, -0.3, -0.4, -0.5],
        "NES": [2.2, 1.9, 1.6, -1.6, -1.9, -2.2],
        "pval": [1e-5] * 6,
        "padj": [1e-4, 5e-4, 2e-3, 3e-3, 1e-2, 2e-2],
        "size": [690, 420, 250, 180, 90, 50],
        "leadingEdge": ["A;B"] * 6,
    })
    fig = plot_gsea_dot(res, "Hallmarks", "Pathway", "Ctrl", "Trt", GseaDotParams())
    # Legend axes is created before colorbar
    lax = fig.axes[2]
    ys = [c.get_offsets()[0, 1] for c in lax.collections]
    radii = [0.5 * np.sqrt(c.get_sizes()[0]) for c in lax.collections]
    gaps = [(ys[i] - radii[i]) - (ys[i + 1] + radii[i + 1]) for i in range(len(ys) - 1)]
    assert np.allclose(gaps, 9.0)
    assert all(g > 0 for g in gaps)
    plt.close(fig)


def test_heatmap_top_n_and_clustering():
    from ezplot.plots.heatmap import DISTANCE_METRICS, LINKAGE_METHODS, TOP_RANK_MODES

    pr = project()
    for top_by in TOP_RANK_MODES:
        fig = pr.heatmap_fig(HeatmapParams(top_n=8, top_by=top_by))
        assert len(fig.axes[0].get_yticklabels()) == 8
        assert "(n = 8)" in fig._suptitle.get_text()
        plt.close(fig)

    for method in LINKAGE_METHODS:
        for metric in DISTANCE_METRICS:
            fig = pr.heatmap_fig(HeatmapParams(top_n=10, method=method, metric=metric))
            assert len(fig.axes[0].get_yticklabels()) == 10
            plt.close(fig)


def test_heatmap_directly_from_de_matrix():
    pr = Project()
    pr.load(EXAMPLE / "example_DE_matrix.csv", "de.csv")
    assert pr.heatmap is None and pr.has_de_expr and pr.has_heatmap
    assert pr.heatmap_source_used() == "de"

    hm_de = pr.heatmap_matrix(HeatmapParams(source="de", fdr=0.05, lfc=1.0))
    pr_both = project()
    assert set(hm_de.index) == set(pr_both.heatmap.index)

    fig = pr.heatmap_fig(HeatmapParams(source="de", fdr=0.05, lfc=1.0, top_n=12))
    assert len(fig.axes[0].get_yticklabels()) == 12
    assert "(n = 12)" in fig._suptitle.get_text()
    plt.close(fig)

    # Empty filter branch renders cleanly without raising
    fig_empty = pr.heatmap_fig(HeatmapParams(source="de", fdr=1e-50, lfc=50.0))
    assert "(n = 0)" in fig_empty._suptitle.get_text()
    plt.close(fig_empty)


def test_gsea_mountain_single_and_multi():
    from ezplot.plots.gsea_mountain import GseaMountainParams

    pr = Project()
    pr.load(EXAMPLE / "example_DE_matrix.csv", "de.csv")

    # Single-pathway mountain plot (default top pathway + leading-edge gene labels)
    fig1 = pr.gsea_mountain_fig("hallmarks", GseaMountainParams(label_genes=5))
    assert len(fig1.axes) == 3
    assert len(export.figure_bytes(fig1, "png", 120)) > 1000
    assert len(export.figure_bytes(fig1, "svg", 120)) > 1000
    plt.close(fig1)

    # Multi-pathway mountain plot without bottom ranking panel
    fig_multi = pr.gsea_mountain_fig(
        "hallmarks",
        GseaMountainParams(
            pathways=[
                "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            ],
            show_ranking=False,
        ),
    )
    assert len(fig_multi.axes) == 2
    assert len(fig_multi.legends) == 1
    assert len(fig_multi.legends[0].get_texts()) == 3
    plt.close(fig_multi)


def _tall_figure(w=5.0, h=5.0):
    fig = plt.figure(figsize=(w, h))
    fig.add_axes([0, 0, 1, 1]).plot([0, 1], [0, 1])
    return fig


def test_figure_bytes_no_decompression_bomb(monkeypatch):
    # Real case: a ~750-gene heatmap at 300 dpi (or most heatmaps at 1200 dpi) is more than
    # 2 x PIL's MAX_IMAGE_PIXELS. Shrink the limit so a small figure reproduces it cheaply.
    monkeypatch.setattr(Image, "MAX_IMAGE_PIXELS", 100_000)   # bomb error above 200k px
    fig = _tall_figure()                                       # 5 x 5 in at 100 dpi ~ 250k px
    raw = _io.BytesIO()
    fig.savefig(raw, format="png", dpi=100, bbox_inches="tight", pad_inches=0.08)
    try:
        Image.open(_io.BytesIO(raw.getvalue()))
        raise AssertionError("expected PIL to reject this image size")
    except Image.DecompressionBombError:
        pass                                                   # the old code path failed here

    for fmt, width in (("png", None), ("png", 120), ("tiff", None), ("tiff", 120)):
        data = export.figure_bytes(fig, fmt, 100, width)
        assert Image.MAX_IMAGE_PIXELS == 100_000               # guard restored afterwards
        monkeypatch.setattr(Image, "MAX_IMAGE_PIXELS", None)
        img = Image.open(_io.BytesIO(data))
        img.load()
        monkeypatch.setattr(Image, "MAX_IMAGE_PIXELS", 100_000)
        assert img.format == ("PNG" if fmt == "png" else "TIFF")
        assert round(img.info["dpi"][0]) == 100
        if width:
            assert img.width == round(width / 25.4 * 100)
        else:
            assert img.width * img.height > 200_000
    plt.close(fig)


def test_figure_bytes_caps_raster_size(monkeypatch):
    # Agg cannot render 2**16 px per side; huge figures are rendered at a lower DPI instead
    monkeypatch.setattr(export, "RASTER_MAX_PX", 300)
    fig = _tall_figure(3.0, 6.0)                               # 6 in at 100 dpi = 600 px
    for fmt in ("png", "tiff"):
        img = Image.open(_io.BytesIO(export.figure_bytes(fig, fmt, 100)))
        assert max(img.size) <= 300 * 1.1                      # tight bbox may add the pad
        assert round(img.info["dpi"][0]) == 50
    plt.close(fig)


def test_bundle_survives_failing_items():
    import pandas as pd

    def broken():
        raise RuntimeError("plot exploded")

    def half_built():
        plt.figure()
        raise ValueError("failed after creating a figure")

    class BadTable:
        def to_csv(self, **kw):
            raise OSError("cannot serialise")

    open_before = set(plt.get_fignums())
    data = export.bundle(
        {"good": _tall_figure, "broken": broken, "half": half_built},
        {"ok": pd.DataFrame({"a": [1, 2]}), "bad": BadTable()},
        {"x": 1}, ("png", "pdf"), 100,
        errors=["GSEA GO Biological Process: MemoryError: out of memory"],
    )
    assert set(plt.get_fignums()) == open_before               # nothing leaked
    assert len(data) > 1000
    z = zipfile.ZipFile(_io.BytesIO(data))
    assert z.testzip() is None
    names = set(z.namelist())
    assert {"figures/good.png", "figures/good.pdf", "tables/ok.tsv", "ezplot_parameters.json",
            "errors.txt"} <= names
    assert not any("broken" in n or "half" in n or "bad" in n for n in names - {"errors.txt"})
    txt = z.read("errors.txt").decode()
    for needle in ("figures/broken", "RuntimeError: plot exploded", "figures/half",
                   "ValueError", "tables/bad.tsv", "OSError", "GO Biological Process"):
        assert needle in txt
    import json
    assert len(json.loads(z.read("ezplot_parameters.json"))["errors"]) == 4


def test_bundle_without_errors_has_no_errors_file():
    data = export.bundle({"good": _tall_figure}, {}, {}, ("png",), 80)
    names = zipfile.ZipFile(_io.BytesIO(data)).namelist()
    assert "errors.txt" not in names and "figures/good.png" in names
