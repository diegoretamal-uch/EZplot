"""EZplot Shiny app: load iDEP files, preview and tune the plots, export for a journal."""

import base64
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from shiny import App, reactive, render, ui  # noqa: E402

from . import __version__, export, genesets, io, style  # noqa: E402
from .plots.gsea_dot import GseaDotParams  # noqa: E402
from .plots.heatmap import HeatmapParams  # noqa: E402
from .plots.volcano import VolcanoParams  # noqa: E402
from .project import GseaSettings, Project  # noqa: E402

EXAMPLE = Path(__file__).parent / "example"
CSS = """
.ez-drop .form-group{margin-bottom:.25rem}
.ez-badge{display:inline-block;font-size:.78rem;padding:.1rem .45rem;border-radius:.6rem;
  background:var(--bs-secondary-bg);margin:.1rem .2rem .1rem 0}
.ez-badge.ok{background:#d9f2e6;color:#11603b}.ez-badge.warn{background:#fdeccf;color:#7a4b00}
.ez-preview img{max-width:100%;height:auto;border:1px solid var(--bs-border-color);
  border-radius:.4rem;background:#fff}
.ez-empty{color:var(--bs-secondary-color);padding:3rem 1rem;text-align:center}
.ez-dir{font-weight:600}.ez-small{font-size:.82rem;color:var(--bs-secondary-color)}
input[type=color].ez-color{width:2.6rem;height:1.9rem;padding:0;border:none;background:none;
  vertical-align:middle}
.ez-colors label{margin-right:.4rem}
"""
JS = """
$(document).on('input change', 'input[type=color].ez-color', function(e){
  Shiny.setInputValue(e.target.id, e.target.value);
});
$(document).on('shiny:connected', function(){
  $('input[type=color].ez-color').each(function(){ Shiny.setInputValue(this.id, this.value); });
});
"""


def color_input(id, label, value):
    return ui.span(ui.tags.label(label, **{"for": id}),
                   ui.tags.input(type="color", id=id, value=value, class_="ez-color"),
                   class_="ez-colors me-2")


def fig_img(fig, dpi=110):
    try:
        png = export.figure_bytes(fig, "png", dpi)
    finally:
        plt.close(fig)
    uri = "data:image/png;base64," + base64.b64encode(png).decode()
    return ui.div(ui.img(src=uri), class_="ez-preview")


def empty(msg):
    return ui.div(ui.markdown(msg), class_="ez-empty")


def size_inputs(prefix, w, h):
    return ui.layout_columns(
        ui.input_numeric(f"{prefix}_w", "Width (in)", w, min=3, max=30, step=0.5),
        ui.input_numeric(f"{prefix}_h", "Height (in, 0 = auto)" if not h else "Height (in)",
                         h, min=0, max=40, step=0.5),
    )


def download_row(id):
    return ui.div(ui.download_button(id, "Download this plot", class_="btn-sm btn-primary"),
                  ui.span(" format, DPI and width are set in the Export tab", class_="ez-small"),
                  class_="mt-2")


# ── UI ────────────────────────────────────────────────────────

sidebar = ui.sidebar(
    ui.h6("1 · Load iDEP files"),
    ui.div(ui.input_file("files", None, multiple=True, accept=[".csv", ".tsv", ".txt", ".gmt"],
                         button_label="Browse…", placeholder="Drop files here"),
           class_="ez-drop"),
    ui.output_ui("loaded"),
    ui.input_action_link("example", "Try it with example data"),
    ui.hr(),
    ui.h6("2 · Experiment"),
    ui.output_ui("experiment"),
    width=340,
    open="always",
)

volcano_tab = ui.nav_panel(
    "Volcano",
    ui.layout_columns(
        ui.div(ui.output_ui("volcano_out"), download_row("dl_volcano")),
        ui.accordion(
            ui.accordion_panel(
                "Thresholds & labels",
                ui.layout_columns(
                    ui.input_numeric("v_fdr", "FDR cutoff", 0.05, min=0, max=1, step=0.01),
                    ui.input_numeric("v_lfc", "log₂FC cutoff", 1, min=0, step=0.25)),
                ui.input_radio_buttons("v_mode", "Gene labels",
                                       {"significant": "All significant", "top": "Top N",
                                        "custom": "Custom list", "none": "None"}, inline=True),
                ui.panel_conditional("input.v_mode === 'top'",
                                     ui.input_numeric("v_top", "N (by FDR)", 20, min=1, max=200)),
                ui.panel_conditional("input.v_mode === 'custom'",
                                     ui.input_text_area("v_genes", "Genes (one per line or comma)",
                                                        rows=4)),
            ),
            ui.accordion_panel(
                "Design",
                ui.input_text("v_title", "Title", VolcanoParams.title),
                ui.div(color_input("v_neg", "Negative", style.NEG),
                       color_input("v_pos", "Positive", style.POS),
                       color_input("v_ns", "Not significant", style.NS), class_="mb-2"),
                ui.layout_columns(
                    ui.input_numeric("v_xlim", "x limit (0 = auto)", 0, min=0, step=0.5),
                    ui.input_numeric("v_ymax", "y max (0 = auto)", 0, min=0, step=1)),
                ui.layout_columns(
                    ui.input_numeric("v_ptsize", "Point size", 70, min=5, max=300, step=5),
                    ui.input_numeric("v_lfs", "Label font size", 13, min=5, max=30)),
                size_inputs("v", VolcanoParams.width, VolcanoParams.height),
            ),
            open="Thresholds & labels",
        ),
        col_widths=(8, 4),
    ),
)

heatmap_tab = ui.nav_panel(
    "Heatmap",
    ui.layout_columns(
        ui.div(ui.output_ui("heatmap_out"), download_row("dl_heatmap")),
        ui.accordion(
            ui.accordion_panel(
                "Options",
                ui.input_text("h_title", "Title ({n} = number of genes)", HeatmapParams.title),
                ui.input_checkbox("h_z", "Scale each gene to z-scores (as in the template)", True),
                ui.input_numeric("h_vlim", "Colour scale limit (±)", 2, min=0.5, max=10, step=0.5),
                ui.layout_columns(ui.input_checkbox("h_crow", "Cluster genes", True),
                                  ui.input_checkbox("h_ccol", "Cluster samples", True)),
                ui.layout_columns(ui.input_checkbox("h_genes", "Show gene names", True),
                                  ui.input_checkbox("h_ital", "Italic genes", True)),
                ui.input_select("h_cmap", "Colour map",
                                {"RdYlBu_r": "Red-yellow-blue (template)", "RdBu_r": "Red-blue",
                                 "vlag": "Blue-white-red (soft)", "viridis": "Viridis",
                                 "PiYG_r": "Pink-green"}),
                ui.input_text("h_legend", "Legend title", HeatmapParams.legend_title),
                ui.layout_columns(
                    ui.input_numeric("h_gfs", "Gene font", 15, min=4, max=30),
                    ui.input_numeric("h_sfs", "Sample font", 14, min=4, max=30)),
                size_inputs("h", HeatmapParams.width, 0),
            ),
            open=True,
        ),
        col_widths=(8, 4),
    ),
)

gsea_tab = ui.nav_panel(
    "GSEA",
    ui.layout_columns(
        ui.div(
            ui.input_select("g_coll", "Collection", {}),
            ui.output_ui("gsea_out"),
            download_row("dl_gsea"),
            ui.hr(),
            ui.h6("Results table"),
            ui.output_data_frame("gsea_table"),
            ui.download_button("dl_gsea_tsv", "Download full table (.tsv)", class_="btn-sm mt-2"),
        ),
        ui.accordion(
            ui.accordion_panel(
                "Which gene sets to show",
                ui.layout_columns(
                    ui.input_numeric("g_top", "Top N", 10, min=1, max=60),
                    ui.input_numeric("g_q", "FDR q cutoff", 0.05, min=0, max=1, step=0.01)),
                ui.input_radio_buttons("g_dir", "Direction",
                                       {"both": "Both", "positive": "Positive NES",
                                        "negative": "Negative NES"}, inline=True),
                ui.input_checkbox("g_idep", "Mark sets also found by iDEP (†)", True),
            ),
            ui.accordion_panel(
                "Design",
                ui.input_text("g_title", "Title", GseaDotParams.title),
                ui.input_text("g_sub", "Subtitle", GseaDotParams.subtitle),
                ui.input_text("g_foot", "Footnote (empty = automatic)", ""),
                ui.layout_columns(
                    ui.input_numeric("g_wrap", "Wrap names at", 50, min=15, max=120),
                    ui.input_numeric("g_qmax", "Colour max −log₁₀q", 4, min=2, max=20)),
                ui.input_select("g_cmap", "Colour map",
                                {"viridis": "Viridis (template)", "plasma": "Plasma",
                                 "cividis": "Cividis (colour-blind safe)", "Reds": "Reds"}),
                size_inputs("g", GseaDotParams.width, 0),
            ),
            ui.accordion_panel(
                "GSEA settings",
                ui.layout_columns(
                    ui.input_numeric("g_min", "Min set size", 15, min=1),
                    ui.input_numeric("g_max", "Max set size", 500, min=5)),
                ui.layout_columns(
                    ui.input_numeric("g_nperm", "Permutations", 1000, min=100, max=100000, step=100),
                    ui.input_numeric("g_seed", "Random seed", 42)),
                ui.p("Small p-values beyond the permutation floor are refined with "
                     "multilevel splitting, as in fgsea.", class_="ez-small"),
                ui.input_action_button("gsea_apply", "Apply and re-run", class_="btn-sm"),
            ),
            open="Which gene sets to show",
        ),
        col_widths=(8, 4),
    ),
)

export_tab = ui.nav_panel(
    "Export",
    ui.layout_columns(
        ui.card(
            ui.card_header("Figure files"),
            ui.input_checkbox_group("x_fmt", "Formats (for 'Download all')", export.FORMATS,
                                    selected=["png", "pdf"]),
            ui.input_select("x_one", "Format for single-plot downloads", export.FORMATS),
            ui.input_select("x_dpi", "Resolution (PNG/TIFF)", {"300": "300 dpi (journals)",
                            "600": "600 dpi", "1200": "1200 dpi (line art)", "150": "150 dpi (slides)"}),
            ui.input_select("x_width", "Printed width (PNG/TIFF)", list(export.WIDTHS_MM)),
            ui.p("PDF and SVG are vector files with editable text; they scale without loss.",
                 class_="ez-small"),
        ),
        ui.card(
            ui.card_header("Everything at once"),
            ui.p("A zip with every plot in the chosen formats, all GSEA tables that have been "
                 "run, and a parameters file so the figures can be reproduced."),
            ui.output_ui("x_status"),
            ui.input_checkbox("x_runall", "Also run GSEA for collections not run yet", False),
            ui.download_button("dl_all", "Download all (.zip)", class_="btn-primary"),
        ),
    ),
)

help_tab = ui.nav_panel(
    "Help",
    ui.markdown(f"""
### What to download from iDEP
| EZplot needs | Where in iDEP | Used for |
|---|---|---|
| **DE results matrix** (has columns like `CTRL-TRT_log2FC`, `CTRL-TRT_adjPval` and one column per sample) | *DEG* tab → download the full results table | Volcano, GSEA |
| **DEG_Heatmap_Data.csv** | *DEG* tab → heatmap → download data | Heatmap |
| `sig_pathways*.csv` *(optional)* | *Pathway* tab → download | Marks sets iDEP also found |
| a `.gmt` file *(optional)* | MSigDB, Enrichr, … | Extra gene-set collection |

Drop all files at once; EZplot recognises each one by its columns.

### Direction of the comparison
Pick your control and treatment groups in the sidebar; every plot is oriented so that
positive = higher in the treatment. iDEP's log2FC column name does not reliably say which
group its sign refers to, so EZplot reads that from the per-sample expression in the same
file. Only if the file has no sample columns does it ask you.

### GSEA and NES
iDEP's exported "NES" column actually holds the raw enrichment score (ES), and NES cannot be
recovered from that file. EZplot therefore re-runs preranked GSEA on the full log2FC ranking
against {genesets.VERSION} (Hallmark, GO BP/MF/CC; mouse and human), with the same statistics
as fgsea. NES agree with fgsea to about 0.01.

EZplot {__version__}
"""),
)

app_ui = ui.page_sidebar(
    sidebar,
    ui.navset_card_underline(volcano_tab, heatmap_tab, gsea_tab, export_tab, help_tab,
                             id="tabs"),
    ui.head_content(ui.tags.style(CSS), ui.tags.script(JS)),
    title="EZplot",
    fillable=False,
)


# ── server ────────────────────────────────────────────────────

def server(input, output, session):
    proj = Project()
    rev = reactive.value(0)               # bumped whenever proj changes
    struct = reactive.value(0)            # bumped when files or sample groups change
    files_state = reactive.value([])      # [(name, kind, error)]

    def bump(structure=False):
        with reactive.isolate():
            rev.set(rev() + 1)
            if structure:
                struct.set(struct() + 1)

    def load_files(items):
        with reactive.isolate():
            out = list(files_state())
        for path, name in items:
            try:
                kind = proj.load(path, name)
                err = None if kind not in ("unknown", "deg_list") else io.KIND_LABELS[kind]
            except io.InputError as e:
                kind, err = "unknown", str(e)
            out = [x for x in out if x[0] != name] + [(name, kind, err)]
        files_state.set(out)
        bump(structure=True)

    @reactive.effect
    @reactive.event(input.files)
    def _():
        load_files([(f["datapath"], f["name"]) for f in input.files()])

    @reactive.effect
    @reactive.event(input.example)
    def _():
        load_files([(EXAMPLE / f, f) for f in
                    ("example_DE_matrix.csv", "example_DEG_Heatmap_Data.csv")])

    @render.ui
    def loaded():
        items = files_state()
        if not items:
            return ui.p("No files yet.", class_="ez-small")
        return ui.div(*[ui.div(ui.span(io.KIND_LABELS.get(k, k) if not e else e,
                                       class_="ez-badge " + ("warn" if e else "ok")),
                               ui.span(n, class_="ez-small"))
                        for n, k, e in items])

    # ── experiment card ──
    @render.ui
    def experiment():
        struct()
        groups = proj.group_names
        if not groups:
            return ui.p("Load a DE matrix or heatmap file first.", class_="ez-small")
        members = {g: [s for s, gg in proj.groups.items() if gg == g] for g in groups}
        sentence = []
        if proj.de is not None:
            sentence = [ui.output_ui("direction")]
        with reactive.isolate():
            return ui.TagList(
                ui.layout_columns(
                    ui.input_select("control", "Control", groups, selected=proj.control),
                    ui.input_select("treatment", "Treatment", groups, selected=proj.treatment)),
                ui.p(*[ui.span(f"{g}: ", ui.tags.b(len(m)), f" ({', '.join(m)})", ui.br())
                       for g, m in members.items()], class_="ez-small"),
                *sentence,
                ui.accordion(ui.accordion_panel(
                    "Labels, colours, species",
                    ui.layout_columns(
                        ui.input_text("control_label", "Control label", proj.labels[0]),
                        ui.input_text("treatment_label", "Treatment label", proj.labels[1])),
                    ui.div(color_input("control_color", "Control", proj.control_color),
                           color_input("treatment_color", "Treatment", proj.treatment_color),
                           class_="mb-2"),
                    ui.input_select("species", "Gene sets for species", genesets.SPECIES,
                                    selected=proj.species),
                    ui.input_text_area("group_edit", "Fix sample groups (one 'sample = group' per line)",
                                       rows=3, placeholder="WT5R = WT"),
                    ui.input_action_button("group_apply", "Apply groups", class_="btn-sm"),
                ), open=False),
            )

    @render.ui
    def direction():
        rev()
        out = [ui.p(proj.direction_sentence(), class_="ez-dir mt-2 mb-1")]
        if proj.de.expr.empty or proj.direction_confidence < 0.6:
            # the file alone can't say which group iDEP's log2FC sign refers to
            with reactive.isolate():
                out.append(ui.input_select(
                    "up_group", "In this iDEP file, a positive log₂FC means higher in:",
                    proj.group_names, selected=proj.up_group))
        return ui.TagList(*out)

    @reactive.effect
    @reactive.event(input.control, input.treatment, ignore_init=True)
    def _():
        c, t = input.control(), input.treatment()
        if c == t:
            others = [g for g in proj.group_names if g != c]
            if others:
                t = others[0] if c != proj.control else proj.control
                ui.update_select("treatment" if c != proj.control else "control",
                                 selected=t)
                return
        if (c, t) != (proj.control, proj.treatment):
            old_labels = proj.labels
            proj.control, proj.treatment = c, t
            proj.control_label = proj.treatment_label = None
            if old_labels != proj.labels:
                ui.update_text("control_label", value=proj.labels[0])
                ui.update_text("treatment_label", value=proj.labels[1])
            bump()

    @reactive.effect
    @reactive.event(input.up_group, ignore_init=True)
    def _():
        if input.up_group() in proj.group_names and input.up_group() != proj.up_group:
            proj.up_group = input.up_group()
            bump()

    @reactive.effect
    def _():
        try:
            cl, tl = input.control_label(), input.treatment_label()
            cc, tc = input.control_color(), input.treatment_color()
            sp = input.species()
        except Exception:      # noqa: BLE001 - inputs not rendered yet
            return
        changed = False
        for attr, val in (("control_label", cl), ("treatment_label", tl),
                          ("control_color", cc), ("treatment_color", tc), ("species", sp)):
            if val and getattr(proj, attr) != val:
                setattr(proj, attr, val)
                changed = True
        if changed:
            bump()

    @reactive.effect
    @reactive.event(input.group_apply)
    def _():
        for line in input.group_edit().splitlines():
            if "=" in line:
                s, g = (x.strip() for x in line.split("=", 1))
                if s in proj.groups and g:
                    proj.groups[s] = g
        names = proj.group_names
        if proj.control not in names:
            proj.control = names[0]
        if proj.treatment not in names or proj.treatment == proj.control:
            proj.treatment = next((g for g in names if g != proj.control), None)
        proj.detect_direction()
        bump(structure=True)

    # ── volcano ──
    def volcano_params():
        genes = [g for chunk in input.v_genes().replace(",", "\n").splitlines()
                 for g in [chunk.strip()] if g]
        return VolcanoParams(
            title=input.v_title(), fdr=input.v_fdr() or 0.05, lfc=input.v_lfc() or 0,
            label_mode=input.v_mode(), label_top=input.v_top() or 20, label_genes=genes,
            neg_color=input.v_neg() or style.NEG, pos_color=input.v_pos() or style.POS,
            ns_color=input.v_ns() or style.NS, xlim=input.v_xlim() or 0,
            ymax=input.v_ymax() or 0, point_size=input.v_ptsize() or 70,
            label_fontsize=input.v_lfs() or 13, width=input.v_w() or 8.1,
            height=input.v_h() or 7.0)

    @render.ui
    def volcano_out():
        rev()
        if proj.de is None:
            return empty("**Load the iDEP DE results matrix** to draw the volcano plot.")
        return fig_img(proj.volcano(volcano_params()))

    # ── heatmap ──
    def heatmap_params():
        return HeatmapParams(
            title=input.h_title(), zscore=input.h_z(), vlim=input.h_vlim() or 2,
            cmap=_cmap(input.h_cmap()), cluster_rows=input.h_crow(), cluster_cols=input.h_ccol(),
            show_genes=input.h_genes(), italic_genes=input.h_ital(),
            gene_fontsize=input.h_gfs() or 15, sample_fontsize=input.h_sfs() or 14,
            width=input.h_w() or 10.5, height=input.h_h() or 0, legend_title=input.h_legend())

    @render.ui
    def heatmap_out():
        rev()
        if proj.heatmap is None:
            return empty("**Load DEG_Heatmap_Data.csv** from iDEP to draw the heatmap.")
        return fig_img(proj.heatmap_fig(heatmap_params()))

    # ── GSEA ──
    gsea_rev = reactive.value(0)

    @reactive.effect
    def _():
        rev()
        colls = proj.collections()
        with reactive.isolate():
            cur = input.g_coll()
        choices = {k: v[0].replace("GSEA ", "") + (" ✓" if proj.gsea_done(k) else "")
                   for k, v in colls.items()}
        ui.update_select("g_coll", choices=choices,
                         selected=cur if cur in colls else next(iter(colls), None))

    def gsea_settings():
        return GseaSettings(min_size=input.g_min() or 15, max_size=input.g_max() or 500,
                            nperm=input.g_nperm() or 1000, seed=input.g_seed() or 0)

    def run_collection(key):
        with ui.Progress(min=0, max=1) as p:
            p.set(0, message=proj.collections()[key][0], detail="starting")
            proj.run_gsea(key, progress=lambda f, m: p.set(f, detail=m))

    @reactive.effect
    def _():
        # run (or reuse the cached result) whenever the GSEA tab shows a collection
        rev()
        if input.tabs() != "GSEA" or proj.de is None or not input.g_coll():
            return
        key = input.g_coll()
        if key not in proj.collections() or proj.gsea_done(key):
            return
        run_collection(key)
        gsea_rev.set(gsea_rev() + 1)
        bump()

    @reactive.effect
    @reactive.event(input.gsea_apply)
    def _():
        proj.gsea = gsea_settings()
        bump()

    def gsea_params():
        return GseaDotParams(
            title=input.g_title(), subtitle=input.g_sub(), footnote=input.g_foot(),
            q=input.g_q() or 0.05, top=input.g_top() or 10, qmax=input.g_qmax() or 4,
            wrap=input.g_wrap() or 50, cmap=input.g_cmap(), width=input.g_w() or 10,
            height=input.g_h() or 0, direction=input.g_dir())

    def current_gsea_ready():
        rev(), gsea_rev()
        if proj.de is None:
            return None
        key = input.g_coll()
        return key if key and proj.gsea_done(key) else None

    @render.ui
    def gsea_out():
        if proj.de is None:
            return empty("**Load the iDEP DE results matrix**; GSEA ranks all genes by log₂FC.")
        key = current_gsea_ready()
        if key is None:
            return empty("Running GSEA… (Hallmarks: a few seconds; GO Biological Process: "
                         "up to a minute). Results are kept, so switching back is instant.")
        saved = proj.idep_pathways
        if not input.g_idep():
            proj.idep_pathways = None
        try:
            return fig_img(proj.gsea_fig(key, gsea_params()))
        finally:
            proj.idep_pathways = saved

    @render.data_frame
    def gsea_table():
        key = current_gsea_ready()
        if key is None:
            return None
        df = proj.run_gsea(key).copy()
        df["pathway"] = [style.pretty_name(p, None) for p in df["pathway"]]
        for c in ("ES", "NES"):
            df[c] = df[c].round(3)
        for c in ("pval", "padj"):
            df[c] = df[c].map(lambda v: f"{v:.2e}")
        return render.DataGrid(df.drop(columns="leadingEdge"), height="320px", filters=True)

    @render.ui
    def x_status():
        rev(), gsea_rev()
        if proj.de is None:
            return ui.p("GSEA: load the DE matrix first.", class_="ez-small")
        colls = proj.collections()
        done = [v[0].replace("GSEA ", "") for k, v in colls.items() if proj.gsea_done(k)]
        todo = [v[0].replace("GSEA ", "") for k, v in colls.items() if not proj.gsea_done(k)]
        return ui.p(ui.tags.b("GSEA included: "), ", ".join(done) or "none yet", ui.br(),
                    ui.tags.b("Not run with the current settings: "), ", ".join(todo) or "none",
                    class_="ez-small")

    # ── downloads ──
    def width_mm():
        return export.WIDTHS_MM.get(input.x_width())

    def one(fig_fn, stem):
        fmt = input.x_one()
        fig = fig_fn()
        try:
            data = export.figure_bytes(fig, fmt, int(input.x_dpi()), width_mm())
        finally:
            plt.close(fig)
        return data

    @render.download_button(filename=lambda: f"volcano.{input.x_one()}")
    def dl_volcano():
        if proj.de is not None:
            yield one(lambda: proj.volcano(volcano_params()), "volcano")

    @render.download_button(filename=lambda: f"heatmap.{input.x_one()}")
    def dl_heatmap():
        if proj.heatmap is not None:
            yield one(lambda: proj.heatmap_fig(heatmap_params()), "heatmap")

    @render.download_button(filename=lambda: f"gsea_{input.g_coll()}.{input.x_one()}")
    def dl_gsea():
        key = current_gsea_ready()
        if key:
            yield one(lambda: proj.gsea_fig(key, gsea_params()), f"gsea_{key}")

    @render.download_button(filename=lambda: f"gsea_{input.g_coll()}.tsv")
    def dl_gsea_tsv():
        key = current_gsea_ready()
        if key:
            yield proj.run_gsea(key).to_csv(sep="\t", index=False)

    @render.download_button(filename="ezplot_figures.zip")
    def dl_all():
        figs, tables = {}, {}
        if proj.de is not None:
            vp = volcano_params()
            figs["volcano"] = lambda: proj.volcano(vp)
        if proj.heatmap is not None:
            hp = heatmap_params()
            figs["heatmap"] = lambda: proj.heatmap_fig(hp)
        gp = gsea_params()
        if proj.de is not None:
            for key in proj.collections():
                if not proj.gsea_done(key) and input.x_runall():
                    run_collection(key)
                if proj.gsea_done(key):
                    figs[f"gsea_{key}"] = (lambda k=key: proj.gsea_fig(k, gp))
                    tables[f"gsea_{key}"] = proj.run_gsea(key)
        params = {
            "files": [n for n, _, _ in files_state()],
            "control": proj.control, "treatment": proj.treatment,
            "direction": proj.direction_sentence(), "species": proj.species,
            "volcano": asdict(volcano_params()) if proj.de is not None else None,
            "heatmap": asdict(heatmap_params()) if proj.heatmap is not None else None,
            "gsea_plot": asdict(gp), "gsea_settings": asdict(proj.gsea),
            "gene_sets": genesets.VERSION, "ezplot": __version__,
        }
        yield export.bundle(figs, tables, params, input.x_fmt() or ["png"],
                            int(input.x_dpi()), width_mm())


def _cmap(name):
    if name == "vlag":
        from matplotlib.colors import LinearSegmentedColormap
        return LinearSegmentedColormap.from_list("vlag", ["#2369bd", "#f7f7f7", "#a9373b"])
    return name


app = App(app_ui, server)
