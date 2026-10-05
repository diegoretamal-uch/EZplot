"""Everything one analysis needs, shared by the Shiny app and the command line."""

import hashlib
from dataclasses import asdict, dataclass, field

import pandas as pd

from . import genesets, io, style
from .gsea.engine import prerank
from .plots.gsea_dot import GseaDotParams, plot_gsea_dot
from .plots.heatmap import HeatmapParams, plot_heatmap
from .plots.volcano import VolcanoParams, plot_volcano


@dataclass
class GseaSettings:
    min_size: int = 15
    max_size: int = 500
    nperm: int = 1000
    seed: int = 42


@dataclass
class Project:
    de: io.DEMatrix = None
    heatmap: pd.DataFrame = None
    idep_pathways: pd.DataFrame = None
    custom_sets: dict = None
    custom_name: str = "Custom gene sets"
    groups: dict = field(default_factory=dict)       # sample -> group (DE matrix + heatmap)
    control: str = None
    treatment: str = None
    up_group: str = None                             # group higher when raw log2FC > 0
    direction_confidence: float = 0.0
    species: str = "mouse"
    control_label: str = None                        # legend/subtitle text
    treatment_label: str = None
    control_color: str = style.CONTROL
    treatment_color: str = style.TREATMENT
    gsea: GseaSettings = field(default_factory=GseaSettings)
    _gsea_cache: dict = field(default_factory=dict, repr=False)

    # ── loading ───────────────────────────────────────────────
    def load(self, src, name):
        """Load one file; returns its detected kind."""
        if name.lower().endswith(".gmt"):
            self.custom_sets = io.read_gmt(src)
            self.custom_name = name.rsplit(".", 1)[0]
            self._gsea_cache.pop("custom", None)
            return "gmt"
        df = io.read_table(src, name)
        kind = io.detect_kind(df, name)
        if kind == "de_matrix":
            self.set_de(io.parse_de_matrix(df))
        elif kind == "heatmap":
            self.heatmap = io.parse_heatmap(df)
            self._refresh_groups()
        elif kind == "pathways":
            self.idep_pathways = io.parse_pathways(df)
        return kind

    def set_de(self, de):
        self.de = de
        self._gsea_cache.clear()
        if de.species in genesets.SPECIES:
            self.species = de.species
        self._refresh_groups()

    def _refresh_groups(self):
        samples = []
        hints = ()
        if self.de is not None:
            samples += list(self.de.expr.columns)
            hints = self.de.contrast
        if self.heatmap is not None:
            samples += [s for s in self.heatmap.columns if s not in samples]
        self.groups = io.assign_groups(samples, hints=hints)
        names = self.group_names
        if self.control not in names or self.treatment not in names:
            preferred = self.de.contrast[0] if self.de is not None else None
            self.control = io.guess_control(names, preferred)
            self.treatment = next((g for g in names if g != self.control), None)
        self.detect_direction()

    def detect_direction(self):
        if self.de is None:
            return
        up, conf = io.detect_direction(self.de, self.groups)
        if up is None:
            # no expression columns to check against: iDEP's "A-B_log2FC" was positive
            # for genes higher in A in every file we verified
            up = self.de.contrast[0]
        self.up_group, self.direction_confidence = up, conf

    @property
    def group_names(self):
        return list(dict.fromkeys(self.groups.values()))

    @property
    def labels(self):
        return (self.control_label or self.control or "control",
                self.treatment_label or self.treatment or "treatment")

    # ── derived data ──────────────────────────────────────────
    def oriented_genes(self):
        """DE genes with `lfc` positive = higher in treatment."""
        if self.de is None:
            return None
        sign = 1.0 if self.up_group == self.treatment else -1.0
        return self.de.genes.assign(lfc=self.de.genes["raw_lfc"] * sign)

    def direction_sentence(self):
        if self.de is None:
            return ""
        c, t = self.labels
        return f"Positive log₂FC / NES = higher in {t} than in {c}"

    def ranking(self):
        g = self.oriented_genes()
        return pd.Series(g["lfc"].to_numpy(), index=g["symbol"].to_numpy())

    def collections(self):
        """Available GSEA collections: key -> (title, y label)."""
        out = {}
        if self.species in genesets.SPECIES:
            out.update({k: v[:2] for k, v in genesets.COLLECTIONS.items()})
        if self.custom_sets:
            out["custom"] = (f"GSEA {self.custom_name}", "Gene set")
        return out

    def gene_sets(self, key):
        if key == "custom":
            return self.custom_sets
        return genesets.load(key, self.species)

    def _cache_key(self, key):
        r = self.ranking()
        h = hashlib.sha1(pd.util.hash_pandas_object(r, index=True).values.tobytes())
        h.update(repr((key, self.species, asdict(self.gsea))).encode())
        return h.hexdigest()

    def run_gsea(self, key, progress=None):
        ck = self._cache_key(key)
        if ck not in self._gsea_cache:
            s = self.gsea
            self._gsea_cache[ck] = prerank(self.ranking(), self.gene_sets(key), s.min_size,
                                           s.max_size, s.nperm, s.seed, progress=progress)
        return self._gsea_cache[ck]

    def gsea_done(self, key):
        return self.de is not None and self._cache_key(key) in self._gsea_cache

    def gsea_footnote(self, key):
        src = genesets.VERSION if key != "custom" else self.custom_name
        s = self.gsea
        return (f"Preranked GSEA on DESeq2 log$_2$FC from iDEP ({self.labels[1]} vs "
                f"{self.labels[0]}) · {src} · set size {s.min_size}–{s.max_size}")

    # ── figures ───────────────────────────────────────────────
    def volcano(self, p=None):
        p = p or VolcanoParams()
        return plot_volcano(self.oriented_genes(), p)

    def heatmap_fig(self, p=None):
        c, t = self.labels
        roles = {s: ("control" if g == self.control else "treatment" if g == self.treatment else g)
                 for s, g in self.groups.items()}
        labels = {"control": c, "treatment": t}
        colors = {"control": self.control_color, "treatment": self.treatment_color}
        return plot_heatmap(self.heatmap, roles, labels, colors, p or HeatmapParams())

    def gsea_fig(self, key, p=None):
        res = self.run_gsea(key)
        title, ylabel = self.collections()[key]
        p = p or GseaDotParams()
        if not p.footnote:
            p = GseaDotParams(**{**asdict(p), "footnote": self.gsea_footnote(key)})
        idep = set(self.idep_pathways["key"]) if self.idep_pathways is not None else None
        c, t = self.labels
        return plot_gsea_dot(res, title, ylabel, c, t, p, idep_keys=idep,
                             key_func=io.pathway_key)
