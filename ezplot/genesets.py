"""Bundled MSigDB v2026.1 collections (CC BY 4.0, https://www.gsea-msigdb.org)."""

from functools import lru_cache
from pathlib import Path

from .io import read_gmt

DIR = Path(__file__).parent / "genesets"
_V = "v2026.1"
VERSION = f"MSigDB {_V}"


def _files(mouse, human):
    out = {}
    if mouse:
        out["mouse"] = f"{mouse}.{_V}.Mm.symbols.gmt.gz"
    if human:
        out["human"] = f"{human}.{_V}.Hs.symbols.gmt.gz"
    return out


# key -> (title, y-axis label, {species: file}); MSigDB has no KEGG for mouse
COLLECTIONS = {
    "hallmarks": ("GSEA Hallmarks", "Hallmark gene set", _files("mh.all", "h.all")),
    "go_bp": ("GSEA GO Biological Process", "GO biological process", _files("m5.go.bp", "c5.go.bp")),
    "go_mf": ("GSEA GO Molecular Function", "GO molecular function", _files("m5.go.mf", "c5.go.mf")),
    "go_cc": ("GSEA GO Cellular Component", "GO cellular component", _files("m5.go.cc", "c5.go.cc")),
    "reactome": ("GSEA Reactome", "Reactome pathway", _files("m2.cp.reactome", "c2.cp.reactome")),
    "kegg": ("GSEA KEGG MEDICUS", "KEGG pathway", _files(None, "c2.cp.kegg_medicus")),
}
SPECIES = ["mouse", "human"]


def available(species):
    """Collection keys bundled for `species`, in menu order."""
    return [k for k, v in COLLECTIONS.items() if species in v[2]]


@lru_cache(maxsize=None)
def load(key, species):
    return read_gmt(DIR / COLLECTIONS[key][2][species])
