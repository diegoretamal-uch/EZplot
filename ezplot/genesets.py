"""Bundled MSigDB v2024.1 collections (CC BY 4.0, https://www.gsea-msigdb.org)."""

from functools import lru_cache
from pathlib import Path

from .io import read_gmt

DIR = Path(__file__).parent / "genesets"
VERSION = "MSigDB v2024.1"

# key -> (title, y-axis label, {species: file})
COLLECTIONS = {
    "hallmarks": ("GSEA Hallmarks", "Hallmark gene set",
                  {"mouse": "mh.all.v2024.1.Mm.symbols.gmt.gz", "human": "h.all.v2024.1.Hs.symbols.gmt.gz"}),
    "go_bp": ("GSEA GO Biological Process", "GO biological process",
              {"mouse": "m5.go.bp.v2024.1.Mm.symbols.gmt.gz", "human": "c5.go.bp.v2024.1.Hs.symbols.gmt.gz"}),
    "go_mf": ("GSEA GO Molecular Function", "GO molecular function",
              {"mouse": "m5.go.mf.v2024.1.Mm.symbols.gmt.gz", "human": "c5.go.mf.v2024.1.Hs.symbols.gmt.gz"}),
    "go_cc": ("GSEA GO Cellular Component", "GO cellular component",
              {"mouse": "m5.go.cc.v2024.1.Mm.symbols.gmt.gz", "human": "c5.go.cc.v2024.1.Hs.symbols.gmt.gz"}),
}
SPECIES = ["mouse", "human"]


@lru_cache(maxsize=None)
def load(key, species):
    return read_gmt(DIR / COLLECTIONS[key][2][species])
