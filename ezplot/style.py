"""Shared look: bundled font, palettes and pathway-name humanisation.

The font (Liberation Sans, SIL OFL 1.1, Arial-metric) ships with the package so figures
render identically on Windows, macOS, Linux and in the browser build.
"""

import re
import textwrap
from pathlib import Path

import matplotlib
from matplotlib import font_manager

FONT_DIR = Path(__file__).parent / "fonts"
FONT = "Liberation Sans"

# template colours (sampled from the reference figures EZplot reproduces)
CONTROL = "#18c7c9"
TREATMENT = "#f28b82"
NEG = "#0072B2"
POS = "#D55E00"
NS = "#c9ccd1"
ZERO_LINE = "#1f77d4"
TEXT = "#222222"
MUTED = "#666666"

_registered = False


def setup():
    """Register the bundled font and set the rcParams every plot relies on."""
    global _registered
    if not _registered:
        for f in FONT_DIR.glob("*.ttf"):
            font_manager.fontManager.addfont(str(f))
        _registered = True
    matplotlib.rcParams.update({
        "font.family": [FONT, "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT,
        "mathtext.it": f"{FONT}:italic",
        "mathtext.bf": f"{FONT}:bold",
        "axes.unicode_minus": True,
        "svg.fonttype": "none",   # editable text in Illustrator/Inkscape
        "pdf.fonttype": 42,       # TrueType embedded, journals accept it
        "savefig.facecolor": "white",
    })


# words that keep their capitalisation when MSigDB names are humanised
WORDS = {
    "mhc": "MHC", "ii": "II", "i": "I", "atp": "ATP", "adp": "ADP", "dna": "DNA",
    "rna": "RNA", "mrna": "mRNA", "rrna": "rRNA", "trna": "tRNA", "ncrna": "ncRNA",
    "nadh": "NADH", "nadph": "NADPH", "ecm": "ECM", "er": "ER", "camp": "cAMP",
    "cgmp": "cGMP", "gtpase": "GTPase", "atpase": "ATPase", "mapk": "MAPK",
    "erk1": "ERK1", "erk2": "ERK2", "jak": "JAK", "stat": "STAT", "stat3": "STAT3",
    "stat5": "STAT5", "tgf": "TGF", "kras": "KRAS", "p53": "p53", "e2f": "E2F",
    "g2m": "G2/M", "mtorc1": "mTORC1", "pi3k": "PI3K", "akt": "AKT", "mtor": "mTOR",
    "wnt": "Wnt", "il2": "IL-2", "il6": "IL-6", "myc": "Myc", "uv": "UV",
    "nfkb": "NF-κB", "tnfa": "TNFα", "cd4": "CD4", "cd8": "CD8", "nk": "NK",
    "iv": "IV", "iii": "III", "gpi": "GPI", "snrna": "snRNA", "snorna": "snoRNA",
    "t": "T", "b": "B", "g": "G", "rho": "Rho", "ras": "Ras", "ca2": "Ca²⁺",
    "coa": "CoA", "ap": "AP", "u2": "U2",
}
PHRASES = [
    ("interferon alpha", "interferon-α"), ("interferon gamma", "interferon-γ"),
    ("epithelial mesenchymal", "epithelial–mesenchymal"), ("TGF beta", "TGF-β"),
    (" up$", " (up)"), (" dn$", " (down)"), (" early$", " (early)"),
    (" late$", " (late)"), ("^Il ", "IL-"),
]


def pretty_name(name, width=44, max_lines=2):
    """HALLMARK_TNFA_SIGNALING_VIA_NFKB -> 'TNFα signaling via NF-κB' (wrapped)."""
    name = re.sub(r"^(HALLMARK|GOBP|GOMF|GOCC|KEGG|REACTOME|WP)[_ ]", "", name.strip())
    out = []
    for w in re.split(r"[_ ]+", name.lower()):
        parts = w.split("-")
        out.append("-".join(WORDS.get(p, p) for p in parts))
    s = " ".join(out)
    for a, b in PHRASES:
        s = re.sub(a, b, s)
    s = s[:1].upper() + s[1:]
    if width is None:
        return s
    lines = textwrap.wrap(s, width)
    if len(lines) > max_lines:
        rest = " ".join(lines[max_lines - 1:])
        lines = lines[:max_lines - 1] + [textwrap.shorten(rest, width, placeholder="…")]
    return "\n".join(lines)

