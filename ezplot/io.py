"""Readers for iDEP exports, file-type sniffing, sample groups and contrast direction.

iDEP files handled:
  de_matrix  - the DEG tab's full results download: symbol, ensembl_ID, ..., "<A>-<B>_log2FC",
               "<A>-<B>_adjPval", then per-sample processed (log) expression
  heatmap    - DEG_Heatmap_Data.csv: rows = genes, columns = samples, row-centred values
  pathways   - sig_pathways*.csv from the Pathway tab (only used as a cross-check)
  gmt        - user gene-set collection

The sign of iDEP's log2FC column is not reliably given by its name ("A-B_log2FC" can be
negative for genes that are higher in B), so the direction is
inferred from the expression columns of the same file; see detect_direction().
"""

import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

ENSEMBL_SUFFIX = re.compile(r"\s*ENS[A-Z]*G\d{6,}(\.\d+)?$")
CONTROL_WORDS = re.compile(
    r"ctrl|control|ctl|wt|wild|vehicle|veh|sham|mock|untreated|naive|baseline|dmso|scr",
    re.I)
SPECIES_PREFIX = [("ENSMUSG", "mouse"), ("ENSRNOG", "rat"), ("ENSDARG", "zebrafish"),
                  ("ENSG", "human")]


class InputError(ValueError):
    """A file could not be understood; the message is shown to the user as-is."""


def read_table(src, name=None):
    """Read a csv/tsv/txt from a path, bytes or file-like; separator is sniffed."""
    if isinstance(src, (str, Path)):
        name = name or str(src)
        raw = Path(src).read_bytes()
    elif isinstance(src, bytes):
        raw = src
    else:
        raw = src.read()
    text = raw.decode("utf-8-sig", errors="replace")
    first = text.split("\n", 1)[0]
    sep = "\t" if first.count("\t") > first.count(",") else ","
    try:
        return pd.read_csv(io.StringIO(text), sep=sep, low_memory=False)
    except Exception as e:  # noqa: BLE001 - surface any parser error to the user
        raise InputError(f"Could not read {name or 'file'} as a table: {e}") from e


def clean_symbol(s):
    """iDEP glues the Ensembl ID to ambiguous symbols ('Gm31520ENSMUSG00000121703')."""
    if not isinstance(s, str):
        return s
    out = ENSEMBL_SUFFIX.sub("", s).strip()
    return out or s.strip()


def contrast_columns(df):
    lfc = [c for c in df.columns if str(c).endswith("_log2FC")]
    padj = [c for c in df.columns if str(c).endswith("_adjPval")]
    return lfc, padj


def detect_kind(df, name=""):
    cols = {str(c) for c in df.columns}
    lfc, padj = contrast_columns(df)
    if lfc and padj:
        return "de_matrix"
    if "adj.Pval" in cols and ("NES" in cols or "Direction" in cols):
        return "pathways"
    if {"log2FC", "Adjusted_P_value"} <= cols:
        return "deg_list"
    if cols & {"pval", "padj", "NES", "pvalue", "log2FoldChange"}:
        return "unknown"
    numeric = df.select_dtypes("number").columns
    if len(numeric) >= 2 and len(df) <= 5000 and df.shape[1] - len(numeric) <= 2:
        return "heatmap"
    if len(df.columns) <= 2 and df.empty and "pathway" in name.lower():
        return "pathways"
    return "unknown"


KIND_LABELS = {
    "de_matrix": "DE matrix (volcano + GSEA)",
    "heatmap": "Heatmap data",
    "pathways": "iDEP pathways (cross-check)",
    "deg_list": "DEG gene list (not needed: load the full DE matrix instead)",
    "gmt": "Custom gene sets (.gmt)",
    "unknown": "Unrecognised file",
}


# ── sample groups ─────────────────────────────────────────────

def _stem(sample):
    """'KO5R' -> 'KO', 'KO_3' -> 'KO', 'WT-rep2' -> 'WT'."""
    s = re.sub(r"[_.\- ]?(rep|r|s|n)?\d+[A-Za-z]?$", "", str(sample), flags=re.I)
    return s or str(sample)


def assign_groups(samples, hints=()):
    """Map each sample to a group, preferring the longest matching hint prefix."""
    hints = sorted({h for h in hints if h}, key=len, reverse=True)
    out = {}
    for s in samples:
        match = next((h for h in hints if str(s).upper().startswith(h.upper())), None)
        out[s] = match or _stem(s)
    return out


def guess_control(groups, preferred=None):
    """Pick the control group: name keywords first, then the caller's preference."""
    groups = list(dict.fromkeys(groups))
    for g in groups:
        if CONTROL_WORDS.search(g):
            return g
    if preferred in groups:
        return preferred
    return groups[0] if groups else None


# ── DE matrix ─────────────────────────────────────────────────

@dataclass
class DEMatrix:
    genes: pd.DataFrame            # symbol, ensembl, raw_lfc, padj, baseMean
    expr: pd.DataFrame             # genes x samples, processed expression (may be empty)
    contrast: tuple                # (A, B) as written in "<A>-<B>_log2FC"
    lfc_col: str
    groups: dict = field(default_factory=dict)   # sample -> group
    species: str = "unknown"

    @property
    def group_names(self):
        return list(dict.fromkeys(self.groups.values()))


def parse_de_matrix(df, contrast=None):
    lfc_cols, padj_cols = contrast_columns(df)
    if not lfc_cols:
        raise InputError("No '<A>-<B>_log2FC' column found; is this iDEP's DE results file?")
    lfc_col = lfc_cols[0]
    if contrast is not None:
        lfc_col = next(c for c in lfc_cols if c.startswith(contrast))
    stem = lfc_col[: -len("_log2FC")]
    padj_col = f"{stem}_adjPval" if f"{stem}_adjPval" in df.columns else padj_cols[0]
    parts = stem.split("-")
    contrast_pair = (parts[0], "-".join(parts[1:])) if len(parts) >= 2 else (stem, "")

    sym_col = next((c for c in df.columns if str(c).lower() == "symbol"), df.columns[0])
    ens_col = next((c for c in df.columns if "ensembl" in str(c).lower()), None)
    genes = pd.DataFrame({
        "symbol": df[sym_col].astype("string").fillna("").map(clean_symbol),
        "ensembl": df[ens_col].astype("string") if ens_col else pd.NA,
        "raw_lfc": pd.to_numeric(df[lfc_col], errors="coerce"),
        "padj": pd.to_numeric(df[padj_col], errors="coerce"),
        "baseMean": pd.to_numeric(df["baseMean"], errors="coerce")
        if "baseMean" in df.columns else np.nan,
    })
    empty = genes["symbol"].str.strip() == ""
    if ens_col:
        genes.loc[empty, "symbol"] = genes.loc[empty, "ensembl"]

    # expression columns: numeric columns after the iDEP marker / search_label column
    known = {sym_col, ens_col, lfc_col, padj_col, "User_ID", "baseMean"} | set(lfc_cols) | set(padj_cols)
    start = 0
    for marker in ("search_label", "Processed data:"):
        if marker in df.columns:
            start = max(start, list(df.columns).index(marker) + 1)
    cand = [c for c in df.columns[start:] if c not in known]
    expr = df[cand].apply(pd.to_numeric, errors="coerce")
    expr = expr.loc[:, expr.notna().any()]   # genes iDEP filtered out are NA in every sample

    groups = assign_groups(expr.columns, hints=contrast_pair)
    return DEMatrix(genes=genes, expr=expr, contrast=contrast_pair, lfc_col=lfc_col,
                    groups=groups, species=detect_species(genes["ensembl"]))


def detect_species(ids):
    ids = pd.Series(ids).dropna().astype(str).head(2000)
    for prefix, sp in SPECIES_PREFIX:
        if len(ids) and ids.str.startswith(prefix).mean() > 0.5:
            return sp
    return "unknown"


def detect_direction(de, groups=None):
    """Return (group_up_when_raw_lfc_positive, confidence in [0, 1]) or (None, 0).

    Correlates iDEP's log2FC with the between-group difference of mean processed
    expression in the same file. Uses the most significant genes, where the signal is
    unambiguous.
    """
    groups = groups or de.groups
    names = list(dict.fromkeys(groups.values()))
    if de.expr.empty or len(names) != 2:
        return None, 0.0
    g1, g2 = names
    s1 = [s for s, g in groups.items() if g == g1 and s in de.expr.columns]
    s2 = [s for s, g in groups.items() if g == g2 and s in de.expr.columns]
    if not s1 or not s2:
        return None, 0.0
    diff = de.expr[s1].mean(axis=1) - de.expr[s2].mean(axis=1)   # g1 - g2
    d = pd.DataFrame({"lfc": de.genes["raw_lfc"], "diff": diff, "padj": de.genes["padj"]}).dropna()
    d = d[d["lfc"].abs() > 0]
    if len(d) < 10:
        return None, 0.0
    d = d.nsmallest(min(500, len(d)), "padj")
    agree = np.mean(np.sign(d["lfc"]) == np.sign(d["diff"]))
    up = g1 if agree >= 0.5 else g2
    return up, float(abs(agree - 0.5) * 2)


def oriented_lfc(de, treatment, up_group):
    """log2FC as treatment vs control: positive = higher in `treatment`."""
    sign = 1.0 if up_group == treatment else -1.0
    return de.genes["raw_lfc"] * sign


# ── heatmap ───────────────────────────────────────────────────

def parse_heatmap(df):
    first = df.columns[0]
    if not pd.api.types.is_numeric_dtype(df[first]):
        df = df.set_index(first)
    df = df.drop(columns=[c for c in df.columns if str(c).lower() in ("gene_id", "ensembl_id")],
                 errors="ignore")
    x = df.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all").dropna(axis=0, how="any")
    if x.shape[1] < 2 or x.empty:
        raise InputError("Heatmap file needs gene rows and at least two numeric sample columns.")
    x.index = [clean_symbol(str(i)) for i in x.index]
    return x


# ── iDEP pathway table (cross-check only) ─────────────────────

def pathway_key(name):
    """Normalise MSigDB ('GOBP_FOO_BAR') and iDEP ('Foo bar ') names to one key."""
    s = re.sub(r"^(HALLMARK|GOBP|GOMF|GOCC)[_ ]", "", str(name).strip().upper())
    return re.sub(r"[^A-Z0-9]+", "_", s).strip("_")


def parse_pathways(df):
    if "adj.Pval" not in df.columns:
        return pd.DataFrame(columns=["pathway", "key", "adj.Pval"])
    name_col = next((c for c in df.columns if "pathway" in str(c).lower()), df.columns[2])
    out = pd.DataFrame({"pathway": df[name_col].astype(str).str.strip(),
                        "adj.Pval": pd.to_numeric(df["adj.Pval"], errors="coerce")})
    out["key"] = out["pathway"].map(pathway_key)
    return out


# ── gene sets ─────────────────────────────────────────────────

def read_gmt(src):
    if isinstance(src, (str, Path)):
        raw = Path(src).read_bytes()
    elif isinstance(src, bytes):
        raw = src
    else:
        raw = src.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    text = raw.decode("utf-8", errors="replace")
    sets = {}
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            sets[parts[0]] = [g for g in parts[2:] if g]
    if not sets:
        raise InputError("No gene sets found; a .gmt has one set per line: name<TAB>description<TAB>genes...")
    return sets
