"""Generate the synthetic example dataset shipped in ezplot/example/ (iDEP file formats).

Mouse genes from the bundled Hallmark sets; 4 WT vs 4 KO samples of simulated log2
expression with planted shifts (oxidative phosphorylation and EMT lower in KO,
interferon-alpha response higher in KO, plus a few strong single genes). Statistics are a
t-test + BH, standing in for DESeq2. Column naming and sign follow iDEP:
"WT-KO_log2FC" is positive for genes higher in WT.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from ezplot import genesets

rng = np.random.default_rng(2024)
out = Path(__file__).resolve().parents[1] / "ezplot" / "example"
out.mkdir(exist_ok=True)

h = genesets.load("hallmarks", "mouse")
bp = genesets.load("go_bp", "mouse")
genes = sorted({g for s in list(h.values()) + list(bp.values())[:1500] for g in s})
genes = list(rng.choice(genes, size=min(9000, len(genes)), replace=False))
n = len(genes)
base = rng.normal(7, 2, n).clip(1, 16)
effect = np.zeros(n)
idx = {g: i for i, g in enumerate(genes)}


def shift(gs, mean, frac=0.6):
    members = [idx[g] for g in gs if g in idx]
    hit = rng.choice(members, size=int(len(members) * frac), replace=False)
    effect[hit] += rng.normal(mean, abs(mean) * 0.5, len(hit))


shift(h["HALLMARK_OXIDATIVE_PHOSPHORYLATION"], -0.45)
shift(h["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"], -0.5)
shift(h["HALLMARK_INTERFERON_ALPHA_RESPONSE"], 0.6)
shift(h["HALLMARK_MYOGENESIS"], -0.3, 0.4)
strong = rng.choice(n, 40, replace=False)
effect[strong] = rng.choice([-1, 1], 40) * rng.uniform(1.2, 5, 40)

wt = base[:, None] + rng.normal(0, 0.3, (n, 4))
ko = base[:, None] + effect[:, None] + rng.normal(0, 0.3, (n, 4))
t = stats.ttest_ind(wt, ko, axis=1, equal_var=True)
p = t.pvalue
o = np.argsort(p)
q = np.empty(n)
q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1].clip(max=1)
lfc = wt.mean(1) - ko.mean(1)          # iDEP sign: WT-KO positive = higher in WT

samples = ["WT1", "WT2", "WT3", "WT4", "KO1", "KO2", "KO3", "KO4"]
expr = np.hstack([wt, ko])
ens = [f"ENSMUSG{100000 + i:011d}" for i in range(n)]
de = pd.DataFrame({"symbol": genes, "ensembl_ID": ens, "User_ID": ens,
                   "baseMean": np.round(2 ** expr.mean(1), 3),
                   "WT-KO_log2FC": lfc, "WT-KO_adjPval": q, "Processed data:": "",
                   "search_label": [f"{g} | {e}" for g, e in zip(genes, ens)]})
de = pd.concat([de, pd.DataFrame(expr, columns=samples)], axis=1).round(6)
de = de.sort_values("baseMean", ascending=False)
de.to_csv(out / "example_DE_matrix.csv", index=False)

sig = (q < 0.05) & (np.abs(lfc) >= 1)
top = np.flatnonzero(sig)
hm = pd.DataFrame(expr[top], columns=samples, index=[genes[i] for i in top])
hm = hm.sub(hm.mean(axis=1), axis=0)   # iDEP exports row-centred values
hm.insert(0, "Gene_ID", [ens[i] for i in top])
hm.to_csv(out / "example_DEG_Heatmap_Data.csv")
print(f"{n} genes, {sig.sum()} DEGs (FDR<0.05, |log2FC|>=1) -> {out}")
