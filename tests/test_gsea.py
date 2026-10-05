import numpy as np
import pandas as pd
import pytest

from ezplot import genesets
from ezplot.gsea.engine import _bh, _es, prerank
from ezplot.project import Project

from .conftest import EXAMPLE


def test_es_matches_bruteforce():
    rng = np.random.default_rng(0)
    n, k = 300, 20
    w = np.sort(np.abs(rng.normal(size=n)))[::-1]
    pos = np.sort(rng.choice(n, k, replace=False))
    hit = np.zeros(n, bool)
    hit[pos] = True
    run = np.cumsum(np.where(hit, w / w[hit].sum(), -1 / (n - k)))
    brute = run[np.argmax(np.abs(run))]
    assert _es(pos[None, :], w, n)[0] == pytest.approx(brute)


def test_bh():
    p = np.array([0.01, 0.04, 0.03, 0.2])
    assert _bh(p) == pytest.approx([0.04, 0.16 / 3, 0.16 / 3, 0.2])


def test_example_planted_signal():
    pr = Project()
    pr.load(EXAMPLE / "example_DE_matrix.csv", "de.csv")
    res = pr.run_gsea("hallmarks").set_index("pathway")
    assert res.loc["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "NES"] < -2
    assert res.loc["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "NES"] < -2
    assert res.loc["HALLMARK_INTERFERON_ALPHA_RESPONSE", "NES"] > 2
    assert (res.loc[["HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                     "HALLMARK_INTERFERON_ALPHA_RESPONSE"], "padj"] < 1e-4).all()


def test_seed_reproducible():
    pr = Project()
    pr.load(EXAMPLE / "example_DE_matrix.csv", "de.csv")
    sets = genesets.load("hallmarks", "mouse")
    a = prerank(pr.ranking(), sets, seed=7)
    b = prerank(pr.ranking(), sets, seed=7)
    pd.testing.assert_frame_equal(a, b)
