import numpy as np
import pandas as pd
import pytest

from ezplot import io
from ezplot.project import Project

from .conftest import EXAMPLE


def test_example_kinds_and_direction():
    pr = Project()
    assert pr.load(EXAMPLE / "example_DE_matrix.csv", "example_DE_matrix.csv") == "de_matrix"
    assert pr.load(EXAMPLE / "example_DEG_Heatmap_Data.csv", "x.csv") == "heatmap"
    assert pr.group_names == ["WT", "KO"]
    assert (pr.control, pr.treatment) == ("WT", "KO")
    # the generator writes iDEP's sign: WT-KO_log2FC positive = higher in WT
    assert pr.up_group == "WT" and pr.direction_confidence > 0.95
    g = pr.oriented_genes()
    expr = pr.de.expr
    top = g.nsmallest(50, "padj").index
    ko_minus_wt = expr.loc[top, ["KO1", "KO2", "KO3", "KO4"]].mean(axis=1) - \
        expr.loc[top, ["WT1", "WT2", "WT3", "WT4"]].mean(axis=1)
    assert ((ko_minus_wt > 0) == (g.loc[top, "lfc"] > 0)).all()


def test_swapping_control_flips_sign():
    pr = Project()
    pr.load(EXAMPLE / "example_DE_matrix.csv", "de.csv")
    a = pr.oriented_genes()["lfc"]
    pr.control, pr.treatment = pr.treatment, pr.control
    assert (pr.oriented_genes()["lfc"] == -a).all()


@pytest.mark.parametrize("name,expected", [("KO5R", "KO"), ("KO_3", "KO"),
                                           ("WT-rep2", "WT"), ("KO10", "KO")])
def test_group_stems(name, expected):
    assert io.assign_groups([name])[name] == expected


def test_clean_symbol():
    assert io.clean_symbol("Gm31520ENSMUSG00000121703") == "Gm31520"
    assert io.clean_symbol("ENSMUSG00000121703") == "ENSMUSG00000121703"
    assert io.clean_symbol("Actb") == "Actb"


def _full_deseq2(tmp_path):
    """The example matrix with the extra columns of the patched iDEP's full DESeq2 export."""
    df = pd.read_csv(EXAMPLE / "example_DE_matrix.csv")
    se = np.random.default_rng(0).uniform(0.1, 2, len(df))
    at = df.columns.get_loc("WT-KO_log2FC") + 1
    df.insert(at, "WT-KO_lfcSE", se)
    df.insert(at + 1, "WT-KO_stat", df["WT-KO_log2FC"] / se)
    df.insert(at + 2, "WT-KO_pvalue", df["WT-KO_adjPval"] / 2)
    path = tmp_path / "Results_full_DESeq2.csv"
    df.to_csv(path, index=False)
    return path


def test_full_deseq2_ranks_by_oriented_stat(tmp_path):
    pr = Project()
    assert pr.load(_full_deseq2(tmp_path), "Results_full_DESeq2.csv") == "de_matrix_stat"
    assert pr.de.has_stat and pr.rank_metric == "stat"
    assert list(pr.de.expr.columns) == ["WT1", "WT2", "WT3", "WT4", "KO1", "KO2", "KO3", "KO4"]
    g = pr.oriented_genes()
    # stat shares log2FC's sign, so it is flipped with it (WT-KO positive = higher in WT)
    assert (np.sign(g["stat"]) == np.sign(g["lfc"])).all()
    assert (pr.ranking().to_numpy() == g["stat"].to_numpy()).all()
    stat_key = pr._cache_key("hallmarks")
    pr.gsea.rank_by = "lfc"
    assert pr.rank_metric == "lfc" and pr._cache_key("hallmarks") != stat_key
    assert (pr.ranking().to_numpy() == g["lfc"].to_numpy()).all()


def test_plain_matrix_falls_back_to_lfc():
    pr = Project()
    pr.load(EXAMPLE / "example_DE_matrix.csv", "de.csv")
    assert pr.gsea.rank_by == "stat" and not pr.de.has_stat
    assert pr.rank_metric == "lfc"
    assert (pr.ranking().to_numpy() == pr.oriented_genes()["lfc"].to_numpy()).all()
