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


def test_pathway_key_matches_msigdb_and_idep():
    assert io.pathway_key("HALLMARK_TNFA_SIGNALING_VIA_NFKB") == \
        io.pathway_key("HALLMARK TNFA SIGNALING VIA NFKB")
    assert io.pathway_key("GOBP_TYPE_B_PANCREATIC_CELL_PROLIFERATION") == \
        io.pathway_key("Type b pancreatic cell proliferation ")
