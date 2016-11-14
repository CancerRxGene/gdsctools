from gdsctools import gdsctools_data, ANOVA
from gdsctools.gdsc import IC50Cluster


def test_IC50Cluster():
    dataset = gdsctools_data("test_v18_clustering.tsv")
    ic50 = IC50Cluster(dataset)


    # In this data set, a drug is reported 3 times (1211) and should appear 
    # as follows:
    assert 1211 in ic50.df.columns
    assert 11211 in ic50.df.columns
    assert 21211 in ic50.df.columns

    assert len(ic50.drugIds) == 860
    assert len(ic50.df) == 50


    an = ANOVA(ic50, dataset)
    an.diagnostics()['feasible_tests'] == 65026










