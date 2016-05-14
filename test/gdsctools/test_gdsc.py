from gdsctools import GDSC, gdsctools_data
from gdsctools.gdsc import IC50Cluster


def test_IC50Cluster():
    ic50 = IC50Cluster(gdsctools_data("test_v18_clustering.tsv"))

    # there are still 6 duplicared drugs not clustered because they 
    # have lots of valid data in each
    assert len(ic50.info()) == 6
    assert 1510 in ic50.info().DRUG_ID.values

    assert "Drug_1032_0.5_IC50" in ic50.df.columns
    assert "Drug_1032_2_IC50" in ic50.df.columns

    ic50.cleanup()

    # In this data set, a drug is reported 3 times (1211) and should appear 
    # as follows:
    assert 1211 in ic50.df.columns
    assert 11211 in ic50.df.columns
    assert 21211 in ic50.df.columns
