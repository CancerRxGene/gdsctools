from gdsctools import gdsctools_data, ANOVA, GDSC
from gdsctools.gdsc import IC50Cluster
import os


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


def test_gdsc():

    import tempfile
    tempdir = tempfile.mkdtemp()
    compdir = tempdir + os.sep + "company_packages"
    tissuedir = tempdir + os.sep + "tissue_packages"

    pathtoGF = os.path.split(gdsctools_data("GF_BRCA_v17.csv.gz"))[0]
    ic50 = gdsctools_data('IC50_v17.csv.gz')
    DD = gdsctools_data("test_drug_decode2.csv")

    gg = GDSC(ic50, DD, pathtoGF+'/GF_*_v17.csv.gz')
    gg.company_directory = compdir
    gg.tissue_directory = tissuedir
    gg.analyse()
    assert gg.companies == ['COMPANY_A', 'COMPANY_B']
    gg.create_data_packages_for_companies()
    gg.create_summary_pages()







