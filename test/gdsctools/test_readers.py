from gdsctools.readers import GenomicFeatures, IC50, PANCAN, DrugDecode
from gdsctools.readers import Reader
from easydev import TempFile
from gdsctools import ic50_test
import pandas as pd
 
from gdsctools.datasets import testing


def test_readers():
    a = Reader()

    try:
        a = Reader('stupido')
        assert False
    except:
        assert True
    try:
        a = Reader(1)
        assert False
    except:
        assert True


def test_read_ic50():
    # -------------------------------- functionalities
    r = IC50(ic50_test)
    # we can also instanciate from a valid dataframe
    r = IC50(r)

    print(r)
    r.hist()
    r.plot_ic50_count()
    r.cosmicIds

    f = TempFile()
    r.to_csv(f.name)
    f.delete()

    # columns may be duplicated
    r = IC50(ic50_test)
    df = pd.concat([r.df, r.df['Drug_999_IC50']], axis=1)
    # create new instance that should raise an error
    try:
        IC50(df)
        assert False
    except:
        assert True

    # ---------------------------------------- different IC50 formats
    # test all files available
    for key in testing.keys() :
        filename = testing[key].location
        if filename.startswith('ic50_test'):
            ic = IC50(filename)
    # some specific checks:
    #ic = IC50(testing['ic50_test_header_drug_prefix_only'].location)
    #assert ic.df.shape == (2,2)
    #assert all(ic.df.columns == ['1','2'])
    ic = IC50(testing['ic50_test_header_no_drug_prefix'].location)
    assert ic.drugIds == ['1', '2']

    ic = IC50(testing['ic50_test_header_drug_prefix_only'].location)
    assert ic.drugIds == ['Drug_1_IC50', 'Drug_2_IC50']

    ic = IC50(testing['ic50_test_header_mixed_drug_prefix'].location)
    assert ic.drugIds == ['Drug_1_IC50', 'Drug_2_IC50']


def test_read_gf():
    # Reads a default file
    r = GenomicFeatures()

    # we can also instanciate from another GenomicFeatures instance
    r = GenomicFeatures(r)
 
    # we can also instanciate from a valid dataframe
    r = GenomicFeatures(r.df)

    print(r)
    r.features
    r.tissues
    r.plot()
    r.drop_tissue_in('breast')
    r.drop_tissue_in(['skin', 'bone'])
    r.keep_tissue_in(['cervix', 'lung'])
    assert r.shift == 2

    assert len(r.unique_tissues) == 2

    gf1 = GenomicFeatures()

    gf2 = GenomicFeatures(testing.genomic_features_csv)
    to_drop = [x for x in gf1.df.index if x not in gf2.df.index]

    gf1.drop_cosmic(to_drop)
    gf1.features = gf2.features

    assert gf2 == gf1

    gf = GenomicFeatures(testing.genomic_features_bare_csv)
    assert gf.shift == 1

def _test_pancan_reader_rdata():
    r = PANCAN()
    len(r.df)


def test_drugs():
    r1 = DrugDecode(testing.drug_test_csv)
    r1.drugIds
    r2 = DrugDecode(testing.drug_test_tsv)
    r2.drugIds
    assert r1 == r2
