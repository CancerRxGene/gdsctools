from gdsctools.readers import GenomicFeatures, IC50, PANCAN, DrugDecoder, COSMIC
from easydev import TempFile
from gdsctools import ic50_test, drug_test
import pandas as pd


def test_cosmic():
    r = COSMIC(924100)
    r.on_web()


def test_read_ic50():
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
    assert len(r.features) == 382
    assert len(r.unique_tissues) == 2

def _test_pancan_reader_rdata():
    r = PANCAN()
    len(r.df)


def test_drugs():
    r = DrugDecoder(drug_test)
    r.drugIds
