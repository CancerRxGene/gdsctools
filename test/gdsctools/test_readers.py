from gdsctools.readers import GenomicFeatures, IC50, PANCAN
from easydev import TempFile

def get_data(filename='test2.tsv'):
    import os
    try:
        r = IC50(filename)
    except:
        r = IC50('test/gdsctools' + os.sep + filename)
    return r


def test_read_ic50():
    r = IC50(get_data())
    # we can also instanciate from a valid dataframe
    r = IC50(r)

    print(r)
    r.hist()
    r.plot_ic50_count()
    r.cosmicIds

    f = TempFile()
    r.to_csv(f.name)
    f.delete()

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


def test_pancan_reader_rdata():
    r = PANCAN()
    len(r.df)
