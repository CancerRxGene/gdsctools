from gdsctools.reader import GenomicFeatures, IC50

def get_data(filename='test2.tsv'):
    import os
    try:
        r = IC50(filename)
    except:
        r = IC50('test/gdsctools' + os.sep + filename)
    return r


def test_read_ic50():
    r = IC50(get_data())
    print(r)

def test_read_gf():
    r = GenomicFeatures()
    print(r)
