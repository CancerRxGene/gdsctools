from gdsctools.tissues import TCGA, TCGA_GDSC1000, TCGA_2_GDSC, Tissues


def test_TCGA():
    t = TCGA()
    assert t['ACC'] == 'Adrenocortical Carcinoma'

def test_tcga_gdsc1000():
    assert 'BLCA' in TCGA_GDSC1000

def test_TCGA_2_GDSC():
    assert TCGA_2_GDSC['MB'] == 'nervous_system'


def test_Tissues():
    t = Tissues()
    assert 'leukemia' in t.v17
    assert 'leukemia' in t.v18
