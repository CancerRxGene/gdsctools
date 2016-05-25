from gdsctools import ANOVA, ic50_test
from easydev import assert_list_almost_equal
from gdsctools.boxplots import BoxPlots


def test_get_boxplot_data():
    an = ANOVA(ic50_test)
    odof = an._get_one_drug_one_feature_data(1047,'TP53_mut')

    bb = BoxPlots(odof)

    data = bb._get_boxplot_data(mode='msi')
    assert data[1] == ['***MSI-stable neg', '***MSI-stable pos',
                  '**MSI-unstable neg',  '**MSI-unstable pos']
    expected = [2.0108071495663922e-47, 0.0012564798887037905]
    assert_list_almost_equal([data[2][0], data[2][1]], expected)  



