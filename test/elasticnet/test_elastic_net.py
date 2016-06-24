from gdsctools.elastic_net import *
from gdsctools import IC50, GenomicFeatures, gdsctools_data
from nose.tools import assert_almost_equal
from math import log


def test_elastic_net():

    ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
    gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))
    en = ElasticNet(ic, gf)
    en.elastic_net(1047, show=True)


    # Tuning of the alpha parameters
    res = en.tune_alpha(1047, N=40, l1_ratio=0.1) 


    #assert_almost_equal(-log(res[3]), 3.5424386046062248, 7)


    res = en.plot_weight(1047, 0.01, l1_ratio=0.5)
    res = en.plot_importance(1047, 0.01, l1_ratio=0.5)

    from gdsctools import ic50_test
    gd = ElasticNet(ic50_test, gf)
    df = gd.elastic_all(0.1, stacked=True)


    #gd.plot_cindex(999, [0.1])

