from gdsctools.elastic_net import *
from gdsctools import IC50, GenomicFeatures, gdsctools_data



def test_elastic_net():

    ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
    gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))
    en = ElasticNet(ic, gf)
    en.elastic_net(1047, show=True)
