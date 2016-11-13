"""
Tuning alpha (elastic net case)
=========================================

Elastic net requires tune an alpha parameter.
"""



#####################################################
#
from gdsctools import *
ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))
en = GDSCElasticNet(ic, gf)
en.tune_alpha(1047, alpha_range=(-3.5,-1), N=40, l1_ratio=0.1)


