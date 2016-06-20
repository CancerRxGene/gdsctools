"""
Plot weights resulting from an Elastic Net analysis
====================================================

Note that only the 50 most important weigths are shown 
"""
#####################################################
# We look at the" effect of the alpha parameter on 
# the weights returned by the elastic net analysis
from gdsctools import *

#####################################################
# First we alpha=0.01
gd = elastic_net.ElasticNet(ic50_v17, gf_v17)
gd.plot_weight(1047, alpha=0.1, fontsize=9)

#####################################################
# increasing alpha 
gd.plot_weight(1047, alpha=0.5, fontsize=9)

#####################################################
# decreasing alpha
gd.plot_weight(1047, alpha=0.01, fontsize=9)

