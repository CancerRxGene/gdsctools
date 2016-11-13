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
gd = GDSCElasticNet(ic50_v17, gf_v17)
drugid = 1047

#####################################################
# Find best model and corresponding alpha
res = gd.runCV(drugid, n_folds=10)
best_alpha = res.alpha

#####################################################
# Plot weights of best model
best_model = gd.get_model(alpha=best_alpha)
gd.plot_weight(drugid, model=best_model)

#####################################################
# increasing alpha 
model1 = gd.get_model(alpha=best_alpha*10.)
gd.plot_weight(drugid, model=model1, fontsize=9)

#####################################################
# decreasing alpha
model2 = gd.get_model(alpha=best_alpha/10.)
gd.plot_weight(drugid, model=model2, fontsize=9)

