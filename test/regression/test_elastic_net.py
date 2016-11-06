from gdsctools.regression import *
from gdsctools import IC50, GenomicFeatures, gdsctools_data
from nose.tools import assert_almost_equal
from math import log


def test_regression_elastic_net():
    ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
    gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))
    gd = GDSCElasticNet(ic, gf)

    drugid = 1047

    # automatic CV to get best model
    results = gd.runCV(drugid, n_folds=5)
    bestalpha = results.alpha
    best_model = gd.get_model(alpha=bestalpha)

    # Some plotting
    gd.plot_weight(drugid, best_model)
    gd.plot_importance(drugid, model=best_model)

    # manual fit
    scores = gd.fit(drugid, alpha=bestalpha)

    results = gd.tune_alpha(drugid, alpha_range=(-3.5,-1))

    res = gd.check_randomness(drugid, N=10)
    
    #res = gd.boxplot(drugid, model=best_model,
    #    n=10,bx_vert=False)
