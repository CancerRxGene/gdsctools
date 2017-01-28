from gdsctools.regression import *
from gdsctools import IC50, GenomicFeatures, gdsctools_data
from nose.tools import assert_almost_equal
from math import log


def test_regression_elastic_net():
    ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
    gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))
    gd = GDSCElasticNet(ic, gf, verbose=True)
    print(gd)
    drugid = 1047

    # automatic CV to get best model
    results = gd.runCV(drugid, kfolds=5)
    bestalpha = results.alpha
    results.coefficients
    best_model = gd.get_model(alpha=bestalpha)

    # Some plotting
    gd.plot_weight(drugid, best_model)
    #gd.plot_weight(drugid, best_model, Nmax=4)

    gd.plot_importance(drugid, model=best_model)
    gd.plot_importance(drugid, model=None)

    # manual fit
    scores = gd.fit(drugid, alpha=bestalpha)

    results = gd.tune_alpha(drugid, alpha_range=(-3.5,-1))

    res = gd.check_randomness(drugid, N=10)

    res = gd.boxplot(drugid, model=best_model, n=10, bx_vert=False)
    res = gd.boxplot(drugid, model=best_model, n=10, bx_vert=True)

    gd.dendogram_coefficients()
    gd.dendogram_coefficients(stacked=False)


def test_ridge():
    ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
    gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))
    gd = GDSCRidge(ic, gf, verbose=True)
    gd.runCV(1047)


def test_lasso():
    ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
    gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))
    gd = GDSCLasso(ic, gf, verbose=True)
    gd.runCV(1047).alpha
