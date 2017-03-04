from gdsctools.anova import ANOVA
from gdsctools import gdsctools_data
import pandas as pd
from easydev import assert_list_almost_equal
from gdsctools import ic50_test

import warnings
warnings.filterwarnings(action="ignore")


def assert_almost_equal(x, delta=1e-6):
    import pytest
    pytest.approx(x,abs=delta)

def test_anova_one_drug_one_feature():

    an = ANOVA(ic50_test)
    # test 1 drug
    drug_id = 999
    df = an.anova_one_drug_one_feature(
        drug_id=drug_id,
        feature_name='ABCB1_mut', show=True)
    df['DRUG_ID'] = drug_id
    df['DRUG_NAME'] = drug_id
    df['DRUG_TARGET'] = drug_id
    df['ANOVA_MEDIA_pval'] = -1

    control = {'DRUG_ID': {1: drug_id},
        'DRUG_NAME': {1: drug_id},
        'DRUG_TARGET': {1: drug_id},
        'FEATURE': {1: 'ABCB1_mut'},
        'ANOVA_FEATURE_pval': {1: 0.86842684367357359},
        'FEATURE_IC50_T_pval': {1: 0.48586107208790896},
        'FEATURE_IC50_effect_size': {1: 0.31407773405409201},
        'FEATURE_delta_MEAN_IC50': {1: -0.36105662590553411},
#        'FEATUREneg_Glass_delta': {1: 0.31296252976074801},
#        'FEATUREneg_IC50_sd': {1: 1.1536736560173895},
        'FEATUREneg_logIC50_MEAN': {1: 2.8007757068403043},
#        'FEATUREpos_Glass_delta': {1: 0.53754504818376181},
#        'FEATUREpos_IC50_sd': {1: 0.67167696386648634},
        'FEATUREpos_logIC50_MEAN': {1: 2.4397190809347702},
        'ANOVA_MSI_pval': {1: 0.14598946672374763},
        'N_FEATURE_neg': {1: 370},
        'N_FEATURE_pos': {1: 5},
        'ANOVA_MEDIA_pval': {1: -1},
        'ANOVA_TISSUE_pval': {1: 3.2808255732569986e-06}}
    control = pd.DataFrame(control)
    for this in control.columns:
        if this in ["FEATURE_delta_MEAN_IC50", "FEATUREneg_logIC50_MEAN", 
            "FEATUREpos_logIC50_MEAN"]:
            pass
        else:
            assert_almost_equal(df[this].values[0],control[this].values[0])


def test_compare_formula_vs_gdsc():
    # TISSUE + MSI + feature 
    an = ANOVA(ic50_test)
    drug_id = 999
    df1 = an.anova_one_drug_one_feature(
        drug_id=drug_id, feature_name='ABCB1_mut')
    an.settings.regression_formula = "Y ~ C(tissue) + C(msi) + feature"
    df2 = an.anova_one_drug_one_feature(
        drug_id=drug_id, feature_name='ABCB1_mut')

    assert all(df1.fillna(1).values[0] == df2.fillna(1).values[0])

    # MSI + feature only
    an = ANOVA(ic50_test)
    an.settings.analysis_type = "breast"
    df1 = an.anova_one_drug_one_feature(
        drug_id=drug_id, feature_name='ABCB1_mut')
    an.settings.regression_formula = "Y ~ C(msi) + feature"
    df2 = an.anova_one_drug_one_feature(an.drugIds[2], an.feature_names[0])
    assert all(df1.fillna(1) == df2.fillna(1))


def test_anova_one_drug():
    # test entire drug across all fearures
    an = ANOVA(ic50_test)
    df = an.anova_one_drug(999)

# This also test ANOVAResults
def test_anova_all():
    an = ANOVA(ic50_test)
    # slow, let us cut the features to keep only ten-ish values 
    # this is a trick but would be nice to have this in the API
    features = an.features.df
    features = features[features.columns[0:12]]

    an = ANOVA(ic50_test, features)
    results = an.anova_all()
    results # test __repr__
    results.volcano()
    results.barplot_effect_size()

    assert_almost_equal(results.df['ANOVA_FEATURE_FDR'].sum()- 
            10312.23061065521)

    an2 = ANOVA(ic50_test, features, set_media_factor=True)
    results = an2.anova_all()
    assert_almost_equal( results.df['ANOVA_FEATURE_FDR'].sum() - 10238.529313503008)


def test_odof_with_without_media():

    gdsc = ANOVA(ic50_test)
    _res = gdsc.anova_one_drug_one_feature(1047, 'TP53_mut')
    odof = gdsc._get_one_drug_one_feature_data(1047, 'TP53_mut')
    dd1 = gdsc._get_anova_summary(gdsc.data_lm, output='dict', odof=odof)

    assert_almost_equal(dd1['feature']- 1.5750735472022118e-58, delta=1e-64)
    assert_almost_equal(dd1['msi']-  0.025902887791637515, delta=1e-10)

    # Change in Nov 2016
    #assert_almost_equal(dd1['tissue'], 5.541879283763767e-44, delta=1e-50)
    assert_almost_equal(dd1['tissue']- 1.0258702741509e-44, delta=1e-50)

    gdsc = ANOVA(ic50_test, set_media_factor=True)
    _res = gdsc.anova_one_drug_one_feature(1047, 'TP53_mut')
    odof = gdsc._get_one_drug_one_feature_data(1047, 'TP53_mut')
    dd2 = gdsc._get_anova_summary(gdsc.data_lm, output='dict', odof=odof)

    assert_almost_equal(dd2['feature']- 2.9236500715529455e-58, delta=1e-64)
    assert_almost_equal(dd2['media']- 0.7762487502315283, delta=1e-10)
    assert_almost_equal(dd2['msi']- 0.023777744527686766, delta=1e-10)
    assert_almost_equal(dd2['tissue']- 1.5729157319290974e-44, delta=1e-50)


def test_anova_summary():
    an = ANOVA(ic50_test)
    # by default  regression includes + msi + feature
    drug_id = 999
    df = an.anova_one_drug_one_feature(drug_id, 'ASH1L_mut')

    x = an.anova_pvalues
    x = [x["tissue"], x["msi"], x["feature"]]
    y = [3.210453608523738e-06, 0.14579091345305398, 0.5430736275249095]
    assert_list_almost_equal(x, y, deltas=1e-10)

    an.settings.analysis_type = 'COREAD' # something different from PANCAN
    df = an.anova_one_drug_one_feature(drug_id, 'ASH1L_mut')
    x = an.anova_pvalues
    x = [x["msi"], x["feature"]]
    y = [0.262294448831941, 0.30599483315087317]
    assert_list_almost_equal(x, y, deltas=1e-10)

    # now remove also the MSI factor, in which case the tissue must also be
    # removed !
    an.settings.include_MSI_factor = False
    an.settings.analysis_type = "COREAD"
    df = an.anova_one_drug_one_feature(drug_id, 'ASH1L_mut')
    x = [an.anova_pvalues["feature"]]
    y = [0.21266050833611852]
    assert_list_almost_equal(x, y, deltas=1e-10)

    assert (df.N_FEATURE_neg == 365).all()


def test_set_cancer_type():
    an = ANOVA(gdsctools_data("IC50_v17.csv.gz"))
    an.set_cancer_type("breast")
    assert_list_almost_equal([an.ic50.df.sum().sum()], [27721.255627472943])


def test_multicore():

    an = ANOVA(ic50_test)
    results = an.anova_all()

    an2 = ANOVA(ic50_test)
    results2 = an2.anova_all(multicore=2)

    all(results.df.fillna(0) == results2.df.fillna(0))


