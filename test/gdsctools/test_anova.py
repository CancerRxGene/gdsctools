from gdsctools.anova import ANOVA
import pandas as pd
from easydev import assert_list_almost_equal
from nose.tools import assert_almost_equal
from gdsctools import ic50_test

def test_anova_one_drug_one_feature():

    an = ANOVA(ic50_test)

    # test 1 drug
    drug_id = 'Drug_999_IC50'
    df = an.anova_one_drug_one_feature(
        drug_id=drug_id,
        feature_name='ABCB1_mut', show=True)

    control = {'DRUG_ID': {1: drug_id},
        'DRUG_NAME': {1: drug_id},
        'DRUG_TARGET': {1: drug_id},
        'FEATURE': {1: 'ABCB1_mut'},
        'ANOVA_FEATURE_pval': {1: 0.86842684367357359},
        'FEATURE_IC50_T_pval': {1: 0.48586107208790896},
        'FEATURE_IC50_effect_size': {1: 0.31407773405409201},
        'FEATURE_delta_MEAN_IC50': {1: -0.36105662590553411},
        'FEATUREneg_Glass_delta': {1: 0.31296252976074801},
        'FEATUREneg_IC50_sd': {1: 1.1536736560173895},
        'FEATUREneg_logIC50_MEAN': {1: 2.8007757068403043},
        'FEATUREpos_Glass_delta': {1: 0.53754504818376181},
        'FEATUREpos_IC50_sd': {1: 0.67167696386648634},
        'FEATUREpos_logIC50_MEAN': {1: 2.4397190809347702},
        'ANOVA_MSI_pval': {1: 0.14598946672374763},
        'N_FEATURE_neg': {1: 370},
        'N_FEATURE_pos': {1: 5},
        'ANOVA_TISSUE_pval': {1: 3.2808255732569986e-06}}
    control = pd.DataFrame(control)

    assert_list_almost_equal(df,control)


def test_anova_one_drug():
    # test entire drug across all fearures
    an = ANOVA(ic50_test)
    df = an.anova_one_drug('Drug_999_IC50')


def test_anova_all():
    an = ANOVA(ic50_test)
    # slow, let us cut the features to keep only tenish
    # this is a trick but would be nice to have this in the API
    features = an.features.df
    features = features[features.columns[0:12]]

    an = ANOVA(ic50_test, features)
    results = an.anova_all()
    assert_almost_equal( results.df['ANOVA_FEATURE_FDR'].sum(), 
            10312.23061065521, 6)

    an2 = ANOVA(ic50_test, features, set_media_factor=True)
    results = an2.anova_all()
    assert_almost_equal( results.df['ANOVA_FEATURE_FDR'].sum(), 
            10238.529313503008, 6)


def test_odof_with_without_media():

    gdsc = ANOVA(ic50_test)
    _res = gdsc.anova_one_drug_one_feature('Drug_1047_IC50', 'TP53_mut')
    dd1 = gdsc._get_anova_summary(gdsc.data_lm, output='dict')
    assert_list_almost_equal([dd1['feature'], dd1['msi'], dd1['tissue']], 
        [1.5750735472022118e-58,  0.025902887791637515, 
            5.541879283763767e-44])

    gdsc = ANOVA(ic50_test, set_media_factor=True)
    _res = gdsc.anova_one_drug_one_feature('Drug_1047_IC50', 'TP53_mut')
    dd2 = gdsc._get_anova_summary(gdsc.data_lm, output='dict')
    assert_list_almost_equal([dd2['feature'], dd2['media'], dd2['msi'],
    dd2['tissue']], [ 2.9236500715529455e-58, 0.7762487502315283,
         0.023777744527686766, 1.5729157319290974e-44])





def test_odof():
    pass
    # gdsc.anova_one_drug_one_feature('Drug_1013_IC50', 'BCR-ABL_mut').TOut[65]: 
    #                                   1
    #ANOVA_FEATURE_pval              0.176635
    #ANOVA_MEDIA_pval                0.720394
    #ANOVA_MSI_pval                  0.720394




def test_anova_summary():
    an = ANOVA(ic50_test)
    # by default  regression includes + msi + feature
    drug_id = 'Drug_999_IC50'

    df = an.anova_one_drug_one_feature(drug_id, 'ASH1L_mut')

    x = an.anova_pvalues
    y = [3.210453608523738e-06, 0.14579091345305398, 0.5430736275249095, None]
    assert_list_almost_equal(x, y)

    
    an.settings.analysis_type = 'COREAD' # something different from PANCAN
    df = an.anova_one_drug_one_feature(drug_id, 'ASH1L_mut')
    x = an.anova_pvalues
    y = [0.262294448831941, 0.30599483315087317, None]
    assert_list_almost_equal(x, y)

    # now remove also the MSI factor
    an.settings.include_MSI_factor = False
    df = an.anova_one_drug_one_feature(drug_id, 'ASH1L_mut')
    x = an.anova_pvalues
    y = [0.21266050833611852, None]
    assert_list_almost_equal(x, y)

    assert (df.N_FEATURE_neg == 365).all()







