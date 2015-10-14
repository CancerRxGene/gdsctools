from gdsctools.anova import GDSC_ANOVA
from gdsctools.reader import IC50
import pandas as pd
from easydev import assert_list_almost_equal


def test_anova_one_drug_one_feature():

    r = IC50('test2.tsv')
    an = GDSC_ANOVA(r.ic50, r.features)

    # test 1 drug
    drug_id = 'Drug_999_IC50'
    df = an.anova_one_drug_one_feature(
        drug_id=drug_id,
        feature_name='ABCB1_mut')

    control = {'Drug id': {1: drug_id},
        'Drug name': {1: drug_id},
        'Drug Target': {1: drug_id},
        'FEATURE': {1: 'ABCB1_mut'},
        'FEATURE_ANOVA_pval': {1: 0.86842684367357359},
        'FEATURE_IC50_T_pval': {1: 0.48586107208790896},
        'FEATURE_IC50_effect_size': {1: 0.31407773405409201},
        'FEATURE_deltaMEAN_IC50': {1: -0.36105662590553411},
        'FEATUREneg_Glass_delta': {1: 0.31296252976074801},
        'FEATUREneg_IC50_sd': {1: 1.1536736560173895},
        'FEATUREneg_logIC50_MEAN': {1: 2.8007757068403043},
        'FEATUREpos_Glass_delta': {1: 0.53754504818376181},
        'FEATUREpos_IC50_sd': {1: 0.67167696386648634},
        'FEATUREpos_logIC50_MEAN': {1: 2.4397190809347702},
        'MSI_ANOVA_pval': {1: 0.14598946672374763},
        'N_FEATURE_neg': {1: 370},
        'N_FEATURE_pos': {1: 5},
        'Tissue_ANOVA_pval': {1: 3.2808255732569986e-06},
        'log max.Conc.tested': {1: None},
        'log max.Conc.tested2': {1: None}}
    control = pd.DataFrame(control)

    assert_list_almost_equal(df,control)


def test_anova_one_drug():
    # test entire drug across all fearures
    r = IC50('test2.tsv')
    an = GDSC_ANOVA(r.ic50, r.features)
    df = an.anova_one_drug('Drug_999_IC50')


def test_anova_all():
    r = IC50('test2.tsv')
    an = GDSC_ANOVA(r.ic50, r.features)
    df = an.anova_all()


# Need a test file with a mutation that can be ignored beacause number of pos is
# <2

# need test includeMSI_factor set to False
# need test analyseType != PANCAN
