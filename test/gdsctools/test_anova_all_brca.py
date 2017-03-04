from gdsctools import ANOVA, gdsctools_data, IC50
from gdsctools import ANOVAReport
from nose.tools import assert_almost_equal
import numpy as np


def test_anova_brca():

    an1 = ANOVA(gdsctools_data('IC50_v17.csv.gz'))
    an1.set_cancer_type('breast')

    an = ANOVA(an1.ic50, gdsctools_data('GF_BRCA_v17.csv.gz'))
    dfori = an.anova_all()

    df = dfori.df.sum()
    df = df.drop(['DRUG_TARGET', 'DRUG_NAME', 'DRUG_ID', 'FEATURE'])
    df = df.fillna(0)
    totest = df.to_dict()

    exact = {'ANOVA_FEATURE_FDR': 1133416.7761055394,
        'ANOVA_FEATURE_pval': 5824.8201538614458,
        'FEATURE_IC50_T_pval': 5824.8201538614449,
        'FEATURE_IC50_effect_size': 4408.511449781573,
        'FEATURE_delta_MEAN_IC50': 261.11373729866705,
        'FEATURE_neg_Glass_delta': 4487.7401723134735,
        'FEATURE_neg_IC50_sd': 14701.868130868914,
        'FEATURE_neg_logIC50_MEAN': 28701.510736736222,
        'FEATURE_pos_Glass_delta': 6536.8938399490198,
        'FEATURE_pos_IC50_sd': 13362.588398939894,
        'FEATURE_pos_logIC50_MEAN': 28962.624474034845,
        'ANOVA_MSI_pval': 0.0,
        'N_FEATURE_neg': 439196.0,
        'N_FEATURE_pos': 92140.0,
        'ANOVA_TISSUE_pval': 0.0,
        'ASSOC_ID': 68509365.0,
        'index': 68497660.0}

    for k, v in totest.items():
        if k in ['ANOVA_MEDIA_pval']:
            continue
        assert_almost_equal(v, exact[k])

    # test part of the report (summary section)
    r = ANOVAReport(an, dfori)
    totest = r.diagnostics().to_dict()

    exact = {'text': {0: 'Type of analysis',
  1: 'Total number of possible drug/feature associations',
  2: 'Total number of ANOVA tests performed',
  3: 'Percentage of tests performed',
  4: '',
  5: 'Total number of tested drugs',
  6: 'Total number of genomic features used',
  7: 'Total number of screened cell lines',
  8: 'MicroSatellite instability included as factor',
  9: '',
  10: 'Total number of significant associations',
  11: ' - sensitive',
  12: ' - resistant',
  13: 'p-value significance threshold',
  14: 'FDR significance threshold',
  15: 'Range of significant p-values',
  16: 'Range of significant % FDRs'},
 'value': {0: 'breast',
  1: 13780,
  2: 11705,
  3: 84.94,
  4: '',
  5: 265,
  6: 52,
  7: 51,
  8: False,
  9: '',
  10: 27,
  11: 17,
  12: 10,
  13: 0.001,
  14: 25,
  15: '[2.098e-09, 0.0004356]',
  16: '[0.002456 18.89]'}}

 
    assert totest == exact

    #import shutil
    #shutil.rmtree('breast')

