"""
Association between a Drug and a Feature
=========================================

Boxplot association
"""



#####################################################
#
from gdsctools import ANOVA, ic50_test
gdsc = ANOVA(ic50_test)
gdsc.set_cancer_type('breast')
df = gdsc.anova_one_drug_one_feature(1047, 'TP53_mut', show=True)
