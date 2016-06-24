"""
Analyse all associations (drug/feature)
=========================================

Volcano plot (all associations)
"""



#####################################################
#
from gdsctools import ANOVA, ic50_test
gdsc = ANOVA(ic50_test)
results = gdsc.anova_all()
results.volcano()


