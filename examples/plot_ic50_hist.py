"""
Histogram of the IC50 from v17 data set
=========================================

histogram of the IC50
"""



#####################################################
#
from gdsctools import IC50, ic50_v17
r = IC50(ic50_v17)
r.hist()
