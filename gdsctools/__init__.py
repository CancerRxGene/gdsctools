import pkg_resources
try:
    version = pkg_resources.require("gdsctools")[0].version
    __version__ = version
except:
    # update this manually is possible when the version in the
    # setup changes
    version = "0.1"

import easydev

ic50_test = easydev.get_share_file('gdsctools', 'data', 'IC50_10drugs.tsv')




from gdsctools.report import HTMLTable, Report
from gdsctools.readers import IC50, GenomicFeatures
from gdsctools.anova import ANOVA, ANOVAReport
from gdsctools.volcano import VolcanoANOVA
