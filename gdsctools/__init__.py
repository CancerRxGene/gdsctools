"""GDSCTools

::

    from gdsctools import ANOVA, ic50_test, ANOVAReport
    gdsc = ANOVA(ic50_test)
    results = gdsc.anova_all()

    report = ANOVAReport(gdsc, results)
    report.create_html_pages()


Please See documentation on gdsctools.readthedocs.org

"""
import pkg_resources
try:
    version = pkg_resources.require("gdsctools")[0].version
    __version__ = version
except:
    # update this manually is possible when the version in the
    # setup changes
    version = "0.3"

try:
    license = open('../LICENSE', 'r').read()
except:
    license = '3-clause ("Simplified" or "New") BSD'


from gdsctools.report import HTMLTable, Report
from gdsctools.readers import IC50, GenomicFeatures, DrugDecoder
from gdsctools.anova import ANOVA, ANOVAReport 
from gdsctools.settings import ANOVASettings
from gdsctools.volcano import VolcanoANOVA
from gdsctools.datasets import *

from easydev import browser



def gdsctools_help(name=None):
    if name is None:
        browser.browse('http://gdsctools.readthedocs.org')
    else:
        url = "http://gdsctools.readthedocs.org/en/master/references.html"
        #url += "#" + name.__module__ + "."+ name
        try:
            url += '#module-' + name.__module__
        except:
            print("Not a known gdsctools class or function")
        browser.browse(url)

