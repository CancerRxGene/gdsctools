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
import os
import warnings
warnings.simplefilter('ignore', DeprecationWarning)

def set_logger_level(level="INFO"):
    import colorlog as logger
    logger.getLogger().setLevel(level)


try:
    version = pkg_resources.require("gdsctools")[0].version
    __version__ = version
except:
    version = "0.17"

try:
    license = open('../LICENSE', 'r').read()
except:
    license = '3-clause ("Simplified" or "New") BSD'


# To be defined before importing any modules
def gdsctools_data(filename, where=None):
    """Simple utilities to retrieve data sets from gdsctools/share directory"""
    import easydev
    gdsctools_path = easydev.get_package_location('gdsctools')
    share = os.sep.join([gdsctools_path, "gdsctools", 'data'])
    # in the code one may use / or \
    if where:
        filename = os.sep.join([share, where, filename])
    else:
        filename = os.sep.join([share, filename])
    if os.path.exists(filename) is False:
        raise Exception('unknown file %s' % filename)
    return filename

from gdsctools.readers import IC50, GenomicFeatures, DrugDecode, Reader
from gdsctools.anova import ANOVA
from gdsctools.anova_report import ANOVAReport
from gdsctools.anova_results import ANOVAResults
from gdsctools.settings import ANOVASettings
from gdsctools.datasets import *
from gdsctools.volcano import VolcanoANOVA
from gdsctools.cosmictools import COSMICInfo
from gdsctools.tissues import TCGA
from gdsctools.gdsc import GDSC, IC50Cluster
from gdsctools.stats import signed_effects
from gdsctools.regression import GDSCElasticNet, GDSCLasso, GDSCRidge
from gdsctools.omnibem import OmniBEMBuilder
from gdsctools.gdsc1000 import GDSC1000


def gdsctools_help(name=None):
    import easydev
    if name is None:
        easydev.onweb('http://gdsctools.readthedocs.org')
    else:
        url = "http://gdsctools.readthedocs.org/en/master/references.html"
        try:
            url += '#module-' + name.__module__
        except:
            print("Not a known gdsctools class or function")
        easydev.onweb(url)

