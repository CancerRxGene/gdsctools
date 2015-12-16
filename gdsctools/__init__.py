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
import easydev

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


from gdsctools.readers import IC50, GenomicFeatures, DrugDecode
from gdsctools.anova import ANOVA
from gdsctools.anova_report import ANOVAReport 
from gdsctools.anova_results import ANOVAResults
from gdsctools.settings import ANOVASettings
from gdsctools.volcano import VolcanoANOVA
from gdsctools.datasets import *
from gdsctools.cosmictools import COSMICInfo
from gdsctools.tissues import TCGA
from gdsctools.gdsc import GDSC

def gdsctools_data(filename, where=None):
    """Simple utilities to retrieve data sets from gdsctools/share directory"""
    import os
    share = easydev.get_shared_directory_path('gdsctools') 
    share = os.sep.join([share, 'data'])
    # in the code one may use / or \ 
    if where:
        filename = os.sep.join([share, where, filename])
    else:
        filename = os.sep.join([share, filename])
    if os.path.exists(filename) is False:
        raise Exception('unknown file %s' % filename)
    return filename


def gdsctools_help(name=None):
    if name is None:
        easydev.onweb('http://gdsctools.readthedocs.org')
    else:
        url = "http://gdsctools.readthedocs.org/en/master/references.html"
        #url += "#" + name.__module__ + "."+ name
        try:
            url += '#module-' + name.__module__
        except:
            print("Not a known gdsctools class or function")
        easydev.onweb(url)

