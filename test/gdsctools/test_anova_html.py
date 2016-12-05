from gdsctools.anova import ANOVA
from gdsctools.anova_report import ANOVAReport
from gdsctools import ic50_test
import numpy as np


def test_html():
    # same as above, could factorise
    an = ANOVA(ic50_test)
    features = an.features.df
    features = features[features.columns[0:30]]
    an = ANOVA(ic50_test, features)
    #an.settings.include_media_factor = False

    an.settings.pvalue_threshold = np.inf
    results = an.anova_all()

    assert len(results.df) == 302 

    r = ANOVAReport(gdsc=an, results=results)
    assert len(r.get_significant_set()) == 3
    # long but should cover everthinh.
    try:
        r.create_html_pages()
    except Exception as err:
        raise err
    finally:
        import shutil
        #shutil.rmtree('html_gdsc_anova')


