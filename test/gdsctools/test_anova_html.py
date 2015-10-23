from gdsctools.anova import GDSC_ANOVA, GDSC_ANOVA_Results
from gdsctools.reader import IC50
import pandas as pd
from easydev import assert_list_almost_equal
from . import tools



def test_html():
    # same as above, could factorise
    an = GDSC_ANOVA(tools.get_data())
    features = an.features.df
    features = features[features.columns[0:13]]
    an = GDSC_ANOVA(tools.get_data(), features)
    df = an.anova_all()

    r = GDSC_ANOVA_Results(df, 
            ic50=tools.get_data())
    r.settings.fdr_threshold= 90
    # long but should cover everthinh.
    r.create_html_pages()

