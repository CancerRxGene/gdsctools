from gdsctools.anova import GDSC_ANOVA, GDSC_ANOVA_Results
from . import tools



def test_html():
    # same as above, could factorise
    an = GDSC_ANOVA(tools.get_data())
    features = an.features.df
    features = features[features.columns[0:30]]
    an = GDSC_ANOVA(tools.get_data(), features)
    df = an.anova_all()

    r = GDSC_ANOVA_Results(gdsc=an, result=df)
    # long but should cover everthinh.
    r.create_html_pages()

