from gdsctools.anova import ANOVA, ANOVAReport
from . import tools



def test_html():
    # same as above, could factorise
    an = ANOVA(tools.get_data())
    features = an.features.df
    features = features[features.columns[0:30]]
    an = ANOVA(tools.get_data(), features)
    df = an.anova_all()

    r = ANOVAReport(gdsc=an, results=df)
    # long but should cover everthinh.
    try:
        r.create_html_pages()
    except Exception as err:
        raise err
    finally:
        import shutil
        shutil.rmtree('html_gdsc_anova')


