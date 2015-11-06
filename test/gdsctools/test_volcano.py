from gdsctools.volcano import VolcanoANOVA
from . import tools
from  gdsctools.anova import ANOVA


def test_volcano_plot():

    an = ANOVA(tools.get_data())
    an.features.df = an.features.df[an.features.df.columns[0:10]]
    an = ANOVA(tools.get_data(), genomic_features=an.features.df)

    results = an.anova_all()

    results.drop('log max.Conc.tested2', axis=1, inplace=True)
    results.drop('log max.Conc.tested', axis=1, inplace=True)
        
    v = VolcanoANOVA(results)
    v.settings.savefig = False
    v.volcano_plot_all_drugs()
    v.volcano_plot_all()
