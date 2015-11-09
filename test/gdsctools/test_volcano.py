from gdsctools.volcano import VolcanoANOVA
from . import tools
from  gdsctools.anova import ANOVA


def test_volcano_plot():

    an = ANOVA(tools.get_data())
    an.features.df = an.features.df[an.features.df.columns[0:10]]
    an = ANOVA(tools.get_data(), genomic_features=an.features.df)

    results = an.anova_all()

    # try the constructors
    v = VolcanoANOVA(results.df)
    v = VolcanoANOVA(results)

    # the selector metho
    v.df = v.selector(v.df)

    v.settings.savefig = False

    # some of the plotting
    v.volcano_plot_all_drugs()
    v.volcano_plot_all_features()
    v.volcano_plot_all()


    v._get_fdr_from_pvalue_interp(1e-10)
    v._get_pvalue_from_fdr(50)
    v._get_pvalue_from_fdr([50,60])
