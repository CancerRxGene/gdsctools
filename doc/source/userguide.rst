Quick Start
=============

Data format
------------



Example
---------


You can analyse a given drug for a given genomic feature

.. plot::
    :include-source:

    from gdsctools import reader, anova
    r = reader.IC50('ANOVA_input.txt')
    an = anova.GDSC_ANOVA(r.ic50, r.features)
    an.anova_one_drug_one_feature('Drug_1_IC50', 'TP53_mut',
        show_boxplot=True)

Or analyse a given drug across all features:

.. plot::
    :include-source:

    from gdsctools import reader, anova
    r = reader.IC50('ANOVA_input.txt')
    an = anova.GDSC_ANOVA(r.ic50, r.features)
    df = an.anova_one_drug('Drug_1_IC50') # no plots are generated here.

Or analyse a all drugs across all features. This takes a long depending on the
number of drugs and features (30 minutes for 250 drugs and 1000 features):

.. plot::

    from gdsctools import reader, anova
    r = reader.IC50('ANOVA_input.txt')
    an = anova.GDSC_ANOVA(r.ic50, r.features)
    drugs = an.drugs[0:5] + ['Drug_29_IC50']
    dfall = an.anova_all(drugs=drugs)
    an.volcano_plot_one_drug(dfall, 'Drug_29_IC50')


Output and visualisation
==========================


HTML report
==============



Standalone applications
==========================



Notebooks
===============
