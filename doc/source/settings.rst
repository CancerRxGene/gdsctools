.. _settings:

Settings
===========

Overview
-----------

When using the :class:`~gdsctools.anova.ANOVA` instance or you wish to create an
HTML report (see :ref:`html`), all tunable settings are accessible from an
attribute called :attr:`settings`::

    from gdsctools import ANOVA, ic50_test
    gdsc = ANOVA(ic50_test)
    gdsc.settings

This :attr:`settings` attribute is an instance of :class:`gdsctools.settings.ANOVASettings`. Here are some of the settings that can be changed before an analysis.


:regression related:

- includeMSI_factor: if true (default) include :term:`MSI` in the regression
- analysis_type:  PANCAN means use all data. Otherwise, cancer specific given a valid tissue name.

:filtering:

- featFactorPopulationTh: across cell lines, a feature must have at least  3 positives and 3 negatives features (e.g., 0, 1 or 2)
- MSIfactorPopulationTh: across cell lines, MSI count (neg and pos) must be at least 2 (e.g., 0, or 1).
- minimum_nonna_ic50 :  Minimum number of IC50 required to perform an analysis
  for a given drug (at least 6).

:multiple testing correction:

- pval_correction_method: type of p-values correction method used (default to
  BH)
- FDR_threshold   FDR threshold used in volcano plot and significant hits


:volcano plots:

- pvalue_threshold    np.inf  Used to select significant hits see ANOVAReport
- effect_threshold    0   Used in the volcano plot. See VolcanoPlot

:others:

- equal_var_ttest: Assume equal variance in the t-test
- low_memory  False   Faster (20%) if set to false but uses about 1Gb per run


.. note:: Some settings will be set automatically when calling some functions.
    For instance, if you call :meth:`anova.ANOVA.set_cancer_type` to a single
    tissue, then the analysis_type will be set to the tissue's name. If there 
    are not enough positive or negative MSI, the MSI factor will ignored.
