# coding=utf-8
# -*- python -*-
#
#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.com>
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
"""Code related to the ANOVA analysis to find associations between drug IC50s
and genomic features"""
import pandas as pd
import numpy as np
from easydev import AttrDict
import easydev
from gdsctools import version


__all__ = ['ANOVASettings']

class ANOVASettings(AttrDict):
    """All settings used in :class:`gdsctools.anova.ANOVA` analysis

    This class behaves as a dictionary but values for a given
    key (setting) can be accessed and changed easily like an
    attribute:

    ::

        >>> from gdsctools import ANOVASettings
        >>> s = ANOVASettings()
        >>> s.FDR_threshold
        25
        >>> s.FDR_threshold = 20

    It is the responsability of the users or developers to check the validity
    of the settings by calling the :meth:`check` method. Note, however, that
    this method does not perform exhaustive checks.

    Finally, the method :meth:`to_html` creates an HTML table that can
    be added to an HTML report.

    .. note:: **for developers** a key can be changed or accessed to as if
       it was an attribute. This prevents some functionalities (such as copy()
       or property) to be used normaly hence the creation of the
       :meth:`check` method to check validity of the values rather than
       using properties.

    Here are the current values available.


    ========================= ================ ========================================
    Name                      Default          Description
    ========================= ================ ========================================
    include_MSI_factor        True             Include MSI in the regression
    feature_factor_threshold  3                Discard association where a
                                               genomic feature has less than 3
                                               positives or 3 negatives values
                                               (e.g., 0, 1 or 2)
    MSI_factor_threshold      2                Discard association where a MSI
                                               count has less than 2 positives
                                               or 2 negatives values (e.g., 0,
                                               or 1).
    analysis_type             PANCAN           Type of analysis. PANCAN means
                                               use all data. Otherwise, you must
                                               provide a valid tissue name to
                                               be found in the Genomic Feature
                                               data set.
    pvalue_correction_method    fdr            Type of p-values correction
                                               method used. Could be *fdr*,
                                               *qvalue*  or one accepted
                                               by
                                               :class:`~gdsctools.stats.MultipleTesting`

    pvalue_correction_level   global           Apply pvalue correction 
                                               globally. Can also be set to 
                                               'drug_level' to apply 
                                               corrections at drug level
                                               only.
    equal_var_ttest           True             Assume equal variance in the
                                               t-test
    minimum_nonna_ic50        6                Minimum number of IC50 required
                                               to perform an analysis for a
                                               given drug.
    fontsize                  25               Used in some plots for labels
    FDR_threshold             25               FDR threshold used in volcano
                                               plot and significant hits
    pvalue_threshold          0.001            Used to select significant hits see
                                               :class:`~gdsctools.anova.ANOVAReport`
    directory                 html_gdsc_anova  Directory where images and HTML
                                               documents are saved.
    savefig                   False            Save the figure or not (PNG format)
    effect_threshold          0                Used in the volcano plot. See
                                               :class:`~gdsctools.volcano.VolcanoPlot`
    ========================= ================ ========================================


    There are parameters dedicated to the regression method. Note that only
    regression_formula is effective right now.

    ======================= ========= =========================================
    Name                    Default   Description
    ======================= ========= =========================================
    regression_method       OLS       Regression method amongst OLS. NOT USED
                                      YET.
    regression_alpha        0.01      Fraction of penalty included. If 0,
                                      equivalent to OLS. NOT USED YET.
    regression_L1_wt        0.5       Fraction of the penalty given to L1
                                      penalty term. Must be between 0 and 1.
                                      If 0, equivalent to Ridge. If 1
                                      equivalent to Lasso. NOT USED YET.
    regression_formula      auto      if auto, use standard regression from
                                      GDSCTools (see link_formula_)
                                      otherwise any valid regression formula
                                      can be used.
    ======================= ========= =========================================

    .. seealso:: :ref:`settings` or
        gdsctools.readthedocs.org/en/master/settings.html#filtering
        decrease the number of significant hits.


    .. _link_formula: http://gdsctools.readthedocs.io/en/master/anova_parttwo.html#regression-analysis
    """
    def __init__(self, **kargs):
        super(ANOVASettings, self).__init__(**kargs)

        ## ANALYSIS ---------------------------------
        # include MSI as a co-factor
        self.include_MSI_factor = True
        # number of positive samples required to perform the test
        self.feature_factor_threshold = 3
        # How many MSI samples must be present to perform the test
        self.MSI_factor_threshold = 2
        self.include_media_factor = False

        self.analysis_type = 'PANCAN'
        self.pvalue_correction_method = 'fdr'   # or qvalue
        self.pvalue_correction_level = 'global'   # or qvalue
        self.equal_var_ttest = True
        self.minimum_nonna_ic50 = 6

        # Visualisation and HTML related ---------------------
        self.fontsize = 25
        self.FDR_threshold = 25
        self.pvalue_threshold = 0.001
        self.directory = 'html_gdsc_anova'
        self.savefig = False
        self.effect_threshold = 0 # use in volcano
        self.volcano_additional_FDR_lines = [0.01, 0.1, 10]
        self.volcano_FDR_interpolation = True

        # ----------------------- regression related
        self.regression_method = 'OLS' # can be ElasticNet, LAsso, Ridge
        self.regression_alpha = 0.01
        self.regression_L1_wt = 0.5

        self.regression_formula = "auto"

        # uses statsmodels package
        # The fraction of the penalty given to the L1 penalty term. Must be
        # between 0 and 1 (inclusive). If 0, the fit is ridge regression. If
        # 1, the fit is the lasso.
        self.version = version
        self.animate = True

        for k, v in kargs.items():
            self[k] = v

    def check(self):
        """Checks the values of the parameters

        This may not be exhaustive. Right now, checks

         - MSI factor is boolean.
         - Regression.method in OLS/Ridge/Lasso/ElasticNet
         - FDR thresohld in [0,1]
         - pvalues_threshold in [0,inf[
         - effect_threshold in [0,inf[
         - pvalue_correction_method
         - etc
        """
        inrange = easydev.check_range
        inlist = easydev.check_param_in_list
        # check validity of the settings
        inlist(self.include_MSI_factor, [False, True], 'MSI')
        inrange(self.feature_factor_threshold, 0, np.inf)
        inrange(self.MSI_factor_threshold, 0, np.inf)

        # all those methods are from statsmodels.stats.multitest.multipletests
        inlist(self.pvalue_correction_method, ['bonferroni', 'sidak',
            'holm-sidak', 'simes-hochberg', 'hommel', 'fdr_bh',
            'fdr_tsbj', 'fdr_tskby', 'fdr'],
            'pvalue correction method')
        inlist(self.equal_var_ttest, [True, False], 'equal_var_ttest')
        inrange(self.minimum_nonna_ic50, 0, np.inf)
        inrange(self.FDR_threshold, 0, 100)
        inrange(self.pvalue_threshold, 0, np.inf)
        inrange(self.effect_threshold, 0, np.inf)

        # for now, if MSI is False, this cannot be a PANCAN analysis
        # but a cancer specific analysis
        if self.include_MSI_factor is False:
            assert self.analysis_type != 'PANCAN', \
                'If MSI factor is not included, the analysis must be cancer'+\
                ' specific (i.e., a tissue must be set.'

        valid_reg_meth = ['OLS', 'ElasticNet', 'Lasso', 'Ridge']
        inlist(self.regression_method, valid_reg_meth)

        inlist(self.pvalue_correction_level, ['global', 'drug_level'])

    def to_html(self):
        """Convert the sets of parameters into a nice HTML table"""
        data = self.copy()
        data['volcano_additional_FDR_lines'] = \
                str(data['volcano_additional_FDR_lines'])

        settings = pd.DataFrame(data, index=[0]).transpose()

        settings.reset_index(inplace=True)
        settings.columns = ['name', 'value']

        html = settings.to_html(header=True, index=False)
        return html

    def __str__(self):
        txt = ''
        for k in sorted(self.keys()):
            txt += '- %s: %s\n' % (k, self[k])
        return txt
