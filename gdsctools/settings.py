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

    When you change a value, you can check its validity by calling the
    :meth:`check`  method. This is not very thorough right now but should
    help.

    Finally, the method :meth:`to_html` creates an HTML table that can
    be added to HTML report.

    .. note:: **for developers** a key can be changed or accessed to as if
       it was an attribute. This prevents some functionalities (such as copy()
       or property) to be used effectively normaly hence the creation of the
       :meth:`check` method to check validity of the values rather than
       using properties.

    Here are the current values used:

    ======================= ========= =========================================
    Name                    Default   Description
    ======================= ========= =========================================
    includeMSI_factor       True      Include MSI in the regression
    featFactorPopulationTh  3         Discard association where a
                                      genomic feature has less than 3
                                      positives or 3 negatives values
                                      (e.g., 0, 1 or 2)
    MSIfactorPopulationTh   2         Discard association where a MSI
                                      count has less than 2 positives
                                      or 2 negatives values (e.g., 0,
                                      or 1).
    analysis_type           PANCAN    Type of analysis. PANCAN means
                                      use all data. Otherwise, you must
                                      provide a valid tissue name to
                                      be found in the Genomic Feature
                                      data set.
    pval_correction_method  fdr       Type of p-values correction
                                      method used. Could be *fdr*,
                                      *qvalue*  or one accepted
                                      by
                                      :class:`~gdsctools.stats.MultipleTesting`
    equal_var_ttest         True      Assume equal variance in the
                                      t-test
    minimum_nonna_ic50      6         Minimum number of IC50 required
                                      to perform an analysis for a
                                      given drug.
    fontsize                20        Used in some plots for labels
    FDR_threshold           25        FDR threshold used in volcano
                                      plot and significant hits
    pvalue_threshold        np.inf    Used to select significant hits
                                      see
                                      :class:`~gdsctools.anova.ANOVAReport`
    directory               html_gdsc Directory where images and HTML documents
                                      are saved.
    savefig                 False     Save the figure or not (PNG format)
    effect_threshold        0         Used in the volcano plot. See
                                      :class:`~gdsctools.volcano.VolcanoPlot`

    low_memory              False     Faster (20%) if set to false but
                                      uses about 1Gb per run
    ======================= ========= =========================================

    """
    def __init__(self, **kargs):
        super(ANOVASettings, self).__init__(**kargs)

        ## ANALYSIS ---------------------------------
        # include MSI as a co-factor
        self.includeMSI_factor = True
        # number of positive samples required to perform the test
        self.featFactorPopulationTh = 3
        # How many MSI samples must be present to perform the test
        self.MSIfactorPopulationTh = 2
        self.analysis_type = 'PANCAN'
        self.pval_correction_method = 'fdr'   # or qvalue
        self.equal_var_ttest = True
        self.minimum_nonna_ic50 = 6

        self.low_memory = False # False means this will be 10-20 faster

        # Visualisation and HTML related ---------------------
        self.fontsize = 20
        self.FDR_threshold = 25
        self.pvalue_threshold = np.inf
        self.directory = 'html_gdsc_anova'
        self.savefig = False
        self.effect_threshold = 0 # use in volcano

    def check(self):
        """Checks the values of the parameters"""
        inrange = easydev.check_range
        inlist = easydev.check_param_in_list
        # check validity of the settings
        inlist(self.includeMSI_factor, [False, True], 'MSI')
        inrange(self.featFactorPopulationTh, 0, np.inf)
        inrange(self.MSIfactorPopulationTh, 0, np.inf)

        # all those methods are from statsmodels.stats.multitest.multipletests
        inlist(self.pval_correction_method, ['bonferroni', 'sidak',
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
        if self.includeMSI_factor is False:
            assert self.analysis_type != 'PANCAN', \
                'If MSI factor is not included, the analysis must be cancer'+\
                ' specific (i.e., a tissue must be set.'

    def to_html(self):
        """Convert the sets of parameters into a nice HTML table"""
        settings = pd.DataFrame(self, index=[0]).transpose()
        settings.reset_index(inplace=True)
        settings.columns = ['name', 'value']
        html = settings.to_html(header=True, index=False)
        return html

    def __str__(self):
        txt = ''
        for k,v in self.items():
            txt += '- %s: %s\n' % (k,v)
        return txt

