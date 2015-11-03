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
from statsmodels.stats import multitest
import easydev


__all__ = ['MultipleTesting']


class MultipleTesting(object):
    """


    """
    def __init__(self, method=None):
        self.valid_methods = ['bonferroni', 'sidak',
            'holm-sidak', 'simes-hochberg', 'hommel', 'fdr_bh',
            'fdr_tsbh', 'fdr_tsbky', 'fdr']

        self._method = 'fdr_bh'
        if method is not None:
            self.method = method

    def _get_method(self):
        return self._method
    def _set_method(self, method):
        easydev.check_param_in_list(method, self.valid_methods)
        if method == 'fdr':
            method = 'fdr_bh'
        self._method = method
    method = property(_get_method, _set_method, doc="")

    def get_corrected_pvalues(self, pvalues, method=None):
        if method is not None:
            self.method = method
        corrections = multitest.multipletests(pvalues,
                method=self.method)
        return corrections

    def plot_comparison(self, pvalues):

        import pylab
        pylab.clf()
        for method in self.valid_methods:
            print(method)
            pv = self.get_corrected_pvalues(pvalues, method=method)[1]
            pylab.plot(pvalues, pv, 'o-', label=method.replace("_","\_"))
        pylab.legend()

