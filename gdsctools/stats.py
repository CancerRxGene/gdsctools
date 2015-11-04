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
import numpy as np

__all__ = ['MultipleTesting', 'cohens']


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


def cohens(x, y):
    r"""Effect size metric through Cohen's *d* metric

    :param x: first vector
    :param y: second vector
    :return: absolute effect size value

    The Cohen's effect size *d* is defined as the difference
    between two means divided by a standard deviation of the data.

    .. math::

        d = \frac{\bar{x}_1 - \bar{x}_2}{s}

    For two independent samples, the *pooled standard deviation* is used
    instead, which is defined as:

    .. math::

        s = \sqrt{  \frac{(n_1-1)s_1^2 + (n_2-1)s_2^2}{n_1+n_2-2} }


    A Cohen's *d* is frequently used in estimating sample sizes for
    statistical testing: a lower *d* value indicates the necessity of
    larger sample sizes, and vice versa.

    .. note:: we return the absolute value

    :references: https://en.wikipedia.org/wiki/Effect_size

    """
    x = np.array(x)
    y = np.array(y)

    Nx = len(x) - 1.  # note the dot to cast to float
    Ny = len(y) - 1.
    # mean difference:
    md = np.abs(x.mean() - y.mean())
    # here, we want same as in R that is unbiased variance
    # so we use ddof = 1
    xv = x.var(ddof=1)
    yv = y.var(ddof=1)
    csd = Nx * xv + Ny * yv
    csd /= Nx + Ny  # make sure this is float
    csd = np.sqrt(csd)

    return md / csd


def glass(x, y):
    r"""Return Effect size through Glass :math:`\Delta` estimator

    :param x: first sample
    :param y: second sample
    :return: 2 values (one or each sample)

    The Glass effect size is computed as 

    .. math::  


        \Delta = \frac{\bar{x}_1-\bar{x}_2}{\sigma_i}

    .. note:: the standard deviation is the unbiased one (divided by N-1)

    where :math:`\sigma` is the standard deviation of either group

    """
    x = np.array(x)
    y = np.array(y)

    # mean difference: 
    md = np.abs(x.mean() - y.mean())

    # here, we want same as in R that is unbiased variance
    # so we use ddof = 1
    g1 = md / x.std(ddof=1)
    g2 = md / y.std(ddof=1)

    return g1, g2
