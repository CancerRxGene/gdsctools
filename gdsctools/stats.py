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
from gdsctools.qvalue import QValue


__all__ = ['MultipleTesting', 'cohens', "signed_effects"]


def multiple_correction(pvalues, method='fdr'):
    mt = MultipleTesting(method=method)
    values = mt.get_corrected_pvalues(pvalues, method=None)
    return values


class MultipleTesting(object):
    """This class eases the computation of multiple testing corrections

    The method implemented so far are based on statsmodels or a local
    implementation of **qvalue** method.

    ================    =============================================
    method name         Description
    ================    =============================================
    bonferroni          one-step correction
    sidak               one-step correction
    holm-sidak          step down method using Sidak adjustments
    holm                step down method using Bonferroni adjustments
    simes-hochberg      step up method (independent)
    hommel              close method based on Simes tests (non
                        negative)
    fdr_bh              FDR Benjamini-Hochberg (non-negative)
    fdr_by              FDR Benjamini-Yekutieli (negative)
    fdr_tsbky           FDR 2-stage Benjamini-Krieger-Yekutieli
                        non negative
    frd_tsbh            FDR 2-stage Benjamini-Hochberg'
                        non-negative
    fdr                 same as fdr_bh
    qvalue              see :class:`~gdsctools.qvalue.QValue` class
    ================    =============================================


    .. seealso:: :mod:`gdsctools.qvalue`.

    .. seealso:: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2907892/

    """
    def __init__(self, method=None):
        """.. rubric:: Constructor

        :param method: default to **fdr** that is the FDR Benjamini-Hochberg
            correction.
        """

        #: set of valid methods
        self.valid_methods = ['bonferroni', 'sidak', 'fdr_by',
            'holm-sidak', 'simes-hochberg', 'hommel', 'fdr_bh',
            'fdr_tsbh', 'fdr_tsbky', 'fdr', 'qvalue']

        self._method = 'fdr'
        if method is not None:
            self.method = method
        # parameter of the multiple test (e.g. used if method is bonferroni
        self.alpha = 0.1

    def _get_method(self):
        return self._method
    def _set_method(self, method):
        easydev.check_param_in_list(method, self.valid_methods)
        if method == 'fdr':
            method = 'fdr_bh'
        self._method = method
    method = property(_get_method, _set_method, doc="get/set method")

    def get_corrected_pvalues(self, pvalues, method=None):
        """Return corrected pvalues

        :param list pvalues: list or array of pvalues to correct.
        :param method: use the one defined in the constructor by default
            but can be overwritten here
        """
        if method is not None:
            self.method = method

        pvalues = np.array(pvalues)

        if self.method == 'qvalue':
            qv = QValue(pvalues)
            corrections = qv.qvalue()
            return corrections
        else:
            corrections = multitest.multipletests(pvalues,
               alpha=self.alpha, method=self.method)[1]
            return corrections

    def plot_comparison(self, pvalues, methods=None):
        """Simple plot to compare the pvalues correction methods

        .. plot::
            :include-source:
            :width: 80%

            from gdsctools.stats import MultipleTesting
            mt = MultipleTesting()
            pvalues = [1e-10, 9.5e-2, 2.2e-1, 3.6e-1, 5e-1, 6e-1,8e-1,9.6e-1]
            mt.plot_comparison(pvalues,
                methods=['fdr_bh', 'qvalue', 'bonferroni', 'fdr_tsbh'])

        .. note:: in that example, the qvalue and FDR are identical, but
            this is not true in general.

        """
        if methods is None:
            methods = self.valid_methods

        import pylab
        pylab.clf()
        for method in methods:
            pv = self.get_corrected_pvalues(pvalues, method=method)
            pylab.plot(pvalues, pv, 'o-', label=method.replace("_","\_"))
        pylab.legend(loc='best')
        pylab.ylabel('corrected pvalues')
        pylab.grid()
        pylab.ylim([0, 1.05])


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


def signed_effects(df):
    import numpy as np
    _colname_deltas = 'FEATURE_delta_MEAN_IC50'
    _colname_effect_size = 'FEATURE_IC50_effect_size'
    deltas = df[_colname_deltas]
    effects = df[_colname_effect_size]
    signed_effects = list(np.sign(deltas) * effects)
    return signed_effects
