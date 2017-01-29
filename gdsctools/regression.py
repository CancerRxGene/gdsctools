# -*- python -*-
# -*- coding utf-8 -*-
#
#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#  Copyright (c) 2016 - Institut Pasteur
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.com>
#  File author(s): Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
"""Look for IC50 vs and genomic features associations using Regression methods"""
import itertools
import warnings

import pandas as pd
import pylab
import numpy as np

from easydev import Progress

from gdsctools.models import BaseModels
from gdsctools.boxswarm import BoxSwarm

from sklearn.linear_model import enet_path
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import linear_model # must use the module rather than classes to


__all__ = ["Regression", 'GDSCRidge', "GDSCLasso", "GDSCElasticNet", 
           "RegressionCVResults"]


"""book keeping

from statsmodels.formula.api import OLS

if self.settings.regression_method == 'ElasticNet':
    self.data_lm = OLS(odof.Y, df.values).fit_regularized(
                       alpha=self.settings.regression_alpha,
                       L1_wt=self.settings.regression_L1_wt)
elif self.settings.regression_method == 'OLS':
    self.data_lm = OLS(odof.Y, df.values).fit()
elif self.settings.regression_method == 'Ridge':
    self.data_lm = OLS(odof.Y, df.values).fit_regularized(
                       alpha=self.settings.regression_alpha, L1_wt=0)
elif self.settings.regression_method == 'Lasso':
    self.data_lm = OLS(odof.Y, df.values).fit_regularized(
                       alpha=self.settings.regression_alpha, L1_wt=1)

"""


class RegressionCVResults(object):
    """Simple data structure to hold some results of the regression analysis

    - :attr:`model`
    - :attr:`kfold`: number of folds used
    - :attr:`Rp`
    - :attr:`alpha`: best alpha parameter
    - :attr:`ln_alpha`  best alpha parameter (log scale)
    """
    def __init__(self, model, Rp, kfold=None):
        self.model = model
        self.Rp = Rp
        self.kfold = kfold
    def _get_alpha(self):
        return self.model.alpha_
    alpha = property(_get_alpha)
    def _get_ln_alpha(self):
        return pylab.log(self.alpha)
    ln_alpha = property(_get_ln_alpha)
    def _get_coefficients(self):
        return self.model.coef_
    coefficients = property(_get_coefficients)
    def __str__(self):
        txt = "Best alpha on %s folds: %s (%.2f in log scale); Rp=%s" %\
                  (self.kfold, self.alpha, self.ln_alpha, self.Rp)
        return txt


class Regression(BaseModels):
    """Base class for all Regression analysis



    In the :class:`gdsctools.anova.ANOVA` case, the regression is based on the 
    OLS method and is computed for a given drug and a given feature (:term:`ODOF`). 
    Then, the analysis is repeated across all features for a
    given drug (:term:`ODAF`) and finally extended to all drugs (:term:`ADAF`). 
    So, there is one test for each combination of drug and feature.

    Here, all features for a given drug are taken together to perform a
    Regression analysis (:term:`ODAF`). The regression algorithm implemented so 
    far are:

    - Ridge
    - Lasso
    - ElasticNet
    - LassoLars

    Based on tools from the scikit-learn library.

    """
    def __init__(self, ic50, genomic_features=None,
            verbose=False):
        """.. rubric:: Constructor

        :param ic50: an IC50 file
        :param genomic_features: a genomic feature file

        see :ref:`data` for help on the input data formats.

        """
        super(Regression, self).__init__(ic50, genomic_features,
            verbose=verbose, set_media_factor=False)
        self.scale = False

    def _get_one_drug_data(self, name, randomize_Y=False):
        """Returns X and Y for a given drug, dropping NAs

        :param name: drug name
        :param randomize_Y: randomize Y

        - drops NA
        - drops TISSUE_FACTOR
        - drops MSI factor

        """
        Y = self.ic50.df[name]
        Y.dropna(inplace=True)
        X = self.features.df.ix[Y.index].copy()
        try:X = X.drop('TISSUE_FACTOR', axis=1)
        except:pass
        try: X = X.drop('MSI_FACTOR', axis=1)
        except:pass

        if self.scale is True:
            columns = X.columns
            # cast is essential here otherwise ValueError is raised
            X = preprocessing.scale(X.astype(float))
            X = pd.DataFrame(X, columns=columns)

        if randomize_Y:
            Y = Y.copy()
            pylab.shuffle(Y.values)
        return X, Y

    def _fit_model(self, drug_name, model):
        """call fit method of a model given a drug name

        Save the current X, Y, model fitter in _X, _Y and _model attributes
        """
        X, Y = self._get_one_drug_data(drug_name)
        model.fit(X, Y)
        return model

    def plot_importance(self, drug_name, model=None, fontsize=11,
            max_label_length=35, orientation="vertical"):
        """Plot the absolute weights found by a fittd model.


        :param str drug_name:
        :param model: a model
        :param int fontsize: (defaults to 11)
        :param max_label_length: 35 by default
        :param orientation: orientation of the plot (vertical or horizontal)

        :return: the dataframe with the weights (may be empty)

        .. note:: if no weights are different from zeros, no plots are
            created.

        """
        X, Y = self._get_one_drug_data(drug_name)
        if model is None:
            model = self.get_best_model(drug_name)
        model.fit(X, Y)

        df = pd.DataFrame({'name': X.columns, 'weight': model.coef_})
        df = df.set_index("name")

        df = df[df['weight'] != 0]

        if len(df):
            barplot(df, "weight", orientation=orientation, max_label_length=max_label_length,
                fontsize=fontsize)
            # sometimes there is only a few weights close to zero. Set the ylim
            # to 1 if it is below 0.1. 
            if pylab.ylim()[1] < 0.1:
                pylab.ylim([0,1])
            if len(df) < 5:
                pylab.xlim(-5,5)
        return df

    def _print(self, txt):
        if self.verbose:
            print(txt)

    def get_best_model(self, drug_name, kfolds=10, alphas=None, l1_ratio=0.5):
        """Return best model fitted using a CV

        :param drug_name:
        :param kfolds:
        :param alphas:
        :param l1_ratio:

        """
        self._print("Running CV to estimate best alpha.")
        results = self.runCV(drug_name, kfolds=kfolds, alphas=alphas,
                             l1_ratio=l1_ratio)
        best_alpha = results.alpha
        model = self.get_model(alpha=best_alpha)
        self._print("Using alpha=%s." % model.alpha)
        return model

    def plot_weight(self, drug_name, model=None, fontsize=12,
            figsize=(10,7), max_label_length=20, Nmax=40):
        """Plot the elastic net weights

        :param drug_name: the drug identifier
        :param alpha:
        :param l1_ratio:

        Large alpha values will have a more stringent effects on the
        weigths and select only some of them or maybe none. Conversely,
        setting alphas to zero will keep all weights.

        .. plot::
            :include-source:

            from gdsctools import *
            ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
            gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))
            en = GDSCElasticNet(ic, gf)
            model = en.get_model(alpha=0.01)
            en.plot_weight(1047, model=model)

        """
        X, Y = self._get_one_drug_data(drug_name)
        if model is None:
            model = self.get_best_model(drug_name)
        model.fit(X, Y)

        df = pd.DataFrame({'name': X.columns, 'weight': model.coef_})
        df = df.set_index("name").sort_values("weight")
        df = df[df != 0].dropna()

        # split the data keeping only 50 best weights at most
        if len(df) > Nmax:
            # figure out the threshold in absolute value so that the
            # we keep the Nmax strongest weights irrespective of the sign
            threshold = df.abs().sort_values(by="weight").values[-Nmax:][0,0]

            df1 = df.query("weight<=0 and abs(weight)>=@threshold").copy()
            df2 = df.query("weight>=0 and abs(weight)>=@threshold").copy()
        else:
            df1 = df[df.weight<0].copy()
            df2 = df[df.weight>=0].copy()

        df1.index = [this[0:max_label_length] for this in df1.index]
        df2.index = [this[0:max_label_length] for this in df2.index]

        # We also want some symmetry so as many red as blue so that the span
        # of positive and negative is equivalent
        N = len(df2) - len(df1)
        if N > 0:
            # more red LHS than blue RHS
            for i in range(1, N+1):
                label = "_dummy%s" % i
                df1.loc[label, "weight"] = 0
            df1.index = [x if not x.startswith("_dummy") else ""
                         for x in df1.index]
        elif N < 0:
            # more blue RHS than red LHS
            for i in range(1, abs(N)+1):
                label = "_dummy%s" % i
                df2.loc[label, "weight"] = 0
            df2.index = [x if not x.startswith("_dummy") else ""
                         for x in df2.index]
            df2.sort_values(by="weight", ascending=True, inplace=True)


        f, (ax, ax2) = pylab.subplots(1,2, sharey=True, figsize=(10,7))
        ff = pylab.gcf()
        ff.set_facecolor('white')
        self.df1 = df1
        self.df2 = df2

        if len(df1):
            df1.plot(y="weight", kind="bar",  width=1, lw=1, ax=ax,
                color="b", legend=False, fontsize=fontsize, figsize=figsize)
        if len(df2):
            df2.plot(y="weight", kind="bar",  width=1, lw=1, ax=ax2,
                color="r", sharey=True, legend=False, fontsize=fontsize,
                figsize=figsize)
        if len(df1) == 0 and len(df2) == 0:
            pylab.xlim([-5,5])
            pylab.ylim([-1,1])
            pylab.grid(True)
            return df

        # hide the spines between ax and ax2
        ax.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax.yaxis.tick_left()
        ax2.yaxis.tick_right()
        ax2.tick_params(labelleft='off')

        d = 0.02 # diagonal lines
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((1 -d, 1 + d), (-d, +d), **kwargs)
        ax.plot((1 -d, 1 + d), (1-d, 1+d), **kwargs)

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot(( -d,  d), (1-d, 1+d), **kwargs)
        ax2.plot(( -d,  d), (-d, d), **kwargs)
        ax.grid()
        ax2.grid()
        # x0, y0, width_x, width_y
        ax.set_position([0.06,0.3,0.425,0.6])
        ax2.set_position([0.50,0.3,0.425,0.6])
        return df

    def _get_rpearson(self, Y_pred, Y_test):
        self.a = Y_pred
        self.b = Y_test
        if Y_pred.std() == 0:
            Rp = 0
        else:
            Rp = np.corrcoef(Y_pred, Y_test)[0,1]
        if abs(Rp) <1e-10:
            Rp = 0
        return Rp

    def fit(self, drug_name, alpha=1, l1_ratio=0.5, kfolds=10,
                    show=True, tol=1e-3, normalize=False,
                    shuffle=False, perturbation=0.01, randomize_Y=False):
        """Run Elastic Net with a cross validation for one value of alpha

        :param drug_name: the drug to analyse
        :param float alpha: note that theis alpha parameter corresponds to the
            lambda parameter in glmnet R package.
        :param float l1_ratio: This is the lasso penalty parameter.
            Note that in scikit-learn, the l1_ratio correspond
            to the alpha  parameter in glmnet R package. l1_ratio set to 0.5
            means that there is a weight equivalent for the Lasso and Ridge
            effects.
        :param int kfolds: defaults to 10
        :param shuffle: shuffle the indices in the KFold
        :return: kfolds scores for each fold. The score is the pearson
            correlation.

        .. note:: l1_ratio < 0.01 is not reliable unless sequence of
            alpha is provided.

        .. note:: alpha = 0 correspond to an OLS analysis

        """
        assert kfolds > 1, "kfolds must be larger than 1"
        # Get the data for the requested drug
        X, Y = self._get_one_drug_data(drug_name, randomize_Y=randomize_Y)

        # Get a model
        en = self.get_model(alpha=alpha, l1_ratio=l1_ratio,
                            normalize=normalize, tol=tol)

        # Create a cross validation set of indices for training and testing
        kfold = model_selection.KFold(kfolds, shuffle=shuffle)

        # Store the results
        scores = []
        count = 1
        if show is True:
            pylab.clf()

        for train_index, test_index in kfold.split(Y):
            # Get X training and testing data set
            X_train = X.iloc[train_index]
            X_test = X.iloc[test_index]

            # Get Y training and testing data set
            Y_train = Y.iloc[train_index]
            Y_test = Y.iloc[test_index]

            # Fit model on the training set
            en.fit(X_train, Y_train)

            # now compare the prediction with Y_test. This is the coefficient
            # of determination R^2. See scikit learn doc for details.
            Y_pred = en.predict(X_test)

            scores.append(self._get_rpearson(Y_pred, Y_test))

            if show is True:
                N = len(Y_pred)
                import random
                pylab.plot(Y_test, Y_pred + np.append(Y_pred[1:],
                            Y_pred[0]) * perturbation,
                           "ob", alpha=0.5)
                pylab.xlabel("prediction")
                pylab.ylabel("test values")

            if kfolds == 1 and count == 1:
                break
            else:
                count += 1

        if show:
            pylab.title("Prediction on test set (Pearson correlation=%.2f)" %
                        np.mean(scores))
            pylab.xlabel("observed drug response")
            pylab.ylabel("Predicted drug response")
            pylab.grid(True)

        self.en = en
        self.X = X
        self.Y = Y

        return scores

    def runCV(self, drug_name, l1_ratio=0.5, alphas=None, kfolds=10,
                       verbose=True, shuffle=True, randomize_Y=False, **kargs):
        """Perform the Cross validation to get the best alpha parameter.

        :return: an instance of :class:`RegressionCVResults` that contains
            alpha parameter and Pearson correlation value.
        """
        # Get the data for the requested drug
        X, Y = self._get_one_drug_data(drug_name, randomize_Y=randomize_Y)

        if kfolds > len(X):
            kfolds = len(X)

        # Creates a model
        kfold = model_selection.KFold(kfolds, shuffle=shuffle)
        en = self._get_cv_model(l1_ratio=l1_ratio, alphas=alphas, kfold=kfold, **kargs)

        # Fit the model
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = en.fit(X, Y)
        prediction = model.predict(X)

        Rp = self._get_rpearson(prediction, Y)

        res = RegressionCVResults(model, Rp, kfolds)
        if verbose:
            print(res)
        return res

    def tune_alpha(self, drug_name, alphas=None, N=80, l1_ratio=0.5,
                   kfolds=10, show=True, shuffle=False, alpha_range=[-2.8,0.1],
                   randomize_Y=False):
        """Interactive tuning of the model (alpha).

        This is much faster than :meth:`plot_cindex` but much slower than
        :meth:`runCV`.

        .. plot::
            :include-source:

            from gdsctools import *
            ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
            gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))

            en = GDSCElasticNet(ic, gf)

            en.tune_alpha(1047, N=40, l1_ratio=0.1)

        """
        if alphas is None:
            # logspace returns a vector in natural space that guarantees a
            # uniform spacing in a log space (log10 or ln)
            # -2.8 to 0.5 means alpha from 1.58e-3 to 3.16
            # This is equivalent to log(1.58e-3)=-6.45 to log(3.16)=1.15 in ln
            # scale
            a1, a2 = alpha_range
            alphas = pylab.logspace(a1, a2, N)

        # Let us now do a CV across difference alphas
        all_scores = []
        for alpha in alphas:
            scores = self.fit(drug_name, alpha, l1_ratio=l1_ratio,
                              kfolds=kfolds, shuffle=shuffle,
                         randomize_Y=randomize_Y)
            all_scores.append(scores)

        # We can now plot the results that is the mean scores + error enveloppe
        df = pd.DataFrame(all_scores)

        # we also identify the max correlation and corresponding alpha
        maximum = df.mean(axis=1).max()
        alpha_best = alphas[df.mean(axis=1).argmax()]

        if show is True:
            mu = df.mean(axis=1)
            sigma = df.var(axis=1)
            pylab.clf()
            pylab.errorbar(pylab.log(alphas), mu, yerr=sigma, color="gray")
            pylab.plot(pylab.log(alphas), mu, 'or')
            pylab.axvline(pylab.log(alpha_best), lw=4, alpha=0.5, color='g')
            pylab.title("Mean scores (pearson) across alphas for Kfold=%s" % kfolds)
            pylab.xlabel("ln(alpha)")
            pylab.ylabel("mean score (pearson)")
            pylab.grid()

        results = {"alpha_best":alpha_best, "ln_alpha":pylab.log(alpha_best),
            "maximum_Rp":maximum}
        return results
        #return alphas, all_scores, maximum, alpha_best

    def check_randomness(self, drug_name, kfolds=10, N=10,
            progress=False, nbins=40, show=True, **kargs):
        """Compute Bayes factor between NULL model and best model fitted N times


        :param drug_name:
        :param kfolds:
        :param int N: optimise NULL models and real model N times
        :param progress:
        :param nbins:
        :param show:

        Bayes factor::

            S = sum([s>r for s,r in zip(scores, random_scores)])
            proba = S / len(scores)
            bayes_factor = 1. / (1-proba)

        Interpretation for values of the Bayes factor according to Kass 
        and Raftery (1995).

        ============================  ==================
            Interpretation                  B(1,2)
        ============================  ==================
           Very strong support for 1   < 0.0067
                    Strong support 1   0.0067 to 0.05
              Positive support for 1   0.05 to .33
                  Weak support for 1   0.33 to 1
         No support for either model   1
                  Weak support for 2   1 to 3
              Positive support for 2   3 to 20
                Strong support for 2   20 to 150
           Very strong support for 2   > 150
        ============================  ==================

        references: http://www.socsci.uci.edu/~mdlee/LodewyckxEtAl2009.pdf
            http://www.aarondefazio.com/adefazio-bayesfactor-guide.pdf
        """
        scores = []
        pb = Progress(N)
        for i in range(N):
            # Fit a model using CV
            inter_results = self.runCV(drug_name, kfolds=kfolds,
                                       verbose=False, **kargs)
            scores.append(inter_results.Rp)
            if progress:
                pb.animate(i+1)

        random_scores = []
        pb = Progress(N)
        for i in range(N):
            # Fit a model using CV
            inter_results = self.runCV(drug_name, kfolds=kfolds,
                                randomize_Y=True, verbose=False, **kargs)
            random_scores.append(inter_results.Rp)
            if progress:
                pb.animate(i+1)

        from scipy.stats import ttest_ind
        ttest_res = ttest_ind(scores, random_scores)
        results = { "scores": scores,
                    "random_scores": random_scores,
                    "ttest_pval": ttest_res.pvalue}

        # Compute the log of the Bayes factor to avoid underflow as communicated
        # by M.Menden.
        S = sum([s>r for s,r in zip(scores, random_scores)])
        proba = S / len(scores)
        if proba == 1:
            # Set the maximum instead of infinite
            # bayes_factor = np.inf
            bayes_factor = 1. / (1./len(scores))
        else:
            bayes_factor = 1. / (1-proba)
        results['bayes_factor'] = bayes_factor

        M = max(max(scores), max(random_scores)) * 1.2
        m = min(min(scores), min(random_scores)) * 1.2
        if show:
            bins = pylab.linspace(m, M, nbins)
            pylab.clf()
            pylab.hist(scores, bins=bins, color="b", alpha=0.5)
            pylab.hist(random_scores, color="r", alpha=0.5, bins=bins)
            pylab.title("Bayes factor=%(bayes_factor).2f" % results)
            pylab.grid(True)
            pylab.xlabel("Coefficient of correlation Rp")
            pylab.xlabel("#")

        return results

    def dendogram_coefficients(self, stacked=False, show=True, cmap="terrain"):
        """

        shows the coefficient of each optimised model for each drug

        This works for demonstration and small data sets.

        
        """
        drugids = self.drugIds
        from easydev import Progress
        pb = Progress(len(drugids))
        d = {}

        for i, drug_name in enumerate(drugids):
            X, Y = self._get_one_drug_data(drug_name, randomize_Y=False)
            results = self.runCV(drug_name, verbose=False)
            df = pd.DataFrame({'name': X.columns, 'weight': results.coefficients})
            df = df.set_index("name").sort_values("weight")
            d[drug_name] = df.copy()
            pb.animate(i+1)

        # use drugid to keep same order as in the data
        dfall = pd.concat([d[i] for i in drugids], axis=1)
        dfall.columns = drugids

        if show:
            from biokit import heatmap
            h = heatmap.Heatmap(dfall, cmap=cmap)
            h.plot(num=1,colorbar_position="top left")

        if stacked is True:
            dfall = dfall.stack().reset_index()
            dfall.columns = ["feature", "drug", "weight"]
        return dfall

    def boxplot(self, drug_name, model, n=5, minimum_match_per_combo=10,
                bx_vert=True, bx_alpha=0.5, verbose=False, max_label_length=40):
        """Boxplot plot of the most important features.

        :param int n: maximum most important features to be shown (in term
            of coefficient/weight used in the model)
        :param int minimum_match_per_combo: minimum number of alterations
            required for a feature to be shown

        ::

            # assuming there is a drug ID = 29
            r = GDSCLasso()
            model = r.get_best_model(29)
            r.boxplot(29, model)

        In addition, we also show the wild type case where the total number 
        of mutations is 0.


        """
        # The data and fitted model based on provided model
        X, Y = self._get_one_drug_data(drug_name)
        fitted_model = self._fit_model(drug_name, model)
        df = pd.DataFrame({'name': X.columns, 'weight': fitted_model.coef_})
        df = df.set_index("name")

        weights = df.abs().sort_values("weight")[-n:]
        feature_names = weights.index
        features = self.features.df[feature_names].copy() # dont change the data

        sorted_names = features.sum(axis=0).sort_values(ascending=True).index
        features = features[sorted_names].transpose()

        # mutual exclusive sorting
        """for index in range(n):
            new_indices = features.ix[n-index-1].sort_values(ascending=False).index
            features = features.ix[:,new_indices]
        # ignore rows where all features are off

        features = features.ix[:,features.sum(axis=0) >0]
        self._features = features
        self.weights= weights
        pylab.imshow(features.values[::-1], interpolation="None", aspect="auto", cmap="gray_r")
        """

        # barplot
        barplot(weights, "weight")

        #
        indices = {}
        features = features.transpose()
        features['total'] = features.sum(axis=1)

        # This is required
        original_columns = features.columns[:]
        columns = [x.replace("-", "_") for x in features.columns]
        columns = [x.replace("(", "_") for x in columns]
        columns = [x.replace(")", "_") for x in columns]
        columns = [x.replace(",", "_") for x in columns]
        columns = [x.replace(".", "_") for x in columns]
        features.columns = columns
        self._features = features
        # loop over all possible combos
        for n_combo in range(1, len(columns)+1):
            for combo in itertools.combinations(columns, n_combo):
                if "total" in combo:
                    continue
                # using query from pandas, we get indices for which 1 or 2
                # features are on (==1).
                query = " and ".join(["%s==1"] * n_combo)
                query += " and total==%s" % n_combo
                inner_indices = features.query(query % combo).index
                if len(inner_indices) >= minimum_match_per_combo:
                    name = ",".join(combo)
                    indices[name] = inner_indices
                    if verbose:
                        print("Found %s with %s events" % (name,len(inner_indices)))
        # Wild type
        inner_indices = features.query("total==0").index
        indices["Wild Type"] = inner_indices

        _X, Y = self._get_one_drug_data(drug_name)

        names = list(indices.keys())
        means = [Y.ix[indices[name]].dropna().mean() for name in names]
        df = pd.DataFrame({"names":names, "mean": means})
        sorted_names_by_mean = df.sort_values(by="mean")['names']

        data = [Y.ix[indices[name]].dropna() for name in sorted_names_by_mean]

        # the dropna means some subdata/names are now empty
        sorted_names_by_mean = [name for this,name in zip(data, sorted_names_by_mean) if len(this)]
        data = [this for this in data if len(this)]

        for subdata in data:
            subdata.index.name = subdata.index.name[0:max_label_length]

        if len(data):
            bx = BoxSwarm(data, sorted_names_by_mean)
            if bx_vert is False:
                bx.xlabel = "Drug response"
            else:
                bx.ylabel = "Drug response"
            bx.plot(vert=bx_vert, alpha=bx_alpha, widths=0.5)
        return {'weights': weights, "data":data, }


class GDSCRidge(Regression):
    """Same as :class:`Regression`"""
    def __init__(self, ic50, genomic_features=None, verbose=False):
        super(GDSCRidge, self).__init__(ic50, genomic_features,
                                        verbose=verbose)

    def get_model(self,alpha=1, l1_ratio=None, **kargs):
        return linear_model.Ridge(alpha)

    def _get_cv_model(self, alphas=None, kfold=None, l1_ratio=None, **kargs):
        if alphas is None:
            alphas = (0.1, 1.0, 10.0)
        return linear_model.RidgeCV(alphas=alphas, cv=kfold, **kargs)


class GDSCLasso(Regression):
    """See as :class:`Regression`"""
    def __init__(self, ic50, genomic_features=None, verbose=False):
        super(GDSCLasso, self).__init__(ic50, genomic_features,
                                        verbose=verbose)

    def get_model(self, alpha=1, l1_ratio=None, **kargs):
        return linear_model.Lasso(alpha)

    def _get_cv_model(self, alphas=None, kfold=None, l1_ratio=None, **kargs):
        return linear_model.LassoCV(alphas=alphas, cv=kfold, **kargs)


class GDSCLassoLars(Regression):
    """See as :class:`Regression`"""
    def __init__(self, ic50, genomic_features=None, verbose=False):
        super(GDSCLassoLars, self).__init__(ic50, genomic_features,
                                            verbose=verbose)

    def get_model(self, alpha=1, l1_ratio=None, **kargs):
        return linear_model.LassoLarsCV(alpha)

    def _get_cv_model(self, alphas=None, kfold=None, l1_ratio=None, **kargs):
        return linear_model.LassoLarsCV(cv=kfold, **kargs)


class GDSCElasticNet(Regression):
    """Variant of :class:`ANOVA` that handle the association at the drug level
    using an ElasticNet analysis of the IC50 vs Feature matrices.


    As compared to the :class:`GDSCRidge` and :class:`GDSCLasso`

    Here is an example on how to perform the analysis, which is similar to the
    ANOVA API:

    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import GDSCElasticNet, gdsctools_data, IC50, GenomicFeatures
        ic50 = IC50(gdsctools_data("IC50_v5.csv.gz"))
        gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))

        en = GDSCElasticNet(ic50, gf)
        en.enetpath_vs_enet(1047)


    For more information about the input data sets please see
    :class:`~gdsctools.anova.ANOVA`, :mod:`~gdsctools.readers`
    """
    def __init__(self, ic50, genomic_features=None, verbose=False,
            set_media_factor=False):
        """.. rubric:: Constructor

        :param DataFrame IC50: a dataframe with the IC50. Rows should be
            the COSMIC identifiers and columns should be the Drug names
            (or identifiers)
        :param features: another dataframe with rows as in the IC50 matrix
            and columns as features.  The first 3 columns must be named
            specifically to hold tissues, MSI (see format).
        :param verbose: verbosity in "WARNING", "ERROR", "DEBUG", "INFO"

        """
        super(GDSCElasticNet, self).__init__(ic50, genomic_features,
                                             verbose=verbose)

    def get_model(self,  alpha=1, l1_ratio=0.5, tol=1e-3, normalize=False,
                  **kargs):
        return linear_model.ElasticNet(alpha=alpha, l1_ratio=l1_ratio,
                            normalize=normalize, tol=tol)

    def _get_cv_model(self, l1_ratio=0.5, alphas=None, kfold=None, **kargs):
        return linear_model.ElasticNetCV(l1_ratio=l1_ratio, alphas=alphas, cv=kfold, **kargs)

    def plot_cindex(self, drug_name, alphas, l1_ratio=0.5, kfolds=10,
            hold=False):
        """Tune alpha parameter using concordance index


        This is longish and performs the following task. For a set of alpha
        (list), run the elastic net analysis for a given **l1_ratio** with
        **kfolds**. For each alpha, get the CIndex and find the CINdex for
        which the errors are minimum.

        .. warning:: this is a bit longish (300 seconds for 10 folds
            and 80 alphas) on GDSCv5 data set.
        """
        from dreamtools.core.cindex import cindex

        CI_train = {}
        CI_test = {}
        for c in range(kfolds):
            CI_train[c] = []
            CI_test[c] = []

        pb = Progress(len(alphas))

        for i, alpha in enumerate(alphas):
            self.fit(drug_name, alpha=alpha, l1_ratio=l1_ratio,
                             kfolds=kfolds)

            # Look at the results and store cindex
            for kf in range(kfolds):
                x_train = self.kfold_data['x_train'][kf].values
                y_train = self.kfold_data['y_train'][kf].values

                x_test = self.kfold_data['x_test'][kf].values
                y_test = self.kfold_data['y_test'][kf].values

                x_train_pred = self.en.predict(x_train)
                x_test_pred = self.en.predict(x_test)

                CI_test[kf].append(1-cindex(x_test_pred, y_test,
                    [True]*len(y_test)))
                CI_train[kf].append(1-cindex(x_train_pred, y_train,
                    [True] * len(y_train)))
            pb.animate(i)

        mu_train = pd.DataFrame(CI_train).transpose().mean()
        sigma_train = pd.DataFrame(CI_train).transpose().std()

        mu_test = pd.DataFrame(CI_test).transpose().mean()
        sigma_test = pd.DataFrame(CI_test).transpose().std()

        best_alpha = alphas[pd.DataFrame(CI_test).mean(axis=1).argmax()]

        pylab.clf()
        pylab.errorbar(pylab.log(alphas), mu_train, yerr=sigma_train,
                label="train")
        pylab.errorbar(pylab.log(alphas)+.1, mu_test, yerr=sigma_test,
                label="test")
        pylab.plot(pylab.log(alphas), mu_train, 'ob')
        pylab.plot(pylab.log(alphas)+.1, mu_train, 'or')
        pylab.legend()
        pylab.axvline(pylab.log(best_alpha), lw=2, color="purple")

        return best_alpha

    def enetpath_vs_enet(self, drug_name, alphas=None, l1_ratio=0.5, nfeat=5,
                         max_iter=1000, tol=1e-4, selection="cyclic",
                         fit_intercept=False):
        """
        #if X is not scaled, the enetpath and ElasticNet will give slightly
        # different results
        #if X is scale using::

            from sklearn import preprocessing
            xscaled = preprocessing.scale(X)
            xscaled = pd.DataFrame(xscaled, columns=X.columns)
        """
        # 1. use enet to loop over alphas and then plot the coefficients
        # along alpha for each feature

        # Get the data for the requested drug
        X, Y = self._get_one_drug_data(drug_name)

        # Run elasticnet for a bunch of alphas to get the coefficients
        alphas, coeffs, _ = enet_path(X, Y, l1_ratio=l1_ratio, alphas=alphas)

        # estimate the best alpha for later
        best_alpha = self.runCV(drug_name, verbose=False).alpha

        # The final data, coeffs is sorted by coefficients on the smallest alpha
        coeffs = pd.DataFrame(coeffs, index=list(X.columns))
        N = len(alphas)
        coeffs.sort_values(by=N-1, inplace=True)
        results = {"alphas": alphas, "coeffs": coeffs, "best_alpha": best_alpha}

        # the Viewer
        self._plot_enet(results)

        return results

    def _plot_enet(self, data, numfig=1, labels_to_show=5, fontsize=10):
        """Plot the enetpath data

        :param coeffs: a dataframe. The rows are the features. The columns are
            the alpha. The matrix contains the coefficients
        :param coeffs: the alphas. Must match the columns of the coeffs
            dataframe
        :param labels_to_show: number of labels to add on the figure

        """
        pylab.figure(numfig)
        pylab.clf()
        for this in data['coeffs'].iterrows():
            pylab.plot(pylab.log(data['alphas']), this[1], color="gray")
        pylab.axvline(pylab.log(data['best_alpha']), color="r", lw=2)
        pylab.xlabel("ln alpha")
        pylab.ylabel("Coefficients")
        pylab.grid(True)
        # The data is sorted with last rows having largest coeff and last
        # columns the smallest

        if len(data) < labels_to_show:
            labels_to_show = len(data)

        try:
            # recent python3/pandas version
            for text,v in data['coeffs'].iloc[-labels_to_show:,-1].items():
                pylab.text(pylab.log(min(data["alphas"])),v,text, fontsize=fontsize)
        except:
            # old python2/pandas
            for text,v in data['coeffs'].iloc[-labels_to_show:,-1].iteritems():
                pylab.text(pylab.log(min(data["alphas"])),v,text, fontsize=fontsize)


def barplot(df, colname, color="sign",colors=("b", "r"),
            orientation="vertical", Nmax=50, max_label_length=20,
            fontsize=12):

    assert orientation in ["vertical", "horizontal"]
    if orientation=="vertical":
        kind = "bar"
        figsize = (10, 7)
    else:
        kind = "barh"
        figsize = (8, 10)
    df = df.copy()

    pylab.figure(figsize=figsize)
    df.loc[:,'sign'] = df[colname]>0
    df[colname] = abs(df[colname])
    df.index = [this[0:max_label_length] for this in df.index]

    df = df.sort_values(colname, ascending=True)
    df.loc[df['sign'] == True, 'sign'] = 'r'
    df.loc[df['sign'] == False, 'sign'] = 'b'

    colors = "".join(df['sign'])

    if len(df) < Nmax:
        df.plot(y=colname, kind=kind, color=colors, width=1, lw=1,
            title='Importance plot', ax=pylab.gca(),
            fontsize=fontsize, figsize=figsize)
    else:
        df.iloc[-Nmax:].plot(y=colname, kind=kind, color=colors[-50:],
            width=1, lw=1,
            title='Importance plot', ax=pylab.gca(),
            fontsize=fontsize, figsize=figsize)

    import matplotlib.patches as mpatches
    red_patch = mpatches.Patch(color='red', label='positive weights')
    blue_patch = mpatches.Patch(color='blue', label="negative weights")

    ax = pylab.gca()
    if orientation=="horizontal":
        ax.set_position([0.3,0.1,0.6,0.8])
        pylab.legend(handles=[red_patch, blue_patch], loc='lower right')
    elif orientation=="vertical":
        ax.set_position([0.1,0.3,0.8,0.6])
        pylab.legend(handles=[red_patch, blue_patch], loc='upper left')
    pylab.grid()
    pylab.ylabel("Absolute weights", fontsize=fontsize)

