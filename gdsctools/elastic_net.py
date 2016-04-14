# -*- python -*-
# -*- coding utf-8 -*-
#
#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.comWE HERE>
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
import scipy
import pylab
import numpy as np

from statsmodels.formula.api import OLS

from easydev import Progress, AttrDict

from gdsctools.models import BaseModels

from gdsctools.boxplots import BoxPlots
from gdsctools.settings import ANOVASettings
from gdsctools.anova_results import ANOVAResults

from sklearn.linear_model import enet_path
from sklearn import preprocessing
import sklearn.linear_model
from sklearn import cross_validation
from sklearn.linear_model import ElasticNetCV


from easydev.profiler import do_profile

__all__ = ['ElasticNet']



"""

# IC50
ic50 = IC50("/home/cokelaer/Work/2016/MichaElasticNet/GDSC_web_v5.0/drugRes_OMAUC.csv")
Y = ic50.df['1047'].dropna()

# geonimc features
gf = GenomicFeatures("/home/cokelaer/Work/2016/MichaElasticNet/GDSC_web_v5.0/bem.csv")
gf.df.drop("TISSUE_FACTOR", axis=1, inplace=True)
X = gf.df

X = gf.df.ix[Y.index]



train_obs = Y.iloc[0:595]
XTrain_obs = Y.iloc[595:]

train_feat = X.iloc[0:595]
XTrain_feat = X.iloc[595:]


myalphas = [0.08633579, 0.07866596, 0.07167749, 0.06530986, 0.05950791, 0.05422139, 0.04940451, 0.04501555,  \
0.0410165, 0.0373727, 0.03405262, 0.03102747, 0.02827108, 0.02575955, 0.02347115, 0.02138603,  0.01948616, \
0.01775506, 0.01617775, 0.01474056, 0.01343105, 0.01223788, 0.0111507, 0.0101601, 0.009257504, 0.008435093, \
0.007685742, 0.007002962, 0.006380838, 0.005813982, 0.005297484, 0.00482687, 0.004398064, 0.004007352, \
0.00365135, 0.003326974, 0.003031415, 0.002762113, 0.002516734, 0.002293154, 0.002089437, 0.001903817, \
0.001734687, 0.001580582, 0.001440168, 0.001312227, 0.001195653, 0.001089434,  0.0009926518, 0.0009044673, \
0.000824117, 0.0007509047, 0.0006841964, 0.0006234143, 0.0005680319, 0.0005175695, 0.00047159, 0.0004296953, \
0.0003915223, 0.0003567406, 0.0003250487, 0.0002961723, 0.0002698612, 0.0002458874, 0.0002240435, \
0.0002041401, 0.0001860048, 0.0001694807, 0.0001544245, 0.0001407058, 0.0001282059, 0.0001168165, \
0.0001064388, 9.698307e-05, 8.836736e-05, 8.051705e-05, 7.336414e-05, 6.684667e-05, 6.090819e-05, 5.549728e-05]


"""

class ElasticNet(BaseModels):
    """ElasticNet analysis of the IC50 vs Feature matrices




    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import ElasticNet
        ic = IC50("IC50_v5.csv.gz")
        gf = GenomicFeatures("genomic_features_v5.csv.gz")

        en = ElasticNet(ic50, gf)
        en.elastic_net('1047')


    """
    def __init__(self, ic50, genomic_features=None,
            drug_decode=None, verbose=True, low_memory=True,
            set_media_factor=False):
        """.. rubric:: Constructor

        :param DataFrame IC50: a dataframe with the IC50. Rows should be
            the COSMIC identifiers and columns should be the Drug names
            (or identifiers)
        :param features: another dataframe with rows as in the IC50 matrix
            and columns as features.  The first 3 columns must be named
            specifically to hold tissues, MSI (see format).
        :param drug_decode: a 3 column CSV file with drug's name and targets
            see :mod:`readers` for more information.
        :param verbose: verbosity in "WARNING", "ERROR", "DEBUG", "INFO"

        The attribute :attr:`settings` contains specific settings related
        to the analysis or visulation.
        """
        super(ElasticNet, self).__init__(ic50, genomic_features,
            drug_decode=drug_decode, verbose=verbose, low_memory=low_memory,
            set_media_factor=set_media_factor)

    def _get_one_drug_data(self, name, scale=True):
        Y = self.ic50.df[name]
        Y.dropna(inplace=True)
        X = self.features.df.ix[Y.index].copy()
        X = X.drop('TISSUE_FACTOR', axis=1)

        if scale is True:
            columns = X.columns
            X = preprocessing.scale(X)
            X = pd.DataFrame(X, columns=columns)
        return X, Y

    def elastic_net(self, drug_name, alpha=1, l1_ratio=0.5, n_folds=10,
                    plot=False, tol=1e-3):
        # by default alpha=1,  # l1_ratio=0.5

        # l1_ratio correspond to alpha in glmnet R pacakge
        # while alpha correspond to the lambda parameter in glmnet

        # specifically l1_ratio=1 is the lasso penalty 

        # l1_ratio < 0.01 is not reliable unless sequence of alpha is provided.
        #alpha = 0 correspond to an OLS

        # Get the data for the requested drug
        xscaled, Y = self._get_one_drug_data(drug_name, scale=True)

        # Create a model of elastic net
        # normalise is False because X is scaled
        en = sklearn.linear_model.ElasticNet(alpha=alpha, l1_ratio=l1_ratio,
                                             normalize=False, tol=tol)

        # Create a cross validation set of indices for training and testing
        if n_folds == 1:
            kf = cross_validation.KFold(len(Y), n_folds=2, shuffle=False)
        else:
            kf = cross_validation.KFold(len(Y), n_folds=n_folds, shuffle=False)

        # 
        scores = []
        count = 1
        self.kfold_data = {'x_test':[], 'y_test':[], 'y_train':[], 'x_train':[]}
        for train_index, test_index in kf:
            # Get X training and testing data sets
            X_train = xscaled.iloc[train_index]
            X_test = xscaled.iloc[test_index]

            # Get Y training and testing data sets
            Y_train = Y.iloc[train_index]
            Y_test = Y.iloc[test_index]

            self.kfold_data['x_test'].append(X_test)
            self.kfold_data['x_train'].append(X_train)

            self.kfold_data['y_test'].append(Y_test)
            self.kfold_data['y_train'].append(Y_train)

            en.fit(X_train, Y_train)
            #en.predict(X_test)

            # now compare the prediction with Y_test 
            score = en.score(X_test, Y_test)
            scores.append(score)

            if plot is True:
                pylab.clf()
                pylab.plot(en.predict(X_test), Y_test, 'o')
                #print en.predict(X_test), Y_test
                #print en.predict(X_test), Y_test
                pylab.xlabel("prediction")
                pylab.ylabel("test values")

            if n_folds == 1 and count == 1:
                break
            else:
                count += 1
        self.en = en
        self.kf = kf
        self.X = xscaled
        self.Y = Y

        return scores

    def plot_cindex(self, drug_name, alphas, l1_ratio=0.5, n_folds=10, hold=False):
        # This is longish (300 seconds with 10 folds and 80 alphas
        # for GDSC v5 data sets.
        from dreamtools.core.cindex import cindex

        CI_train = {}
        CI_test = {}
        for c in range(n_folds):
            CI_train[c] = []
            CI_test[c] = []

        from easydev import Progress
        pb = Progress(len(alphas))

        for i, alpha in enumerate(alphas):
            self.elastic_net(drug_name, alpha=alpha, l1_ratio=l1_ratio,
                             n_folds=n_folds)

            # Look at the first fold only
            for kf in range(n_folds):
                x_train = self.kfold_data['x_train'][kf].values
                y_train = self.kfold_data['y_train'][kf].values

                x_test = self.kfold_data['x_test'][kf].values
                y_test = self.kfold_data['y_test'][kf].values

                x_train_pred = self.en.predict(x_train)
                x_test_pred = self.en.predict(x_test)

                CI_test[kf].append(1-cindex(x_test_pred, y_test, [True]*len(y_test)))
                CI_train[kf].append(1-cindex(x_train_pred, y_train, [True] * len(y_train)))
            pb.animate(i)

        mu_train = pd.DataFrame(CI_train).transpose().mean()
        sigma_train = pd.DataFrame(CI_train).transpose().std()

        mu_test = pd.DataFrame(CI_test).transpose().mean()
        sigma_test = pd.DataFrame(CI_test).transpose().std()

        best_alpha = alphas[pd.DataFrame(CI_test).mean(axis=1).argmax()]

        pylab.clf()
        pylab.errorbar(pylab.log(alphas), mu_train, yerr=sigma_train, label="train")
        pylab.errorbar(pylab.log(alphas)+.1, mu_test, yerr=sigma_test, label="test")
        pylab.plot(pylab.log(alphas), mu_train, 'ob')
        pylab.plot(pylab.log(alphas)+.1, mu_train, 'or')
        pylab.legend()
        pylab.axvline(pylab.log(best_alpha), lw=2, color="purple")

        return best_alpha

    def tune_alpha(self, drug_name, alphas=None, N=100, l1_ratio=0.5,
                   n_folds=10, plot=True):
        """

        .. plot::

            an.tune_alpha("1047", N=100, l1_ratio=0.1)
            an.tune_alpha("1047", N=100, l1_ratio=0.01)
            an.tune_alpha("1047", N=100, l1_ratio=0.001)
            an.tune_alpha("1047", N=100, l1_ratio=0.0001)

        29, 34, 52, 1014, 1015, 1024, 1036, 1047, 1061

        """
        # alphas = 10**-linspace(6,1,100)
        if alphas is None:
            alphas = pylab.logspace(-5,0,N)

        all_scores = []
        median_scores = []
        for alpha in alphas:
            scores = self.elastic_net(drug_name, alpha, l1_ratio=l1_ratio,
                                      n_folds=n_folds)
            median_scores.append(np.mean(scores))
            all_scores.append(scores)

        #pylab.plot(pylab.log(alphas), median_scores, '-o')
        df = pd.DataFrame(all_scores)

        maximum = df.mean(axis=1).max()
        alpha_best = alphas[df.mean(axis=1).argmax()]

        if plot is True:
            mu = df.mean(axis=1)
            sigma = df.std(axis=1)
            pylab.clf()
            pylab.errorbar(pylab.log(alphas), mu, yerr=sigma)
            pylab.plot(pylab.log(alphas), mu, 'or')
            pylab.axvline(pylab.log(alpha_best), lw=4, alpha=0.5, color='g')
            pylab.title("Mean scores across alphas")
            pylab.xlabel("alpha")
            pylab.ylabel("mean score")

        return alphas, all_scores, maximum, alpha_best

    def plot_weight(self, drug_name, alpha, l1_ratio=0.5):
        """

        small alphas will have a more stringent effects on the weigths and
        select only some of them setting others to 0. Setting alphas to zero will
        show all weights

        """
        pylab.figure(1)
        pylab.clf()
        self.elastic_net(drug_name, alpha=alpha)
        df = pd.DataFrame({'name': self.X.columns, 'weight': self.en.coef_})
        df = df.set_index("name").sort_values("weight")
        df.plot(kind="bar",  width=1, lw=1, ax=pylab.gca())

        pylab.figure(2)
        pylab.clf()
        df = abs(df)
        df.sort_values("weight").plot(kind="bar", width=1, lw=1,
                                      title='importance plot', ax=pylab.gca())

        return df.sort_values('weight')

    def elastic_net_cv(self, drug_name, l1_ratio=0.5, alphas=None, n_folds=10):

        # Get the data for the requested drug
        xscaled, Y = self._get_one_drug_data(drug_name)

        en = ElasticNetCV(l1_ratio=l1_ratio, alphas=alphas, cv=n_folds)

        encv = en.fit(xscaled, Y)

        self.encv = encv
        print("Best alpha on %s folds : %s" % (n_folds, encv.alpha_))
        #df.sort_values().plot(kind='bar')
        return encv.alpha_

    def enetpath_vs_enet(self, drug_name, alphas=None, l1_ratio=0.5, nfeat=5,
                         max_iter=1000, tol=1e-4, selection="cyclic", fit_intercept=False):
        """
        #if X is not scaled, the enetpath and ElasticNet will give slightly differen results
        #if X is scale using::

        from sklearn import preprocessing
        xscaled = preprocessing.scale(X)
        xscaled = pd.DataFrame(xscaled, columns=X.columns)
        """
        # 1. use enet to loop over alphas and then plot the coefficients along alpha
        # for each feature

        # Get the data for the requested drug
        xscaled, Y = self._get_one_drug_data(drug_name)

        alphas, coefs1, _ = enet_path(xscaled, Y, l1_ratio=l1_ratio, alphas=alphas)
        pylab.figure(1)
        pylab.clf()
        for this in coefs1:
            pylab.plot(pylab.log(alphas), this)

        self.alphas1 = alphas
        self.coefs1 = coefs1

        # Identify the first 5

        # 2. should be equivalenet to using ElasticNet for each alphas
        coefs2 = []
        # if alphas is None, it will be created automatically from enet_path
        for alpha in alphas:
            # to have same results as in enet_path, normalize must be set to
            # False when X is scaled.
            en = sklearn.linear_model.ElasticNet(l1_ratio=l1_ratio, alpha=alpha,
                            max_iter=max_iter, tol=tol, selection=selection,
                            fit_intercept=fit_intercept)
            res = en.fit(xscaled, Y)
            coefs2.append(res.coef_)
        coefs2 = np.array(coefs2).transpose()
        pylab.figure(2)
        pylab.clf()
        for this in coefs2:
            pylab.plot(pylab.log(alphas), this)

        self.coefs2 = coefs2

        #pylab.plot(-pylab.log(res.coef_))

        pylab.figure(3)
        pylab.clf()
        self.df1 = pd.DataFrame(coefs1.transpose(), columns=xscaled.columns)
        self.df2 = pd.DataFrame(coefs2.transpose(), columns=xscaled.columns)

        (self.df1 == 0).sum().plot()
        (self.df2 == 0).sum().plot()


        self.indices1 = (self.df1 == 0).sum().sort_values().ix[0:nfeat]
        self.indices2 = (self.df2 == 0).sum().sort_values().ix[0:nfeat]
        names1 = self.indices1.index
        names2 = self.indices2.index
        print(names2)

    def anova_one_drug(self, drug_id, animate=True, output='object'):
        """Computes ANOVA for a given drug across all features

        :param str drug_id: a valid drug identifier.
        :param animate: shows the progress bar
        :return: a dataframe

        Calls :meth:`anova_one_drug_one_feature` for each feature.
        """
        # some features can be dropped ??

        # drop first and second columns that are made of strings
        # works under python2 but not python 3. Assume that the 2 first
        #columns are the sample name and tissue feature
        # Then, we keep only cases with at least 3 features.
        # MSI could be used but is not like in original R code.
        features = self.features.df.copy()
        # need to skip the FACTOR to keep only features
        shift = self.features.shift

        features = features[features.columns[shift:]]
        # FIXME what about features with less than 3 zeros ?
        mask = features.sum(axis=0) >= 3

        # TODO: MSI, tissues, name must always be kept
        #
        selected_features = features[features.columns[mask]]

        # scan all features for a given drug
        assert drug_id in self.ic50.df.columns
        N = len(selected_features.columns)
        pb = Progress(N, 10)
        res = {}

    def anova_all(self, animate=True, drugs=None):
        """Run all ANOVA tests for all drugs and all features.

        :param drugs: you may select a subset of drugs
        :param animate: shows the progress bar
        :return: an :class:`~gdsctools.anova_results.ANOVAResults`
            instance with the dataframe
            stored in an attribute called **df**

        Loops over all drugs calling :meth:`anova_one_drug` for each
        drug and concatenating all results together. Note that once all
        data are gathered, an extra column containing the FDR corrections
        is added to the dataframe using :meth:`add_pvalues_correction`
        method. An extra column  named "ASSOC_ID" is also added with
        a unique identifer sorted by ascending FDR.

        .. note:: A thorough comparison with version v17 give the same FDR
            results (difference ~1e-6); Note however that the qvalue results
            differ by about 0.3% due to different smoothing in R and Python.
        """
        # drop DRUG where number of IC50 (non-null) is below 5
        # axis=0 is default but we emphasize that sum is over
        # column (i.e. drug
        vv = (self.ic50.df.isnull() == False).sum(axis=0)
        # FIXME: should be in one_drug_one_feature ??
        drug_names = vv.index[vv >= self.settings.minimum_nonna_ic50]

        # if user provided a list of drugs, use them:
        if drugs is not None:
            # todo: check valifity of the drug names
            drug_names = drugs[:]

        pb = Progress(len(drug_names), 1)
        drug_names = list(drug_names)
        pylab.shuffle(drug_names)
        if animate is True:
            pb.animate(0)

        for i, drug_name in enumerate(drug_names):
            if drug_name in self.individual_anova.keys():
                pass
            else:
                res = self.anova_one_drug(drug_name, animate=False,
                                          output='dataframe')
                self.individual_anova[drug_name] = res
            if animate is True:
                pb.animate(i+1)
        print("\n")
        if len(self.individual_anova) == 0:
            return ANOVAResults()

        df = pd.concat(self.individual_anova, ignore_index=True)

        if len(df) == 0:
            return df
        # sort all data by ANOVA p-values
        try:
            df.sort_values('ANOVA_FEATURE_pval', inplace=True)
        except:
            df.sort('ANOVA_FEATURE_pval', inplace=True)

        # all ANOVA have been computed individually for each drug and each
        # feature. Now, we need to compute the multiple testing corrections
        if self.settings.pvalue_correction_level == 'global':
            df = self.add_pvalues_correction(df)

        # insert a unique identifier as first column
        df.insert(0, 'ASSOC_ID', range(1, len(df) + 1))

        self.df = df
        # order the column names as defined in the __init__ method
        df = df[self.column_names]
        df.reset_index(inplace=True, drop=True)

        results = ANOVAResults()
        results.df = df
        results.settings = ANOVASettings(**self.settings)
        return results

    def add_pvalues_correction(self, df, colname='ANOVA_FEATURE_pval'):
        """Add the corrected pvalues column in a dataframe based on pvalues

        The default method (FDR correction) is stored in
        :attr:`settings.pvalue_correction_method` and can be changed to other
        methods (e.g., *qvalue*)

        .. seealso:: :meth:`anova_all`,
            :class:`~gdsctools.stats.MultipleTesting`
        """
        if len(df) == 0:
            return

        # extract pvalues
        data = df[colname].values

        # set the method and compute new pvalues
        self.multiple_testing.method = self.settings.pvalue_correction_method
        new_pvalues = self.multiple_testing.get_corrected_pvalues(data)
        new_pvalues *= 100
        # insert new columns.
        colname = 'ANOVA_FEATURE_FDR'

        try:
            df.insert(len(df.columns), colname, new_pvalues)
        except:
            # replaces it otherwise
            df[colname] = new_pvalues
        return df

    def reset_buffer(self):
        self.individual_anova = {}



