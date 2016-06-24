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
"""Look for IC50 vs and genomic features associations using Elastic Net"""
import pandas as pd
import pylab
import numpy as np


from easydev import Progress, AttrDict

from gdsctools.models import BaseModels

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
    """Variant of :class:`ANOVA` that handle the association at the drug level
    using an ElasticNet analysis of the IC50 vs Feature matrices.


    In the ANOVA case, regression using OLS are computed for a given drug and a
    given feature. Then, this analysis is repeated across all features for a
    given drug and finally extended to all drugs. However, there is one test
    for each combination of drug and feature.

    Here, all feature for a given drug are taken together to perform an Elastic
    Net analysis.


    Here is an example on how to perform the analysis, which is similar to the
    ANOVA API:

    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import ElasticNet, gdsctools_data, IC50, GenomicFeatures
        ic50 = IC50(gdsctools_data("IC50_v5.csv.gz"))
        gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))

        en = ElasticNet(ic50, gf)
        en.elastic_net(1047, alpha=0.01, show=True)


    For more information about the input data sets please see
    :class:`~gdsctools.anova.ANOVA`, :mod:`~gdsctools.readers`
    """
    def __init__(self, ic50, genomic_features=None,
            drug_decode=None, verbose=False,
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

        """
        super(ElasticNet, self).__init__(ic50, genomic_features,
            drug_decode=drug_decode, verbose=verbose, 
            set_media_factor=set_media_factor)

    def _get_one_drug_data(self, name, scale=True):
        Y = self.ic50.df[name]
        Y.dropna(inplace=True)
        X = self.features.df.ix[Y.index].copy()
        X = X.drop('TISSUE_FACTOR', axis=1)

        if scale is True:
            columns = X.columns
            # cast is essential here otherwise ValueError is raised
            X = preprocessing.scale(X.astype(float))
            X = pd.DataFrame(X, columns=columns)
        return X, Y

    def elastic_net(self, drug_name, alpha=1, l1_ratio=0.5, n_folds=10,
                    show=False, tol=1e-3):
        """Run Elastic Net

        :param drug_name: the drug to analyse
        :param float alpha: note that theis alpha parameter corresponds to the
            lambda parameter in glmnet R package.
        :param float l1_ratio: This is the lasso penalty parameter. 
            Note that in scikit-learn, the l1_ratio correspond 
            to the alpha  parameter in glmnet R package. l1_ratio set to 0.5
            means that there is a weight equivalent for the Lasso and Ridge 
            effects.
        :param int n_folds: defaults to 10 

        .. note:: l1_ratio < 0.01 is not reliable unless sequence of 
            alpha is provided.

        .. note:: alpha = 0 correspond to an OLS analysis

        """
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

        # Store the results 
        scores = []
        count = 1
        self.kfold_data = {'x_test':[], 'y_test':[], 
                'y_train':[], 'x_train':[]}
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

            if show is True:
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

    def plot_cindex(self, drug_name, alphas, l1_ratio=0.5, n_folds=10, 
            hold=False):
        """Tune alpha parameter using concordance index


        This is longish and performs the following task. For a set of alpha
        (list), run the elastic net analysis for a given **l1_ratio** with 
        **n_folds**. For each alpha, get the CIndex and find the CINdex for
        which the errors are minimum. 

        .. warning:: this is a bit longish (300 seconds for 10 folds 
            and 80 alphas) on GDSCv5 data set.
        """
        from dreamtools.core.cindex import cindex

        CI_train = {}
        CI_test = {}
        for c in range(n_folds):
            CI_train[c] = []
            CI_test[c] = []

        pb = Progress(len(alphas))

        for i, alpha in enumerate(alphas):
            self.elastic_net(drug_name, alpha=alpha, l1_ratio=l1_ratio,
                             n_folds=n_folds)

            # Look at the results and store cindex
            for kf in range(n_folds):
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

    def tune_alpha(self, drug_name, alphas=None, N=100, l1_ratio=0.5,
                   n_folds=10, show=True):
        """Another method to tun the alpha parameter.

        This is much faster than :meth:`plot_cindex`.

        .. plot::
            :include-source:

            from gdsctools import *
            ic = IC50(gdsctools_data("IC50_v5.csv.gz"))
            gf = GenomicFeatures(gdsctools_data("genomic_features_v5.csv.gz"))

            en = ElasticNet(ic, gf)

            en.tune_alpha(1047, N=40, l1_ratio=0.1)


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

        if show is True:
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

    def plot_weight(self, drug_name, alpha, l1_ratio=0.5, fontsize=12,
            figsize=(10,7), max_label_length=20):
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
            en = ElasticNet(ic, gf)
            en.plot_weight(1047, 0.01, l1_ratio=0.5)


        """
        self.elastic_net(drug_name, alpha=alpha)
        df = pd.DataFrame({'name': self.X.columns, 'weight': self.en.coef_})
        df = df.set_index("name").sort_values("weight")

        if len(df)<50:
            pylab.figure(1)
            pylab.clf()
            df.plot(kind="bar",  width=1, lw=1, ax=pylab.gca(),
                fontsize=fontsize, figsize=figsize)
        else:
            # split the data
            df1 = df.iloc[0:25]
            df2 = df.iloc[-25:]
            df1.index = [this[0:max_label_length] for this in df1.index]
            df2.index = [this[0:max_label_length] for this in df2.index]

            f, (ax, ax2) = pylab.subplots(1,2, sharey=True, figsize=(10,7))
            ff = pylab.gcf()
            ff.set_facecolor('white')
            df1.plot(y="weight", kind="bar",  width=1, lw=1, ax=ax,
                color="b", legend=False, fontsize=fontsize, figsize=figsize)
            df2.plot(y="weight", kind="bar",  width=1, lw=1, ax=ax2,
                color="r", sharey=True, legend=False, fontsize=fontsize,
                figsize=figsize)

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

    def plot_importance(self, drug_name, alpha, l1_ratio=0.5, fontsize=11,
            max_label_length=35, orientation="vertical"):
        assert orientation in ["vertical", "horizontal"] 
        if orientation=="vertical":
            kind = "bar"
            figsize = (10, 7)
        else:
            kind = "barh"
            figsize = (8, 10)

        self.elastic_net(drug_name, alpha=alpha)
        df = pd.DataFrame({'name': self.X.columns, 'weight': self.en.coef_})
        df = df.set_index("name")#.sort_values("weight")

        pylab.figure(figsize=figsize)

        df.loc[:,'sign'] = df['weight']>0

        df["weight"] = abs(df['weight'])
        df.index = [this[0:max_label_length] for this in df.index]

        df = df.sort_values('weight', ascending=True)
        df.loc[df['sign'] == True, 'sign'] = 'r'
        df.loc[df['sign'] == False, 'sign'] = 'b'

        colors = "".join(df['sign'])

        if len(df)<50:
            df.plot(y="weight", kind=kind, color=colors, width=1, lw=1,
                title='importance plot', ax=pylab.gca(),
                fontsize=fontsize, figsize=figsize)
        else:
            df.iloc[-50:].plot(y="weight", kind=kind, color=colors[-50:], 
                width=1, lw=1,
                title='importance plot', ax=pylab.gca(),
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
        pylab.xlabel("Absolute weight", fontsize=fontsize)

        return df


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
        xscaled, Y = self._get_one_drug_data(drug_name)

        alphas, coefs1, _ = enet_path(xscaled, Y, l1_ratio=l1_ratio, 
                alphas=alphas)
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
            en = sklearn.linear_model.ElasticNet(l1_ratio=l1_ratio, 
                    alpha=alpha, max_iter=max_iter, tol=tol, 
                    selection=selection, fit_intercept=fit_intercept)
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

        try:
            self.indices1 = (self.df1 == 0).sum().sort_values().ix[0:nfeat]
            self.indices2 = (self.df2 == 0).sum().sort_values().ix[0:nfeat]
        except:
            self.indices1 = (self.df1 == 0).sum().sort().ix[0:nfeat]
            self.indices2 = (self.df2 == 0).sum().sort().ix[0:nfeat]
        names1 = self.indices1.index
        names2 = self.indices2.index
        print(names2)



    def elastic_all(self, alpha, stacked=False):
        drugids = self.drugIds
        from easydev import Progress
        pb = Progress(len(drugids))
        d = {}


        for i, drug in enumerate(drugids):
            self.elastic_net(drug, alpha=alpha)
            df = pd.DataFrame({'name': self.X.columns, 'weight': self.en.coef_})
            df = df.set_index("name").sort_values("weight")
            d[drug] = df.copy()
            pb.animate(i+1)

        
        dfall = pd.concat([d[i] for i in d.keys()], axis=1)
        dfall.columns = drugids

        if stacked is True:
            print('here')
            dfall = dfall.stack().reset_index()
            dfall.columns = ["feature", "drug", "weight"]

        # X = weights for a given feature/drug
        # Y = Sum of all weights on a given drug
        return dfall



