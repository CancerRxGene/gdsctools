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

#from gdsctools.stats import MultipleTesting
#from gdsctools import readers
from gdsctools.boxplots import BoxPlots
from gdsctools.settings import ANOVASettings
from gdsctools.anova_results import ANOVAResults


__all__ = ['ANOVA']




# Not that Logging is not used: it is not pickable and prevent
# multicore analysis.
class ANOVA(BaseModels): #Logging):
    """ANOVA analysis of the IC50 vs Feature matrices

    This class is the core of the analysis. It can be used to
    compute

    #. One association between a drug and a feature
    #. The association**S** between a drug and a set of features
    #. All assocations between a set of deugs and a set of features.

    For instance here below, we read an IC50 matrix and compute the
    association for a given drug with a specific feature.

    Note that genomic features are not provided as input but a default
    file is provided with this package that contains 677 genomic
    features for 1001 cell lines. If your IC50 contains unknown cell lines,
    you can provide your own file.

    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import IC50, ANOVA, ic50_test
        ic = IC50(ic50_test)
        an = ANOVA(ic)
        # This is to select a specific tissue
        an.set_cancer_type('breast')
        df = an.anova_one_drug_one_feature('Drug_1047_IC50',
            'TP53_mut', show=True)

    :Details about the anova analysis: In the example above, we perform a
        regression/anova test based on OLS regression. This is done for
        one feature one drug across all cell lines (tissue) in the method
        :meth:`anova_one_drug`. The regression
        takes into account the following factors: tissue, MSI and features.
        The order matters. If there is only one tissue, this factor is
        dropped. If the number of MSI values is less than a pre-defined
        parameter (see :class:`~gdsctools.settings.ANOVASettings`), it is
        dropped. The other
        methods :meth:`anova_one_drug` and :meth:`anova_all` are wrappers
        around :meth:`anova_one_drug_one_feature` to loop over all drugs, and
        loop over all drugs and all features, respectively.

    V17 : 
        gdsc.volcano_FDR_interpolation = False
        gdsc.settings.pvalue_correction_method = 'qvalue'

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
        super(ANOVA, self).__init__(ic50, genomic_features,
            drug_decode=drug_decode, verbose=verbose, low_memory=low_memory,
            set_media_factor=set_media_factor)

        self.sampling = 0
        self.pvalues_features = {}

    def _get_one_drug_one_feature_data(self, drug_name, feature_name,
            diagnostic_only=False):
        """
        return: a dictionary with relevant information. There is also
            a test to see if the data can be analysis or not. This is
            stored ad a boolean value with key called *status*.
        """
        # dictionary  struture to hold results (can set values as attributes)
        dd = AttrDict()

        # select IC50 of a given drug
        # a fast way to select non-NA values from 1 column:
        # dropna is actually faster than a method using a mask.
        #dd.Y = self.ic50.df[drug_name].dropna()
        #indices = dd.Y.index
        #dd.masked_features = self.features.df[feature_name][indices]
        #dd.masked_tissue = self.tissue_factor[indices]
        #dd.masked_msi = self.msi_factor[indices]
        #dd.positive_feature = dd.masked_features.values.sum()
        #dd.negative_feature = len(dd.masked_features) - dd.positive_feature
        #dd.positive_msi = dd.masked_msi.values.sum()
        #dd.negative_msi = len(dd.masked_msi) - dd.positive_msi
        # using a mask instead of indices is 30% slower
        #mask = self.ic50.df[drug_name].isnull()==False
        #dd.masked_features = self.features.df[feature_name][mask]
        #dd.masked_tissue = self.tissue_factor[mask]
        #dd.masked_msi = self.msi_factor[mask]

        # Amother version using a dictionary instead of dataframer is actually
        # 2-3 times faster. It requires to transform the dataframe into a
        # dictionary once for all and dropping the NA as well.
        # Now, the next line takes no time
        dd.Y = self.ic50_dict[drug_name]['Y']

        # an alias to the indices
        indices = self.ic50_dict[drug_name]['indices']
        dd.indices = indices
        # select only relevant tissues/msi/features
        if self.settings.low_memory is True:
            # This line takes 50% of the time
            dd.masked_features = self.features.df.loc[indices, feature_name]
        else:
            dd.masked_features = self.features_dict[drug_name][feature_name]

        dd.masked_tissue = self.tissue_dict[drug_name]
        if self.features.found_msi:
            dd.masked_msi = self.msi_dict[drug_name]
            dd.positive_msi = dd.masked_msi.values.sum()
            dd.negative_msi = len(dd.masked_msi) - dd.positive_msi

        if self.features.found_media:
            dd.masked_media = self.media_dict[drug_name]

        # compute length of pos/neg features and MSI
        dd.positive_feature = dd.masked_features.values.sum()
        dd.negative_feature = len(dd.masked_features) - dd.positive_feature

        # Some validity tests to run the analysis or not
        feature_threshold = self.settings.feature_factor_threshold
        msi_threshold = self.settings.MSI_factor_threshold

        A = self.settings.include_MSI_factor and\
            dd.positive_feature >= feature_threshold and\
            dd.negative_feature >= feature_threshold and\
            dd.negative_msi >= msi_threshold and \
            dd.positive_msi >= msi_threshold

        B = (not self.settings.include_MSI_factor) and\
            dd.positive_feature >= feature_threshold and\
            dd.negative_feature >= feature_threshold
        # We could of course use the mean() and std() functions from pandas or
        # numpy. We could also use the glass and cohens functions from the
        # stats module but the following code is much faster because it
        # factorises the computations of mean and variance
        dd.positives = dd.Y[dd.masked_features.values == 1]
        dd.negatives = dd.Y[dd.masked_features.values == 0]
        dd.Npos = len(dd.positives)
        dd.Nneg = len(dd.negatives)

        # additional information
        dd.feature_name = feature_name
        dd.drug_name = drug_name

        # FIXME is False does not give the same results as == False
        # in the test test_anova.py !!
        if (A == False) and (B == False):
            dd.status = False
            return dd
        else:
            dd.status = True

        if diagnostic_only is True:
            return dd

        # compute mean and std of pos and neg sets;using mean() takes 15us and
        # using the already computed sum and N takes 5us
        pos_sum = dd.positives.sum()
        neg_sum = dd.negatives.sum()
        dd.pos_IC50_mean = pos_sum / dd.Npos
        dd.neg_IC50_mean = neg_sum / dd.Nneg
        dd.delta_mean_IC50 = dd.pos_IC50_mean - dd.neg_IC50_mean

        # note the ddof to agree with R convention.
        dd.pos_IC50_std = dd.positives.std(ddof=1)
        dd.neg_IC50_std = dd.negatives.std(ddof=1)

        dd.pos_IC50_std = np.sqrt(( (dd.positives**2).sum() -
            pos_sum**2/dd.Npos)/(dd.Npos-1.))
        dd.neg_IC50_std = np.sqrt(( (dd.negatives**2).sum() -
            neg_sum**2/dd.Nneg)/(dd.Nneg-1.))

        # Compute Cohens and Glass effect size. Since underlying code
        # has lots in common, we do not use the modules but add
        # the code here below
        md = np.abs(dd.pos_IC50_mean - dd.neg_IC50_mean)
        dd.pos_glass = md / dd.pos_IC50_std
        dd.neg_glass = md / dd.neg_IC50_std

        csd = (dd.Npos - 1.) * dd.pos_IC50_std**2 + \
                (dd.Nneg - 1.) * dd.neg_IC50_std**2
        csd /= dd.Npos + dd.Nneg - 2.  # make sure this is float
        dd.effectsize_ic50 = md / np.sqrt(csd)

        # Note that equal_var is a user parameter and affects
        # results. The ANOVA_results.txt obtained from SFTP
        # have different values meaning that the equal.var param
        # was set to False. Note that pvalue is stored at index 1
        dd.ttest = self._get_ttest(dd.negatives, dd.positives)
        return dd

    def _get_ttest(self, sample1, sample2):
        # this computes the ttest.
        import scipy
        return scipy.stats.ttest_ind(sample1, sample2,
                equal_var=self.settings.equal_var_ttest)[1]

    def anova_one_drug_one_feature(self, drug_id,
            feature_name, show=False,
            production=False, directory='.'):
        """Compute ANOVA and various tests on one drug and one feature

        :param drug_id: a valid drug identifier
        :param feature_name: a valid feature name
        :param bool show: show some plots
        :param str directory: where to save the figure.
        :param bool production: if False, returns a dataframe otherwise
            a dictionary. This is to speed up analysis when scanning
            the drug across all features.

        .. note:: **for developer** this is the core of tha analysis
            and should be kept as fast as possible. 95% of the time is spent
            here.

        .. note:: **for developer** Data used in this function comes from
            _get_one_drug_one_feature_data method, which should also be kept
            as fast as possible.
        """
        if drug_id not in self.drugIds:
            raise ValueError('Unknown drug name %s. Use e.g., %s'
                    % (drug_id, self.drugIds[0]))

        if feature_name not in self.feature_names:
            # we start index at 3 to skip tissue/name/msi
            raise ValueError('Unknown feature name %s. Use e.g. one of %s'
                    % (feature_name, self.feature_names[0:3]))

        # This extract the relevant data and some simple metrics
        # This is now pretty fast accounting for 45 seconds
        # for 265 drugs and 988 features
        odof = self._get_one_drug_one_feature_data(drug_id, feature_name)
        drug_name = self.drug_decode.get_name(drug_id)
        drug_target = self.drug_decode.get_target(drug_id)

        # if the status is False, it means the number of data points
        # in a category (e.g., positive feature) is too low.
        # If so, nothing to do, we return an 'empty' dictionary
        if odof.status is False:
            results = self._odof_dict.copy()
            results['FEATURE'] = feature_name
            results['DRUG_ID'] = drug_id
            results['DRUG_NAME'] = drug_name
            results['DRUG_TARGET'] = drug_target
            results['N_FEATURE_pos'] = odof.Npos
            results['N_FEATURE_neg'] = odof.Nneg
            if production is True:
                # return a dict
                return results
            else:
                # or a dataframe; note that index is not relevant here but
                # required.
                df = pd.DataFrame(results, index=[1])
                return df

        # with the data extract, we can now compute the regression.

        # In R or statsmodels, the regression code is simple since
        # it is based on the formula notation (Y~C(msi)+feature)
        # This is also possible in statsmodels library,  however,
        # this relies on patsy, which is very slow as compared to the
        # statsmodels without formula.
        #### self._mydata = pd.DataFrame({'Y':self.Y,
        ####    'tissue':self.masked_tissue,
        ####       'msi': self.masked_msi, 'feature':self.masked_features})
        #### self.data_lm = ols('Y ~ C(tissue) + C(msi) + feature',
        ####  data=self._mydata, missing='none').fit() #Specify C is category

        # IMPORTANT: the order of the factors in the formula
        # is important. It does not change the total sum of square errors
        # but may change individual effects of the categorical
        # components.

        # Instead of using ols function, we use the OLS one so we cannot
        # use formula. Instead, we need to create manually the input
        # data. In the case of categorical data (tissue), we need to
        # create the dummy variable, which is done in the constructor
        # once for all (slow otherwise).
        if self.settings.analysis_type == 'PANCAN':
            # IMPORTANT: tissues are sorted alphabetically in R aov
            # function. Same in statsmodels but capitalised names
            # are sorted differently. In R, a<b<B<c but in Python,
            # A<B<C<a<b<c. So, 'aero' tissue is before 'Bladder' in R,
            # not in python. Since in a linear regression
            # models, the order of the factor matters and the first
            # factor is used as a reference, we decided to use same
            # convention as in R.
            # see http://statsmodels.sourceforge.net/devel/contrasts.html
            # for a good explanation

            #self._mydata = pd.DataFrame({'Y': odof.Y.copy(),
            #    'tissue':odof.masked_tissue,
            #    'msi':  odof.masked_msi, 'feature': odof.masked_features})
            #self.data_lm2 = ols('Y ~ C(tissue) + C(msi) + feature',
            #    data=self._mydata).fit() #Specify C for Categorical

            # from statsmodels.stats.anova import anova_lm
            # import statsmodels.formula.api as smf
            # df  = pd.DataFrame({'Y': odof.Y.copy(),
            #   'tissue':odof.masked_tissue,'media'
            #    odof.masked_media, 'msi':  odof.masked_msi,
            #   'feature': odof.masked_features})
            # lm = smf.ols('Y~C(tissue)+C(media)+C(msi)+feature',
            #    data=df).fit()
            #  anova_lm(lm)
            # The code above gives same answer as the code in gdsctools
            # but is slower

            # We could use pd.get_dummies but pretty slow
            # instead we create the full matrix in init() method.
            # One issue is that some columns end up with sum == 0
            # and needs to be dropped.
            df = self._tissue_dummies.ix[odof.masked_tissue.index]
            todrop = df.columns[df.values.sum(axis=0) == 0]
            if len(todrop) > 0: # use if since drop() is slow
                df = df.drop(todrop, axis=1)

            # Here we set other variables with dataframe columns' names as
            # expected by OLS.
            if self.settings.include_media_factor == False:
                todrop = [x for x in df.columns if
                        x.startswith('C(media)')]
                df = df.drop(todrop, axis=1)

            df['C(msi)[T.1]'] = odof.masked_msi.values
            df['feature'] = odof.masked_features.values

            self.Y = odof.Y
            self.EV = df.values
            # The regression and anova summary are done here
            #
            """if self.settings.regression_method == 'ElasticNet':
                self.data_lm = OLS(odof.Y, df.values).fit_regularized(
                        alpha=self.settings.regression_alpha,
                        L1_wt=self.settings.regression_L1_wt)
            elif self.settings.regression_method == 'OLS':
                self.data_lm = OLS(odof.Y, df.values).fit()
            elif self.settings.regression_method == 'Ridge':
                self.data_lm = OLS(odof.Y, df.values).fit_regularized(
                        alpha=self.settings.regression_alpha,
                        L1_wt=0)
            elif self.settings.regression_method == 'Lasso':
                self.data_lm = OLS(odof.Y, df.values).fit_regularized(
                        alpha=self.settings.regression_alpha,
                        L1_wt=1)
            """
            # example of computing null model ?
            # Example of computing pvalues ourself
            # with 100 000 samples, we can get a smooth distribution
            # that we can then fit with fitter. good distribution 
            # for the raw data is uniform one but if we take the log10, 
            # we have lots of possible distrob such as beta, exponweib, gamma,
            #....
        elif self.settings.include_MSI_factor is True:
            #self._mydata = pd.DataFrame({'Y': odof.Y,
            #    'msi':  odof.masked_msi, 'feature': odof.masked_features})
            #self.data_lm = ols('Y ~ C(msi) + feature',
            #    data=self._mydata).fit() #Specify C for Categorical
            df = pd.DataFrame()
            df['C(msi)[T.1]'] = odof.masked_msi.values
            df['feature'] = odof.masked_features.values
            df.insert(0, 'Intercept', [1] * (odof.Npos + odof.Nneg))
            #self.data_lm = OLS(odof.Y, df.values).fit()
        else:
            df = pd.DataFrame()
            df['feature'] = odof.masked_features.values
            df.insert(0, 'Intercept', [1] * (odof.Npos + odof.Nneg))
            #self.data_lm = OLS(odof.Y, df.values).fit()
            #self._mydata = pd.DataFrame({'Y': odof.Y,
            #    'feature': odof.masked_features})
            #self.data_lm = ols('Y ~ feature',
            #    data=self._mydata).fit() #Specify C for Categorical

        if self.settings.regression_method == 'ElasticNet':
            self.data_lm = OLS(odof.Y, df.values).fit_regularized(
                    alpha=self.settings.regression_alpha,
                    L1_wt=self.settings.regression_L1_wt)
        elif self.settings.regression_method == 'OLS':
            self.data_lm = OLS(odof.Y, df.values).fit()
        elif self.settings.regression_method == 'Ridge':
            self.data_lm = OLS(odof.Y, df.values).fit_regularized(
                    alpha=self.settings.regression_alpha,
                    L1_wt=0)
        elif self.settings.regression_method == 'Lasso':
            self.data_lm = OLS(odof.Y, df.values).fit_regularized(
                    alpha=self.settings.regression_alpha,
                    L1_wt=1)


        key = drug_id + "__" + feature_name
        if self.sampling and key not in self.pvalues_features.keys():
            # This can be computed for a drug once for all
            # no need to redo it for each feature ?
            # If the length of Y is too small (e.g., < 20) the results may not be
            # great. This can be check zith the errors
            self.samples1 = []
            self.samples2 = []
            self.samples3 = []
            Y = odof.Y.copy()
            N = self.sampling
            pb = Progress(N, 20)
            for i in range(0, N):

                # To get the random distribution, shuffle Y
                # and noise not required
                # To get the noise effects, do not shuffle and set noise to
                # something different from 0
                noise = 0.0
                pylab.shuffle(Y)
                #data_lm = OLS(Y, df.values).fit()
                data_lm = OLS(Y+noise*pylab.randn(len(Y)), df.values).fit()
                anova_pvalues = self._get_anova_summary(data_lm,
                    output='dict')
                try:self.samples1.append(anova_pvalues['msi'])
                except:pass
                self.samples2.append(anova_pvalues['feature'])
                try:self.samples3.append(anova_pvalues['tissue'])
                except:pass
                #pb.animate(i+1)
            import fitter
            ff = fitter.Fitter(-pylab.log10(self.samples2))
            dist = "genexpon"
            ff.distributions = [dist]
            ff.fit()
            self.pvalues_features[key] = {
                'error': ff.df_errors.ix[dist].values[0],
                'params': ff.fitted_param[dist],
                'feature': feature_name,
                'N':len(Y)
            }
            print(self.pvalues_features[key])


        self.anova_pvalues = self._get_anova_summary(self.data_lm,
                 output='dict')

        # Store the pvalues. Note that some may be missing so we use try
        # except, which is faster than if/else
        try:
            tissue_PVAL = self.anova_pvalues['tissue']
        except:
            tissue_PVAL = None

        try:
            MSI_PVAL = self.anova_pvalues['msi']
        except:
            MSI_PVAL = None

        try:
            FEATURE_PVAL = self.anova_pvalues['feature']
        except:
            FEATURE_PVAL = None

        try:
            MEDIA_PVAL = self.anova_pvalues['media']
        except:
            MEDIA_PVAL = None

        if show is True:
            boxplot = BoxPlots(odof, savefig=self.settings.savefig,
                    directory=directory)
            boxplot.boxplot_association(fignum=1)

            # a boxplot to show cell lines effects. This requires
            # the settings.analyse_type to be PANCAN
            if self.settings.analysis_type == 'PANCAN':
                boxplot.boxplot_pancan(fignum=2, mode='tissue')
            if self.settings.include_MSI_factor:
                boxplot.boxplot_pancan(fignum=3, mode='msi')

        results = {'FEATURE': feature_name,
                'DRUG_ID': drug_id,
                'DRUG_NAME': drug_name,
                'DRUG_TARGET': drug_target,
                'N_FEATURE_pos': odof.Npos,
                'N_FEATURE_neg': odof.Nneg,
                'FEATURE_pos_logIC50_MEAN': odof.pos_IC50_mean,
                'FEATURE_neg_logIC50_MEAN': odof.neg_IC50_mean,
                'FEATURE_delta_MEAN_IC50': odof.delta_mean_IC50,
                'FEATURE_pos_IC50_sd': odof.pos_IC50_std,
                'FEATURE_neg_IC50_sd': odof.neg_IC50_std,
                'FEATURE_IC50_effect_size': odof.effectsize_ic50,
                'FEATURE_pos_Glass_delta': odof.pos_glass,
                'FEATURE_neg_Glass_delta': odof.neg_glass,
                'ANOVA_FEATURE_pval': FEATURE_PVAL,
                'ANOVA_TISSUE_pval': tissue_PVAL,
                'ANOVA_MSI_pval': MSI_PVAL,
                'ANOVA_MEDIA_pval': MEDIA_PVAL,
                'FEATURE_IC50_T_pval': odof.ttest # pvalues is in index 1
                }

        # 12% of the time here
        if production is True:
            return results
        else:
            df = pd.DataFrame(results, index=[1])
            return df

    def optimise_elastic_net(self, drug_name, feature_name, N=20, Nalpha=20):
        lwts = pylab.linspace(0, 1, N)
        alphas = pylab.linspace(0, 5, Nalpha)

        mses = np.zeros((N, Nalpha))

        pb = Progress(N)
        for i, lwt in enumerate(lwts):
            for j, alpha in enumerate(alphas):
                self.settings.regression_method = 'ElasticNet'
                self.settings.regression_alpha = alpha
                self.settings.regression_L1_wt = lwt
                odof = self.anova_one_drug_one_feature(drug_name,
                        feature_name)
                anova = self._get_anova_summary(self.data_lm,
                        output='dataframe')
                mses[i,j] = self.data_lm.bic
            pb.animate(i+1)
        return mses

    def optimise_ridge(self, drug_name, feature_name, alphas=None):
        return self._opt_ridge_lasso(drug_name, feature_name,
                'Ridge', alphas=alphas)

    def optimise_lasso(self, drug_name, feature_name, alphas=None):
        return self._opt_ridge_lasso(drug_name, feature_name,
                'Lasso', alphas=alphas)

    def _opt_ridge_lasso(self, drug_name, feature_name, method, alphas=None):

        if alphas is None:
            alphas = pylab.linspace(0,1, 20)

        mses = []
        params = []
        method_buf = self.settings.regression_method
        alpha_buf = self.settings.elastic_net.alpha

        pb = Progress(len(alphas))
        for j, alpha in enumerate(alphas):
            self.settings.regression_method = method
            self.settings.elastic_net.alpha = alpha
            odof = self.anova_one_drug_one_feature(drug_name,
                    feature_name)
            anova = self._get_anova_summary(self.data_lm,
                    output='dataframe')
            #mses.append(anova.ix['Residuals']['Sum Sq'])
            mses.append(anova.ix['tissue']['F value'])
            #mses.append(anova['Sum Sq'].sum())
            pb.animate(j+1)
            params.append(self.data_lm.params)
        self.settings.regression_method = method_buf
        self.settings.elastic_net.alpha = alpha_buf
        return alphas, mses, params

    # no need to optimise anymore
    def _get_anova_summary(self, data_lm, output='dict'):
        # could use this with statsmodels but somehow anova_lm with typ I
        # does not work, which is the one used in R version, so we implement
        # the anova here
        q, r = np.linalg.qr(data_lm.model.data.exog)
        effects = np.dot(q.T, data_lm.model.data.endog)

        # In the regression, the first tissue is dropped hence -1
        # The degree of freedom for tissues is N - 1
        # self.features.tissues contains all tissues even those that
        # were dropped due to lack of pos or neg features. So, we must use
        modes = self._get_analysis_mode()
        Ncolumns = data_lm.model.data.exog.shape[1]
        Ntissue = Ncolumns
        Ntissue -= 1 # remove intercept
        Ntissue -= 1 # remove feature, which are always taken into account
        if 'msi' in modes:
            Ntissue -= 1
        if 'media' in modes:
            Nmedia = len(self._media_dummies.columns)-1
        else:
            Nmedia = 0
        Ntissue -= Nmedia

        # create the W matrix using tissue and MSI if requested
        # default is that the 3 features are used
        if 'tissue' in modes and 'msi' in modes and 'media' in modes:
            dof = [Ntissue, Nmedia, 1, 1]
            indices = ['tissue', 'media', 'msi', 'feature', 'Residuals']
            # 4 stands for intercept + tissue + msi +feature
            arr = np.zeros((5, Ncolumns))
            arr[1, slice(1, Ntissue+1)] = 1
            arr[2, slice(Ntissue+1, Ntissue+Nmedia+1)] = 1
            arr[3, Ntissue + Nmedia + 1] = 1
            arr[4, Ntissue + Nmedia + 2] = 1
            self.arr = arr
        elif 'tissue' in modes and 'msi' in modes:
            dof = [Ntissue, 1, 1]
            indices = ['tissue', 'msi', 'feature', 'Residuals']
            # 4 stands for intercept + tissue + msi +feature
            arr = np.zeros((4, Ncolumns))
            arr[1, slice(1, Ntissue)] = 1
            arr[2, Ntissue + 1] = 1
            arr[3, Ntissue + 2] = 1
            self.arr = arr
        elif 'tissue' not in modes and 'msi' in modes:
            dof = [1, 1]
            indices = ['msi', 'feature', 'Residuals']
            # 3 stands for intercept + msi +feature
            arr = np.zeros((3, Ncolumns))
            arr[1, 1] = 1
            arr[2, 2] = 1
        elif 'tissue' not in modes and 'msi' not in modes:
            dof = [1]
            indices = ['feature', 'Residuals']
            # 3 stands for intercept + msi +feature
            arr = np.zeros((2, Ncolumns))
            arr[1, 1] = 1
        arr[0, 0] = 1                   # intercept

        sum_sq = np.dot(arr, effects**2)[1:] # drop the intercep
        mean_sq = sum_sq / np.array(dof)
        Fvalues = mean_sq / (data_lm.ssr / data_lm.df_resid)
        F_pvalues = scipy.stats.f.sf(Fvalues, dof, data_lm.df_resid)
        sum_sq = np.append(sum_sq, data_lm.ssr)
        mean_sq = np.append(mean_sq, data_lm.mse_resid)
        F_pvalues = np.append(F_pvalues, None)
        Fvalues = np.append(Fvalues, None)
        dof.append(data_lm.model.df_resid)
        #indices.append('Residuals')
        # dataframe is slow, return just the dict of pvalues by default
        if output == 'dataframe':
            anova = pd.DataFrame({'Sum Sq': sum_sq, 'Mean Sq': mean_sq,
                'Df': dof, 'F value': Fvalues, 'PR(>F)': F_pvalues},
                index=indices,
                columns=['Df', 'Sum Sq', 'Mean Sq', 'F value', 'PR(>F)']
                )
            return anova
        elif self.settings.analysis_type == 'PANCAN':
            if self.settings.include_media_factor:
                dd = {'tissue': F_pvalues[0],
                      'media': F_pvalues[1],
                      'msi': F_pvalues[2],
                      'feature': F_pvalues[3]}
            else:
                dd = {'tissue': F_pvalues[0],
                      'msi':F_pvalues[1],
                      'feature':F_pvalues[2]}
            return dd
        elif self.settings.include_MSI_factor is True:
            return {'msi': F_pvalues[0], 'feature': F_pvalues[1]}
        else:
            return {'feature': F_pvalues[0]}

    def _draft(self):
        # using sklearn
        #ols = linear_model.LinearRegression()
        #f = ols.fit(an.dff, an.Y)
        #sse = sum(np.square((f.predict(an.dff).T - an.Y))) /
        #           float(an.dff.shape[0] - an.dff.shape[1])
        # ssr = sum(np.square((f.predict(an.dff).T - an.Y.mean())))
        pass

    def _test(self):
        # for drug1047 and featuer ABCB1_mut
        print("""Analysis of Variance Table

        Response: Y
        Df  Sum Sq Mean Sq F value  Pr(>F)
        TISSUEpattern  26  352.35 13.5517  9.2685 < 2e-16 ***
        MSIpattern      1    5.31  5.3094  3.6313 0.05705 .
        FEATpattern     1    3.19  3.1861  2.1791 0.14028
        Residuals     817 1194.55  1.4621
        """)

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
        # 
        for i, feature in enumerate(selected_features.columns):
            # production True, means we do not want to create a DataFrame
            # for each call to the anova_one_drug_one_feature function
            # Instead, we require dictionaries
            this = self.anova_one_drug_one_feature(drug_id, feature,
                    production=True)
            if this['ANOVA_FEATURE_pval'] is not None:
                res[feature] = this
            if animate is True:
                pb.animate(i+1)

        # if production is False:
        # df = pid.concat(res, ignore_index=True)
        df = pd.DataFrame.from_records(res)
        df = df.T

        df = ANOVAResults().astype(df)
        if len(df) == 0:
            return df

        # append DRUG_NAME/DRUG_TARGET columns
        df = self.drug_decode.drug_annotations(df)

        # TODO: drop rows where ANOVA_FEATURE_PVAL is None
        if output != 'object':
            df = self.add_pvalues_correction(df)
            return df
        else:
            df = self.add_pvalues_correction(df)
            res = ANOVAResults(df)
            res.settings = ANOVASettings(**self.settings)
            return res

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
        else:
            pass

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



"""

Script to compute the pvalues based on esti;ation of the pvalues distribution
for a given drug and feature


an = anova.ANOVA("IC50_v18.csv", "GF_BLCA.csv", "DRUG_DECODE.csv")
an.sampling = 1000

def get_data(an):
    keys = an.pvalues_features.keys()
    a = []
    b = []
    c = []
    for key in keys:
        drug, feature = key.split("__")
        params = an.pvalues_features[key]['params'
        ]
        data =res.df.query("DRUG_ID==@drug and FEATURE==@feature"); pval = data["ANOVA_FEATURE_pval"]
        pval_corr =  scipy.stats.genexpon.sf(-log10(pval), *params)
        signed_effect = np.sign(data['FEATURE_delta_MEAN_IC50']) * data['FEATURE_IC50_effect_size']
        a.append(-log10(pval_corr[0]))
        b.append(signed_effect.values[0])
    return b, a


plot(get_data(an)[0], get_data(an)[1], 'o',  markersize=20, alpha=0.3, color='g')

signed_effects = np.sign(res.df['FEATURE_delta_MEAN_IC50']).values *
                        res.df['FEATURE_IC50_effect_size'].values
pvalues = -log10(res.df["ANOVA_FEATURE_pval"])
plot(signed_effects, pvalues, 'o',  markersize=20, alpha=0.3, color='r')






"""

