import os
import pandas as pd
import scipy
import pylab
import numpy as np

import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.formula.api import OLS
from statsmodels.stats.multitest import fdrcorrection

from easydev import Progress, AttrDict

from gdsctools import boxswarm, reader
try:
    from cno.misc.profiler import do_profile
except:
    pass
from gdsctools import cohens, glass
# See reader module to get the format. The file IC50_input.txt was provided
# by Howard as a test case
#data = reader.IC50()
# Matrix with all IC50 values
#ic50 = data.ic50
# Structure containing the input features to be correlated with drug response
#features = data.features

from gdsctools.report import Report

# TODO: Could inherit from a dataframe ?
class GDSC_ANOVA_Results(object):
    def __init__(self, data=None, input_feature=None, concentrations=None):
        if data is not None and isinstance(data, str):
            self.read_csv(data)
        elif data is not None:
            self.df = data.copy() # assume it is a dataframe.

        r = reader.GenomicFeatures(input_feature)
        self.input_features = r.df
        self.varname_pval = 'FEATURE_ANOVA_pval'
        self.varname_qval = 'ANOVA FEATURE FDR %'

        if concentrations:
            self.conc = pd.read_csv(concentrations, sep='\t')
            newdata = self.conc['Drug id'].apply(lambda x: "Drug_"+str(x)+"_IC50")
            self.conc['Drug id'] = newdata
            self.conc.set_index('Drug id', inplace=True)
            df = self.df.join(self.conc, on='Drug id', how='left')
            self.df = df
            del self.conc

    def _get_pvalue_from_fdr(self, FDR_threshold=20):
        qvals = df[self.varname_qval]
        pvals = df[self.varname_pval]
        pvalue = pvals[qvals<FDR_threshold].max()
        return pvalue

    def volcano_plot(self):
        pass

    def read_csv(self, filename, sep="\t"):
        self.df = pd.read_csv(filename, sep=sep,
                comment="#")

    def _set_sensible_df(self, FDR_threshold=30, pval_threshold=None):
        if pval_threshold is None:
            pval_threshold = np.Inf

        # just an alias
        logand = np.logical_and

        # select sensible data set
        mask1 = self.df['ANOVA FEATURE FDR %'] < FDR_threshold
        mask2 = self.df['FEATURE_ANOVA_pval'] < pval_threshold
        mask3 = self.df['FEATURE_deltaMEAN_IC50'] < 0
        self.sensible_df = self.df[logand(logand(mask1, mask2), mask3)]

        # select resistant data set
        mask3 = self.df['FEATURE_deltaMEAN_IC50'] >= 0
        self.resistant_df = self.df[logand(logand(mask1, mask2), mask3)]

        # count instance per category

    def _get_data(self, df_count_sensible, df_count_resistant):

        # we can drop all columns except one, which is renamed as count
        df1 = df_count_sensible['assoc_id']
        df1.name = 'sens assoc'
        df2 = df_count_resistant['assoc_id']
        df2.name = 'res assoc'

        # Now, we join the two TimeSeries (note that above, we selected only
        # one column so the dataframe was downcast to time series)
        df_count = pd.DataFrame(df1).join(pd.DataFrame(df2), how='outer')
        # and set NA to zero
        df_count.fillna(0, inplace=True)
        # let us add a new column with the total
        df_count['total'] = df_count['sens assoc'] + df_count['res assoc']

        # we want to sort by 'total' column and is equality by the name,
        # which is the index. So let us add the index temporarily as
        # a column, sort, and remove 'name' column afterwards
        df_count['name'] = df_count.index
        df_count.sort(['total', 'name'], ascending=False, inplace=True)
        df_count.drop('name', axis=1, inplace=True)
        return df_count

    def drug_summary(self, FDR_threshold=30, pval_threshold=None,
            show=True, top=50, fontsize=10):
        # get sensible and resistant sub dataframes
        self._set_sensible_df(FDR_threshold, pval_threshold)

        # group by drug
        df_count_sensible = self.sensible_df.groupby('Drug id').count()
        df_count_resistant = self.resistant_df.groupby('Drug id').count()

        df_count = self._get_data(df_count_sensible, df_count_resistant)

        if show is True:
            self._plot(df_count, 'drug', 50)
        return df_count

    def feature_summary(self, FDR_threshold=30, pval_threshold=None,
            show=True, top=50, fontsize=10):
        # get sensible and resistant sub dataframes
        self._set_sensible_df(FDR_threshold, pval_threshold)

        df_count_sensible = self.sensible_df.groupby('FEATURE').count()
        df_count_resistant = self.resistant_df.groupby('FEATURE').count()

        df_count = self._get_data(df_count_sensible, df_count_resistant)

        # let us add another column with the n_altered_samples
        # Seems to be wrong in R code should be sum across row, i.e
        # get a vector of same length as df_count. Here it should be correct
        # but is therefore different from the R code (15oct2015)
        # Not used anyway so commented for now
        ##n_altered_samples = self.input_features[df_count.index].sum(axis=0)
        ##n_altered_samples = pd.DataFrame(n_altered_samples,
        ##    columns=['n_altered_samples'])
        ##df_count = df_count.join(pd.DataFrame(n_altered_samples), how='outer')
        ##df_count = df_count[['n_altered_samples', 'total', 'sens assoc',
        ##'res assoc']]

        #
        #s1 = set(self.sensible_df['FEATURE'])
        #s2 = set(self.resistant_df['FEATURE'])
        #size_domain = len(s1.union(s2))
        #
        #size_total_domain = len(set(self.df['FEATURE']))
        #Gperc = size_domain /  float(size_total_domain)
        if show is True:
            self._plot(df_count, 'feature', 50, fontsize=fontsize)
        return df_count

    def _plot(self, df_count, title_tag, top, fontsize=10, FDR_threshold=30):

        if top > len(df_count):
            top = len(df_count)

        df = df_count.ix[0:top][[u'sens assoc', u'res assoc']]
        labels = list(df.index)
        labels = [x.replace('_', '\_') for x in labels]
        ind = range(0, len(labels))
        ind.reverse()
        data1 = df['sens assoc'].values
        data2 = df['res assoc'].values
        pylab.clf()
        p1 = pylab.barh(ind, data1, height=0.8, color='purple',
            label='sensitivity')
        p2 = pylab.barh(ind, data2, height=0.8, color='orange',
            left=data1, label='resistance')
        ax = pylab.gca()
        self.labels = labels
        ax.set_yticks([x+0.5 for x in ind])
        ax.set_yticklabels(labels, fontsize=fontsize)
        pylab.grid()
        pylab.title(r"Top %s %s most frequently " % (top, title_tag) + \
                    "\nassociated with drug  response", fontsize=15)
        pylab.xlabel(r'Number of significant associations (FDR %s %s %s) '
                    % ("$>$", FDR_threshold, "$\%$"),  fontsize=15)
        pylab.legend(loc='lower right')
        pylab.tight_layout()

        # TODO
        print("saving figure into file to be done")
        print("Saving data into CSV to be done")

    def get_significant_hits(self, concentrations='concentrations.tsv'):
        fdrs = range(5, 50+1, 5)

        significants = []
        significant_meaningful = []
        strong_hits = []
        full_strong_hits = []

        MC1 = self.df['log max.Conc.tested']
        MC2 = self.df['log max.Conc.tested2']
        mask2 = self.df['FEATUREpos_logIC50_MEAN'] < MC1
        mask3 = self.df['FEATUREpos_logIC50_MEAN'] < MC2
        mask4 = self.df['FEATUREneg_logIC50_MEAN'] < MC1
        mask5 = self.df['FEATUREneg_logIC50_MEAN'] < MC2
        maskMC = mask2 + mask3 +mask4 +mask5

        for fdr in fdrs:
            # significant hits
            res = self.df['ANOVA FEATURE FDR %']<fdr
            significants.append(res.sum())

            # meaningful hits
            indices = np.logical_and(self.df['ANOVA FEATURE FDR %']<fdr,
                    maskMC)
            significant_meaningful.append(indices.sum())

            # meaningful strong hits
            mask1 = self.df.ix[indices]['FEATUREpos_Glass_delta'] >= 1
            mask2 = self.df.ix[indices]['FEATUREneg_Glass_delta'] >= 1
            strong_hits.append(np.logical_or(mask1, mask2).sum())

            # meaningful full strong hits
            mask1 = self.df.ix[indices]['FEATUREpos_Glass_delta'] >= 1
            mask2 = self.df.ix[indices]['FEATUREneg_Glass_delta'] >= 1
            full_strong_hits.append(np.logical_and(mask1, mask2).sum())

        data = {'significants': significants,
                'full_strong_hits': full_strong_hits,
               'strong_hits': strong_hits,
               'significant_meaningful': significant_meaningful}

        df = pd.DataFrame(data, columns = ['significants',
            'significant_meaningful', 'strong_hits', 'full_strong_hits'],
            index=fdrs)
        df.columns = ['1) significant', '2) 1 + meaningful',
        '3) 2 + strong', '4) 2+ very strong']

        pylab.clf()
        ax = pylab.gca()
        df.plot(kind='bar', width=.8, colors=['r', 'gray', 'orange', 'black'],
                rot=0, ax=ax)
        pylab.grid()
        # original is 'aquamarine4','cyan2','cornflowerblue    ','aquamarine'),
        return df

    def family_based_drug_summary(self):
        raise NotImplementedError

    def __str__(self):
        self.df.info()
        return ""


class GDSC_ANOVA(object):
    """ANOVA analysis of the IC50 vs Feature matrices

    ::

        from gdsctools import reader, anova
        r = reader.IC50('valid_file.tsv')
        an = GDSC_ANOVA(r.ic50, r.features)
        an.anova_one_drug_one_feature('Drug_1_IC50', 'TP53_mut',
            show_boxplot=True)

    """
    def __init__(self, ic50, features=None):
        """.. rubric:: Constructor

        :param DataFrame IC50: a dataframe with the IC50. Rows should be
            the COSMIC identifiers and columns should be the Drug names
            (or identifiers)
        :param features: another dataframe with rows as in the IC50 matrix
            and columns as features.  The first 3 columns must be named
            specifically to hold tissues, MSI (see format).

        The attribute :attr:`settings` contains specific settings related
        to the analysis or visulation.
        """
        # Reads IC50
        self.ic50 = reader.IC50(ic50)

        # Reads features
        if features is None:
            # Reads default version provided with the package
            self.features = reader.GenomicFeatures()
        else:
            self.features = reader.GenomicFeatures(features)

        # save the tissues
        self.tissue_factor = self.features.df['Tissue Factor Value']

        # and MSI (Microsatellite instability) status of the samples.
        self.msi_factor = self.features.df['MS-instability Factor Value']

        # settings
        self.settings = {
            # include MSI as a co-factor
            'includeMSI_factor': True,
            # number of positive samples required to perform the test
            'featFactorPopulationTh': 3,
            # How many MSI samples must be present to perform the test
            'MSIfactorPopulationTh': 2,
            'analysis_type': 'PANCAN',
            'pval_correction_method': 'fdr',   # or qvalue
            'equal_var_ttest': True,
            'fontsize': 20,
            
            }
        # makes this dict keys accessible as attributes
        self.settings = AttrDict(**self.settings)

        # is it used ?
        self.column_names = [
            'assoc_id', 'FEATURE', 'Drug id', 'Drug name',
            'Drug Target', 'N_FEATURE_pos', 'N_FEATURE_neg',
            'log max.Conc.tested', 'log max.Conc.tested2',
            'FEATUREpos_logIC50_MEAN', 'FEATUREneg_logIC50_MEAN',
            'FEATURE_deltaMEAN_IC50', 'FEATUREpos_IC50_sd',
            'FEATUREneg_IC50_sd', 'FEATURE_IC50_effect_size',
            'FEATUREpos_Glass_delta', 'FEATUREneg_Glass_delta',
            'FEATURE_ANOVA_pval', 'Tissue_ANOVA_pval',
            'MSI_ANOVA_pval', 'FEATURE_IC50_T_pval',
            'ANOVA FEATURE FDR %']

        # skip assoc_id for now
        self._odof_dict = dict([(name, None) for name in self.column_names[1:]])

        # a cache to compute ANOVA
        self.individual_anova = {}


    def _get_one_drug_one_feature_data(self, drug_name, feature_name):
        """

        return: empty dictionary if criteria not fulfilled, otherwise dictionary 
            of relevatn data
        """
        # dictionary  struture to hold results (can set values as attributes)
        dd = AttrDict()

        # select IC50 of a given drug
        # a fast way to select non-NA values from 1 column:
        dd.Y = self.ic50.df[drug_name].dropna()

        # an alias to the indices
        indices = dd.Y.index

        # select only relevant tissues/msi/features
        dd.masked_features = self.features.df[feature_name][indices]
        dd.masked_tissue = self.tissue_factor[indices]
        dd.masked_msi = self.msi_factor[indices]

        # compute length of pos/neg features and MSI 
        dd.positive_feature = dd.masked_features.values.sum()
        dd.negative_feature = len(dd.masked_features) - dd.positive_feature
        dd.positive_msi = dd.masked_msi.values.sum()
        dd.negative_msi = len(dd.masked_msi) - dd.positive_msi

        # Some validity tests to run the analysis or not
        A = self.settings.includeMSI_factor and\
            dd.positive_feature >= self.settings.featFactorPopulationTh and\
            dd.negative_feature >= self.settings.featFactorPopulationTh and\
            dd.negative_msi >= self.settings.MSIfactorPopulationTh and\
            dd.positive_msi >= self.settings.MSIfactorPopulationTh
        B = (not self.settings.includeMSI_factor) and\
            dd.positive_feature >= self.settings.featFactorPopulationTh and\
            dd.negative_feature >= self.settings.featFactorPopulationTh


        # get length final pos/neg
        # use .values to access the data: 4x fastr
        dd.positives = dd.Y.values[dd.masked_features.values==1]
        dd.negatives = dd.Y.values[dd.masked_features.values==0]
        dd.Npos = len(dd.positives)
        dd.Nneg = len(dd.negatives)

        dd.A = A
        dd.B = B
        if (A == False) and (B == False):
            dd.status = False
            return dd
        else:
            dd.status = True
        
        # compute mean and std of pos and neg sets
        dd.pos_IC50_mean = dd.positives.mean()
        dd.neg_IC50_mean = dd.negatives.mean()
        dd.delta_mean_IC50 = dd.pos_IC50_mean - dd.neg_IC50_mean
        dd.pos_IC50_std = dd.positives.std(ddof=1)
        dd.neg_IC50_std = dd.negatives.std(ddof=1)

        # Compute cohens and glass effects
        dd.effectsize_ic50 = cohens.cohens(dd.positives, dd.negatives)
        GLASS_d = glass.glass(dd.positives, dd.negatives)
        dd.pos_glass = GLASS_d[0]
        dd.neg_glass = GLASS_d[0]
        dd.feature_name = feature_name
        dd.drug_name = drug_name


        return dd

    #@do_profile()
    def anova_one_drug_one_feature(self, drug_name,
            feature_name, show_boxplot=False,
            production=False, savefig=False, directory='.'):
        """Compute ANOVA and various tests on one drug and one feature

        :param bool production: if False, returns a dataframe otherwise
            a dictionary. This is to speed up analysis when scanning
            the drug across all features.
        """

        # This extract the relevant data and some simple metrics
        odof = self._get_one_drug_one_feature_data(drug_name, feature_name)

        # if the status is False, it means the number of data points
        # in a category (e.g., positive feature) is too low.
        # If so, nothing to do, we return an 'empty' dictionary
        if odof.status is False:
            results = self._odof_dict.copy()
            results['FEATURE'] = feature_name
            results['Drug id'] = drug_name
            results['Drug name'] = drug_name
            results['Drug Target'] = drug_name
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
        # Note, however, that in statsmodels this is pretty slow because
        # it relies on an underlying code (patsy) that checks and cast
        # lots of data. The code would look like:
        #### self._mydata = pd.DataFrame({'Y':self.Y,
        ####    'tissue':self.masked_tissue,
        ####       'msi': self.masked_msi, 'feature':self.masked_features})
        #### self.data_lm = ols('Y ~ C(tissue) + C(msi) + feature',
        ####  data=self._mydata, missing='none').fit() #Specify C is category

        # Note that order is important... Does not change total sum of square
        # but may change individual effects of the categorical
        # components.

        # Yet, this is slow and we decided to use OLS function instead of
        # the recommended 'ols' api, which means we cannot use formula and
        # need to create the input data sets ourself (get_dummies here
        # below. Besides, the statsmodels.stats.anova_lm does not 
        # work with typ=1 in the version tested or gives slightly different
        # results as compared to R, so we reworte the anova_lm 
        # This looks messier but is faster than ols + anova_lm
        if self.settings.analysis_type == 'PANCAN':
            # First, split tissues into N tissue columns
            # Note that there is no suffix parameter, so we need to do it
            # ourself.:
            self._mydata = pd.DataFrame({'Y': odof.Y,
                'tissue':odof.masked_tissue,
                'msi':  odof.masked_msi, 'feature': odof.masked_features})

            df = pd.get_dummies(odof.masked_tissue)
            columns = ['C(tissue)[T.'+x +']' for x in
                    odof.masked_tissue.unique()]
            Ntissue = len(columns)
            df.columns = columns
            self.ddd = df.copy()
            # Here we set other variables with dataframe columns' names as 
            # expected by OLS
            df['C(msi)[T.1]'] = odof.masked_msi.values
            df['feature'] = odof.masked_features.values
            df.insert(0, 'Intercept', [1] * (odof.Npos + odof.Nneg))

            # Here, we need to get rid of some of the cases to agree with R....
            # ??? why ?? TODO FIXME Could be that aov in R drops
            # some columns if number of positive is not large enough ?
            # dropping bladder across all give correct answers
            # across all drugs and features ?!
            # ?? if we remove skin instead of bladder, same results...
            # ?? if we remove 2 tissues, then we starts to have different 
            # thinkgs??? 
            df = df.drop('C(tissue)[T.Bladder]', axis=1)
            Ntissue -= 1

            self.data_lm = OLS(odof.Y, df.values).fit()

            # SKLearn is also a possiblity
            # works for msi+feature but not if we include tissues ?
            # compared to R lm (not aov), we get the same
            # self.dff = df
            # ols = sklearn.linear_model.LinearRegression(fit_intercept=False)
            # ols.fit(an.dff, an.Y).coef_

        elif self.settings.includeMSI_factor is True:
            self._mydata = pd.DataFrame({'Y': odof.Y,
                'msi':  odof.masked_msi, 'feature': odof.masked_features})
            self.data_lm = ols('Y ~ C(msi) + feature',
                data=self._mydata).fit() #Specify C for Categorical
            Ntissue = 0
        else:
            self._mydata = pd.DataFrame({'Y': odof.Y,
                'feature': odof.masked_features})
            self.data_lm = ols('Y ~ feature',
                data=self._mydata).fit() #Specify C for Categorical
            Ntissue = 0

        self.anova_pvalues = self._get_anova_summary(self.data_lm, 
                Ntissue, output='dict')

        # Identical to R version. Note that equal_var is True
        # is importatn. Note also that the ANOVA_results.txt
        # obtained from SFTP had different values meaning that
        # the equal.var was set to False.
        self.tfit = scipy.stats.ttest_ind(odof.negatives, odof.positives,
                equal_var=self.settings.equal_var_ttest)

        # try/except maybe faster than if/else
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

        # some boxplot including all data
        if show_boxplot:
            self._boxplot(odof, savefig=savefig, directory=directory,
                    fignum=1)

        # a boxplot to show cell lines effects. This requires
        # the settings.analyse_type to be PANCAN
        if show_boxplot and self.settings.analysis_type == 'PANCAN':
            self._boxplot_pancan(odof, directory=directory, 
                    savefig=savefig, fignum=2, mode='tissue')

            if self.settings.includeMSI_factor:
                pylab.figure(3)
                results = self._get_boxplot_data(odof, 'msi')
                if results is None:
                    print("INFO: MSI with at least 2 pos and 2 neg found. " +
                        "No image created.")
                else:
                    data, names, significance = results
                    bb = boxswarm.BoxSwarm(data, names)
                    bb.xlabel = r'%s log(IC50)' % drug_name.replace("_", "\_")
                    bb.title = 'FEATURE/MS-instability interactions'
                    ax = bb.plot(vert=False)

                    # get info from left axis
                    common_ylim = ax.get_ylim()
                    common_ticks = ax.get_yticks()

                    self.ax = ax.twinx()
                    self.ax.set_ylim(common_ylim)
                    self.ax.set_yticks(common_ticks)
                    self.ax.set_yticklabels([len(this) for this in data])

                    pylab.tight_layout()
                    if savefig is True:
                        pylab.savefig(directory+os.sep+'ODOFmsi_{}____{}.png'.format(
                            drug_name, feature_name))
                        pylab.savefig(directory+os.sep+'ODOFmsi_{}____{}.svg'.format(
                        drug_name, feature_name))



        results = {'FEATURE': feature_name,
                'Drug id': drug_name,
                'Drug name': drug_name,
                'Drug Target': drug_name,
                'N_FEATURE_pos': odof.Npos,
                'N_FEATURE_neg': odof.Nneg,
                'log max.Conc.tested': None,
                'log max.Conc.tested2': None,
                'FEATUREpos_logIC50_MEAN': odof.pos_IC50_mean,
                'FEATUREneg_logIC50_MEAN': odof.neg_IC50_mean,
                'FEATURE_deltaMEAN_IC50': odof.delta_mean_IC50,
                'FEATUREpos_IC50_sd': odof.pos_IC50_std,
                'FEATUREneg_IC50_sd': odof.neg_IC50_std,
                'FEATURE_IC50_effect_size': odof.effectsize_ic50,
                'FEATUREpos_Glass_delta': odof.pos_glass,
                'FEATUREneg_Glass_delta': odof.neg_glass,
                'FEATURE_ANOVA_pval': FEATURE_PVAL,
                'Tissue_ANOVA_pval': tissue_PVAL,
                'MSI_ANOVA_pval': MSI_PVAL,
                'FEATURE_IC50_T_pval': self.tfit[1] # pvalues is in index 1
                }

        # 12% of the time here
        if production is True:
            return results
        else:
            df = pd.DataFrame(results, index=[1])
            return df

    def _boxplot_pancan(self, odof, mode, savefig=False, 
            directory='.',
            fignum=1, title_prefix=''):
        assert mode in ['tissue', 'msi']
        drug_name = odof.drug_name.replace("_", "\_")
        feature_name = odof.feature_name.replace("_", "\_")

        results = self._get_boxplot_data(odof, mode)
        if results is None:
            print("INFO: no tissue with at least 2 pos and 2 neg found. " +
                    "No image created.")
            return

        pylab.figure(2)
        pylab.clf()
        data, names, significance = results
        bb = boxswarm.BoxSwarm(data, names)
        bb.xlabel = r'%s log(IC50)' % drug_name
        bb.title = 'FEATURE/Cancer-type interactions'
        ax = bb.plot(vert=False)
        # get info from left axis
        common_ylim = ax.get_ylim()
        common_ticks = ax.get_yticks()

        self.ax = ax.twinx()
        self.ax.set_ylim(common_ylim)
        self.ax.set_yticks(common_ticks)
        self.ax.set_yticklabels([len(this) for this in data])
        pylab.tight_layout()
        if savefig is True:
            filename = directory + os.sep
            filename += 'ODOF_tissue_{}____{}'.format(drug_name, feature_name)
            pylab.savefig(filename + '.png')
            pylab.savefig(filename + '.svg')


    def _boxplot(self, data, savefig=False, directory='.', fignum=1):

        pylab.figure(fignum)
        pylab.clf()
        # aliases
        drug_name = data.drug_name.replace("_", "\_")
        feature_name = data.feature_name.replace("_", "\_")
        fontsize = self.settings.fontsize

        # the plot itself
        boxswarm.boxswarm({'pos': data.positives, 'neg': data.negatives},
                lw=3)
        
        pylab.title('Individual association\n {0} {1}'.format(drug_name,
            feature_name), fontsize=fontsize)
        pylab.ylabel("{0} logIC50".format(drug_name),
                fontsize=fontsize)

        pylab.tight_layout()
        if savefig is True:
            filename = directory + os.sep
            filename += 'ODOF_all_{}____{}'.format(drug_name, feature_name)
            pylab.savefig(filename + '.png')
            pylab.savefig(filename + '.svg')



    #@do_profile()
    def _get_anova_summary(self, data_lm, Ntissue, output='dict'):
        # could use this with statsmodels
        # The only values we want is PR(>F)
        # self.stats = sm.stats.anova_lm(self.data_lm, typ=1)
        # note that typ=1 does not work and typ=2,3 are differnt
        # from R version, that uses type I surely.
        q, r = np.linalg.qr(data_lm.model.data.exog)
        effects = np.dot(q.T, data_lm.model.data.endog)

        # create the W matrix using tissue and MSI if requested
        term_names = ['Intercept']
        design_info = {}
        design_info['Intercept'] = (0, 1, None)
        dof = [] 
        indices = []
        if self.settings.analysis_type == 'PANCAN':
            term_names += ['C(tissue)']
            design_info['C(tissue)'] = (1, 1+Ntissue, None)
            dof.append(Ntissue)
            indices.append('tissue')
        if self.settings.includeMSI_factor is True:
            term_names += ['C(msi)']
            design_info['C(msi)'] = (Ntissue+1, Ntissue+2, None)
            dof.append(1)
            indices.append('msi')
            shift = 1
        else:
            shift = 0
        # Feature are always used.
        term_names += ['feature']
        design_info['feature'] = (Ntissue+1+shift, Ntissue+2+shift, None)
        dof.append(1)
        indices.append('feature')

        Nterms = len(term_names)
        Ncolumns = data_lm.model.data.exog.shape[1] 
        #Ncolumns += 1 # +1 for intercept
        arr = np.zeros((Nterms, Ncolumns))

        # --  
        slices = [slice(*design_info[name]) for name in term_names]
        for i, slice_ in enumerate(slices):
             arr[i, slice_] = 1
        sum_sq = np.dot(arr, effects**2)
        sum_sq = sum_sq[1:]
        mean_sq = sum_sq / np.array(dof)
        Fvalues = mean_sq / (data_lm.ssr / data_lm.df_resid)
        F_pvalues = scipy.stats.f.sf(Fvalues, dof, data_lm.df_resid)

        sum_sq = np.append(sum_sq, data_lm.ssr)
        mean_sq = np.append(mean_sq, data_lm.mse_resid)
        F_pvalues = np.append(F_pvalues, None)
        Fvalues = np.append(Fvalues, None)
        dof.append(data_lm.model.df_resid)
        indices.append('Residuals')
        # dataframe is slow, return just the dict of pvalues by default
        if output == 'dataframe':
            anova = pd.DataFrame({'Sum Sq': sum_sq, 'Mean Sq': mean_sq, 
                'Df': dof, 'F value': Fvalues, 'PR(>F)': F_pvalues}, 
                index=indices, 
                columns=['Df', 'Sum Sq', 'Mean Sq', 'F value', 'PR(>F)'])
            return anova
        elif self.settings.analysis_type == 'PANCAN':
            return {'tissue': F_pvalues[0], 'msi':F_pvalues[1],
                    'feature':F_pvalues[2]}
        elif self.settings.includeMSI_factor is True:
            return {'msi': F_pvalues[0], 'feature':F_pvalues[1]}
        else:
            return {'feature': F_pvalues[0]}

        #return anova


    #98% of time in  method anova_one_drug_one_feature
    def anova_one_drug(self, drug_id, animate=True):
        """Computes ANOVA for a given drug across all features

        :param str drug_id: a valid drug identifier.
        :return: a dataframe


        """
        # some features can be dropped
        # TODO: parameters for settings here

        # drop first and second columns that are made of strings
        # works under python2 but not python 3. Assume that the 2 first
        #columns are the sample name and tissue feature
        # Then, we keep only cases with at least 3 features.
        # MSI could be used but is not like in original R code.
        features = self.features.df.copy()

        features = features[features.columns[3:]]
        mask = features.sum(axis=0) >= 3

        #TODO: MSI, tissues, name must always be kept
        selected_features = features[features.columns[mask]]

        # scan all features for a given drug
        assert drug_id in self.ic50.df.columns
        N = len(selected_features.columns)
        pb = Progress(N, 10)
        res = {}
        # note that we start at idnex 4 to drop sample name, tissue and MSI
        for i,feature in enumerate(selected_features.columns):
            # production True, means we do not want to create a DataFrame
            # for each call to the anova_one_drug_one_feature function
            # Instead, we require dictionaries
            this  = self.anova_one_drug_one_feature(drug_id, feature,
                    production=True)
            if this['FEATURE_ANOVA_pval'] is not None:
                res[feature] = this
            if animate is True:
                pb.animate(i+1)

        # if production is False:
        # df = pid.concat(res, ignore_index=True)
        df = pd.DataFrame.from_records(res)
        df = df.T

        # TODO: drop rows where FEATURE_ANOVA_PVAL is None
        return df

    def anova_all(self, animate=True, drugs=None, features=None):
        """Run all ANOVA tests for all drugs and all features.


        :param drugs: select a subset of drugs
        :param features: select a subset of  features (not implemented yet)

        .. todo:: features


        .. note:: comparison with version contained in this package
            gives same results. FDR (~1e-6) and FEATURE_IC50_T_pval differs
            slighlty (1e-14) especially for FDR variable with large FDR
            close to 1 but nothing to worry about.
        """
        # drop DRUG where number of IC50 (non-null) is below 5
        # axis=0 is default but we emphasize that sum is over column (i.e. drug
        vv = (self.ic50.df.isnull() == False).sum(axis=0)
        drug_names = vv.index[vv >= 6]
        self.drug_names = drug_names

        # if user provided a list of drugs, use them:
        if drugs is not None:
            # todo: check valifity of the drug names
            drug_names = drugs[:]

        N = len(drug_names)
        pb = Progress(N, 1)
        for i, drug_name in enumerate(drug_names):
            # TODO: try/except
            if drug_name in self.individual_anova.keys():
                pass
            else:
                res = self.anova_one_drug(drug_name, animate=False)
                self.individual_anova[drug_name] = res
            if animate is True:
                pb.animate(i+1)
        df = pd.concat(self.individual_anova, ignore_index=True)

        # sort all data by ANOVA p-values
        df = df.sort('FEATURE_ANOVA_pval')

        # all ANOVA have been compute individually for each drug and each
        # feature.
        # Now, we compute the FDR correction
        if self.settings.pval_correction_method == 'fdr':
            data = df['FEATURE_ANOVA_pval'].values
            FDR = fdrcorrection(data)[1]  * 100 # percentage ??
        else:
            raise NotImplementedError
            # should be qvalue correction (see qvalue library in R)
            FDR = [None] * len(df)

        # insert FDR as last column.
        df.insert(len(df.columns), 'ANOVA FEATURE FDR %', FDR)

        # insert a unique identifier as first column
        N = len(df)
        df.insert(0, 'assoc_id', range(1,N+1))
        df = df[self.column_names]

        # save as attribute
        self.anova_df = df
        return df

    def _get_boxplot_data(self, odof, mode='tissue'):
        # should be called by anova_one_drug_one_feature
        # since masked_tissue, masked_ic50 attributes must
        # be populated.
        assert mode in ['tissue', 'msi']

        # Let us use Pandas, this will be easier
        df = pd.DataFrame(
            {'tissue': odof.masked_tissue.values,
             'ic50': odof.Y,
             'feature': odof.masked_features,
             'msi': odof.masked_msi.values})

        if mode == 'tissue':
            df.drop('msi', inplace=True, axis=1)
        elif mode == 'msi':
            df.drop('tissue', inplace=True, axis=1)

        groups = df.groupby(['feature', mode])
        # counts items in each category and fill with NA
        counts = groups.count().unstack().fillna(0)

        # if positive or negative for a combo, is not>=2, drop it
        cc = (counts>=2).all()
        categories = list(cc.unstack().columns[cc])

        groups = df.query(mode + ' in @categories').groupby([mode, 'feature'])

        # TODO; move all this if block into a method
        # figure out the delta between pos and neg
        means = groups.mean().unstack(mode)
        if len(means):
            delta = means.ix[0] - means.ix[1]
            delta.sort()
            significance = {}
            data = []
            names = []
            for category in delta.ix['ic50'].index:
                prefix_query = mode+"==@category"
                neg = df.query(prefix_query+' and feature==0')['ic50']
                pos = df.query(prefix_query+' and feature==1')['ic50']
                # HERE in the original code, equal_var is False. why ?
                res = scipy.stats.ttest_ind(neg, pos, equal_var=False)
                significance[category] = res[1] # p-values
                data.append(neg.values)
                data.append(pos.values)
                if mode == 'tissue':
                    name = category
                elif mode == 'msi':
                    if category == 0:
                        name = 'MSI-stable'
                    elif category == 1:
                        name = 'MSI-unstable'

                for this in [0.05, 0.01, 0.001]:
                    if significance[category] < this:
                        name = '*' + name
                names.append(name + ' neg')
                names.append(name + ' pos')
            return (data, names, significance)
        else:
            return None


def multicore(ic50, maxcpu=4):
    """Using 4 cores, the entire analysis took 15 minutes using
    4 CPUs (16 Oct 2015).

    :param ic50: a filename or :class:`IC50` instance.

    """

    import time
    t1 = time.time()
    master = GDSC_ANOVA(ic50)

    drugs = master.ic50.drugIds

    from easydev import MultiProcessing
    t = MultiProcessing(maxcpu=maxcpu)
    # add all jobs (one per drug)
    for i, drug in enumerate(drugs):
        t.add_job(analyse_one_drug, master, drug)
    t.run()

    # populate the GDSC_ANOVA instance with the results
    for this in t.results:
        drug = this[0]
        result = this[1]
        master.individual_anova[drug] = result

    print("\nTook " + str(time.time() - t1) +  "seconds.")

    return master

def analyse_one_drug(master, drug):
    res = master.anova_one_drug(drug_id=drug, animate=False)
    return (drug, res)


class OneDrugOneFeature(Report):
    def __init__(self, ic50, features=None, drug=None, feature=None):
        self.factory = GDSC_ANOVA(ic50, features=features)
        self.drug = drug
        self.feature = feature
        self.table_class = 'dataframe-summary'

        filename = "{0}____{1}.html".format(self.drug,
                self.feature.replace(" ", "_"))

        super(OneDrugOneFeature, self).__init__(directory='odof',
                filename=filename)

    def run(self):
        df = self.factory.anova_one_drug_one_feature(self.drug,
                self.feature, savefig=True, show_boxplot=True,
                directory=self.report_directory)
        df.insert(0, 'association Id', 'a1')
        return df

    def to_html(self, df):
        for this in ['FEATURE', 'Drug id', 'association Id']:
            df[this] = df[this].apply(lambda x:
                '<a href="{0}.html">{1}</a>'.format(x,x))
        df = df.T
        return df.to_html(escape=False, classes=self.table_class,
                header=False)

    def _create_report(self, onweb=True):
        # generated pictures and results
        print('Generating data, images and HTML')
        df = self.run()

        # Create the table and add it
        html_table = self.to_html(df)
        self.add_section(html_table, 'Individual association analysis')

        section = ""
        for prefix in ['ODOFall', 'ODOFmsi', 'ODOFtissue']:
            tag = "{0}_{1}____{2}.svg".format(prefix, self.drug, self.feature)
            section += '<img src="{0}"></svg>\n'.format(tag)

        self.add_section(section, "Boxplots")


class OneDrug(Report):
    def __init__(self, ic50, drug=None):
        self.factory = GDSC_ANOVA(ic50)
        self.drug = drug
        self.table_class = 'dataframe-summary'
        filename = "{0}.html".format(self.drug)
        super(OneDrug, self).__init__(directory='odof',
                filename=filename)

    def run(self, N=20):
        from volcano_anova import VolcanoANOVA
        va = VolcanoANOVA(alldf)
        va.volcano_plot_one_drug(alldf, 'Drug_1047_IC50')
        df = None
        return df

    def to_html(self, df):
        for this in ['FEATURE', 'Drug id', 'association Id']:
            df[this] = df[this].apply(lambda x:
                '<a href="{0}.html">{1}</a>'.format(x,x))
        return df.to_html(escape=False, classes=self.table_class,
                header=False)

    def _create_report(self, onweb=True):
        # generated pictures and results
        print('Generating data, images and HTML')
        df = self.run()

        # Create the table and add it
        html_table = self.to_html(df)
        self.add_section(html_table, 'Individual association analysis')

        #section = ""
        #for prefix in ['ODOFall', 'ODOFmsi', 'ODOFtissue']:
        #    tag = "{0}_{1}____{2}.svg".format(prefix, self.drug, self.feature)
        #    section += '<img src="{0}"></svg>\n'.format(tag)

        self.add_section(section, "Boxplots")































