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
            'analysisType': 'PANCAN',
            'pval_correction_method': 'fdr',   # or qvalue
            'equal_var_ttest': True,
            'fontsize': 20
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

        # a cache to compute ANOVA
        self.individual_anova = {}


    #@do_profile()
    def anova_one_drug_one_feature(self, drug_id='Drug_1_IC50',
            feature_name='ABCB1_mut', show_boxplot=False,
            production=False):
        """Compute ANOVA and various tests on one drug and one feature

        :param bool production: if False, returns a dataframe otherwise
            a dictionary. This is to speed up analysis when scanning
            the drug across all features.
        """
        # select IC50 of a given drug
        data = self.ic50.df[drug_id]
        mask = data.isnull() == False

        # and the respective features
        self.masked_ic50 = data[mask]
        self.masked_features = self.features.df[feature_name][mask]

        # select only relevant tissues
        self.masked_tissue = self.tissue_factor[self.masked_ic50.index]
        self.masked_msi = self.msi_factor[self.masked_ic50.index]

        # Let us create an alias to indicate the main Y variable to be regressed
        self.Y = self.masked_ic50

        # about 15% of the time in those 4 lines.
        positive_feature = self.masked_features.sum()
        negative_feature = len(self.masked_features) - positive_feature
        positive_msi = self.masked_msi.sum()
        negative_msi = len(self.masked_msi) - positive_msi

        positives = self.masked_ic50[self.masked_features==1]
        negatives = self.masked_ic50[self.masked_features==0]
        Npos = len(positives)
        Nneg = len(negatives)

        # Some validity tests to run the analysis or not
        A = self.settings.includeMSI_factor and\
            positive_feature >= self.settings.featFactorPopulationTh and\
            negative_feature >= self.settings.featFactorPopulationTh and\
            negative_msi >= self.settings.MSIfactorPopulationTh and\
            positive_msi >= self.settings.MSIfactorPopulationTh
        B = (not self.settings.includeMSI_factor) and\
            positive_feature >= self.settings.featFactorPopulationTh and\
            negative_feature >= self.settings.featFactorPopulationTh

        if (A is False) and (B is False):
            drug_name = 'Drug_' + str(drug_id) + '_IC50'
            results = {'FEATURE': feature_name,
                'Drug id': drug_name,
                'Drug name': drug_name,
                'Drug Target': drug_name,
                'N_FEATURE_pos': Npos,
                'N_FEATURE_neg': Nneg,
                'log max.Conc.tested': None,
                'log max.Conc.tested2': None,
                'FEATUREpos_logIC50_MEAN': None,
                'FEATUREneg_logIC50_MEAN': None,
                'FEATURE_deltaMEAN_IC50': None,
                'FEATUREpos_IC50_sd': None,
                'FEATUREneg_IC50_sd': None,
                'FEATURE_IC50_effect_size': None,
                'FEATUREpos_Glass_delta': None,
                'FEATUREneg_Glass_delta': None,
                'FEATURE_ANOVA_pval': None,
                'Tissue_ANOVA_pval': None,
                'MSI_ANOVA_pval': None,
                'FEATURE_IC50_T_pval': None
                }
            if production is True:
                return results
            else:
                # index is not relevant here
                df = pd.DataFrame(results, index=[1])
                return df

        # First, let us create a data frame to hold the relevant data sets
        self.data = pd.DataFrame({'Y':self.Y, 'feature': self.masked_features,
            'msi':self.masked_msi, 'tissue':self.masked_tissue},
            columns=['Y', 'feature', 'msi', 'tissue'])

        # Order is important... Does not change total sum of square
        # but may change individual effects of the categorical
        # components. Not sure how important this is and will need
        # to consider other cases for robustness testing maybe.
        if self.settings.analysisType == 'PANCAN':
            # Note that tissue with less than N? values are dropped
            # This is also the case in R.
            # Possibly ntissue<3 ?
            # See e.g., Drug_1 / ABL2_mut
            #self.data_lm = ols('Y ~ C(tissue) + C(msi) + feature',
            #        data=self.data, missing='none').fit() #Specify C for Categorical

            #
            # This is faster that above but messier
            # The creation of this df represents 20% of the function time
            df = pd.get_dummies(self.data['tissue']) # could use prefix_sep
            # but there is no suffix_sep...
            df.columns = ['C(tissue)[T.'+x +']' for x in df.columns]
            df['C(msi)[T.1]'] = self.data['msi'].values
            df['feature'] = self.data['feature'].values
            df.insert(0, 'Intercept', [1]*len(df))

            # Here, we need to get rid of some of the
            #
            df = df.drop('C(tissue)[T.Bladder]', axis=1)

            self.data_lm = OLS(self.data['Y'], df).fit()

            # this sklean gives same as statsmolde.ols.params
            # and as fast as OLS but faster than 'ols'
            #lmres = LinearRegression(fit_intercept=True).fit(an.data[['msi',
            #    'feature']], an.data['Y'])
            # lmres.coef_


        elif self.settings.includeMSI_factor is True:
            self.data_lm = ols('Y ~ C(msi) + feature',
                data=self.data).fit() #Specify C for Categorical
        else:
            self.data_lm = ols('Y ~ feature',
                data=self.data).fit() #Specify C for Categorical

        # Get those stats from a local version of ANOVA
        # The only values we want is PR(>F)
        #self.stats = sm.stats.anova_lm(self.data_lm, typ=1)

        # Fvalues should be a vector with individual F values
        # df should be a vector with indivudal df
        df_tissue = len([x for x in self.data_lm.model.exog_names if 'C(tissue'
            in x])
        df = [df_tissue, 1, 1] # msi and feature have 1 df each
        endog = self.data_lm.model.data.endog
        exog = self.data_lm.model.data.exog
        q,r = np.linalg.qr(exog)
        effects = np.dot(q.T, endog)

        Nterms = 3 + 1 # msi, tissue, feature + 1 (intercept)
        Ncolumns = sum(df) + 1 # +1 for the intercept
        arr = np.zeros((Nterms, Ncolumns))
        term_names = ['Intercept', 'C(tissue)', 'C(msi)', 'feature']

        design_info = {}
        design_info['Intercept'] = (0, 1, None)
        design_info['C(tissue)'] = (1, 1+df_tissue, None)
        design_info['C(msi)'] = (df_tissue+1, df_tissue+2, None)
        design_info['feature'] = (df_tissue+2, df_tissue+3, None)

        slices = [slice(*design_info[name]) for name in term_names]
        for i,slice_ in enumerate(slices):
             arr[i, slice_] = 1
        sum_sq = np.dot(arr, effects**2)
        sum_sq = sum_sq[1:]
        mean_sq = sum_sq / np.array(df)
        Fvalues = mean_sq / (self.data_lm.ssr / self.data_lm.df_resid)
        F_pvalues = scipy.stats.f.sf(Fvalues, df, self.data_lm.df_resid)
        self.tt = F_pvalues

        #pvalues = pd.DataFrame({'C(tissue)': F_pvalues[0],
        #    'feature':F_pvalues[2], 'C(msi)': F_pvalues[1]})
        self.stats = pd.DataFrame(
                F_pvalues, columns=['PR(>F)'],
                index=['C(tissue)', 'C(msi)', 'feature'])

        # to be used with statsmodels.ols
        pvalues = self.stats['PR(>F)']
        # Then, compute t.test for p-value about feature independence
        dfeat = self.data['feature']
        # Identical to R version. Note that equal_var is True
        # is importatn. Note also that the ANOVA_results.txt
        # obtained from SFTP had different values meaning that
        # the equal.var was set to False.
        self.tfit = scipy.stats.ttest_ind(self.data['Y'][dfeat==0],
                self.data['Y'][dfeat==1],
                equal_var=self.settings.equal_var_ttest)

        # some boxplot including all data
        if show_boxplot:
            pylab.figure(1)
            neg = self.Y[self.masked_features == 0].values
            pos = self.Y[self.masked_features == 1].values
            data = {'pos':pos, 'neg':neg}
            self.data = data
            self.plot3_tuning = boxswarm.boxswarm(data)

        # some data focusing on effect of (1) tissue and (2) MSI
        if show_boxplot:
            pylab.figure(2)
            results = self._get_boxplot_data('tissue')
            if results is None:
                print("INFO: no tissue with at least 2 pos and 2 neg found. " +
                    "No image created.")
            else:
                data, names, significance = results
                bb = boxswarm.BoxSwarm(data, names)
                bb.xlabel = r'%s log(IC50)' % drug_id.replace("_", "\_")
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

            if self.settings.includeMSI_factor:
                pylab.figure(3)
                results = self._get_boxplot_data('msi')
                if results is None:
                    print("INFO: MSI with at least 2 pos and 2 neg found. " +
                        "No image created.")
                else:
                    data, names, significance = results
                    bb = boxswarm.BoxSwarm(data, names)
                    bb.xlabel = r'%s log(IC50)' % drug_id.replace("_", "\_")
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

        #iwith this index: [u'C(tissue)', u'C(msi)', u'feature', u'Residual']
        if 'C(tissue)' in pvalues.index:
            tissue_PVAL = pvalues['C(tissue)']
        else:
            tissue_PVAL = None

        if 'C(msi)' in pvalues.index:
            MSI_PVAL = pvalues['C(msi)']
        else:
            MSI_PVAL = None

        if 'feature' in pvalues.index:
            FEATURE_PVAL = pvalues['feature']
        else:
            FEATURE_PVAL = None

        # STORE value to return
        pos_IC50_mean = positives.mean()
        neg_IC50_mean = negatives.mean()
        delta_mean_IC50 = pos_IC50_mean - neg_IC50_mean

        pos_IC50_std = positives.std(ddof=1)
        neg_IC50_std = negatives.std(ddof=1)


        EFFECTSIZE_IC50 = cohens.cohens(positives, negatives)
        GLASS_d = glass.glass(positives, negatives)
        # compute cohens between IC50 where feature is pos and IC50
        # where feature is
        # negative. same for GLASS

        if drug_id.startswith("Drug"):
            drug_id = int(drug_id.split("_")[1])

        drug_name = 'Drug_' + str(drug_id) + '_IC50'
        results = {'FEATURE': feature_name,
                'Drug id': drug_name,
                'Drug name': drug_name,
                'Drug Target': drug_name,
                'N_FEATURE_pos': Npos,
                'N_FEATURE_neg': Nneg,
                'log max.Conc.tested': None,
                'log max.Conc.tested2': None,
                'FEATUREpos_logIC50_MEAN': pos_IC50_mean,
                'FEATUREneg_logIC50_MEAN': neg_IC50_mean,
                'FEATURE_deltaMEAN_IC50': delta_mean_IC50,
                'FEATUREpos_IC50_sd': pos_IC50_std,
                'FEATUREneg_IC50_sd': neg_IC50_std,
                'FEATURE_IC50_effect_size': EFFECTSIZE_IC50,
                'FEATUREpos_Glass_delta': GLASS_d[0],
                'FEATUREneg_Glass_delta': GLASS_d[1],
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

    #98% of time in  method anova_one_drug_one_feature
    def anova_one_drug(self, drug_id, animate=True):
        """Computes ANOVA for a given drug across all features

        :param str drug_id: a valid drug identifier.
        :return: a dataframe


        """
        # Takes about 10s to run. could be nice to have
        # a caching system.

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


        #df = df[df['FEATURE_ANOVA_pval'].apply(lambda x: x is not None)]

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

    def volcano_plot_all_drugs(self, df, FDR_threshold=20,
            effect_threshold=0):
        """Volcano plot for each drug

        :param df: output of :meth:`anova_all`
        """
        drugs = list(df['Drug id'].unique())
        pb = Progress(len(drugs), 1)
        for i, drug in enumerate(drugs):
            self.volcano_plot_one_drug(df, drug, FDR_threshold=FDR_threshold,
                    effect_threshold=effect_threshold)
            pylab.savefig("volcano_%s.png" % drug)
            pb.animate(i+1)

    def volcano_plot_all_features(self, df, FDR_threshold=20,
            effect_threshold=0):
        """Volcano plot for each feature

        :param df: output of :meth:`anova_all`

        Takes about 10 minutes for 265 drugs and 677 features
        """
        features = list(df['FEATURE'].unique())
        pb = Progress(len(features), 1)
        for i, feature in enumerate(features):
            self.volcano_plot_one_feature(df, feature,
                    FDR_threshold=FDR_threshold,
                    effect_threshold=effect_threshold)
            pylab.savefig("volcano_%s.png" % feature)
            pb.animate(i+1)

    def _get_volcano_global_data(self, df, FDR_threshold=20):
        varname_pval = 'FEATURE_ANOVA_pval'
        varname_qval = 'ANOVA FEATURE FDR %'

        # using all data
        minN = df['N_FEATURE_pos'].min()
        maxN = df['N_FEATURE_pos'].max()
        qvals = df[varname_qval]
        pvals = df[varname_pval]
        fdrlim = pvals[qvals<FDR_threshold].max()
        fdrlim1 = pvals[qvals<10].max()
        fdrlim2 = pvals[qvals<1].max()
        fdrlim3 = pvals[qvals<0.01].max()

        return {'minN': minN, 'maxN':maxN,
                'fdrs': {FDR_threshold:fdrlim,
                    0.01:fdrlim3, 1:fdrlim2,
                    10:fdrlim1}}

    def volcano_plot_one_feature(self,df, feature, FDR_threshold=20,
            effect_threshold=0):
        stats = self._get_volcano_global_data(df, FDR_threshold)
        data = self._get_volcano_sub_data(df, 'FEATURE', feature,
                effect_threshold=effect_threshold, FDR_threshold=FDR_threshold)

        self.volcano_plot(data.signed_effects, data.pvals, data.markersize,
                data.colors, data.annotations, stats, FDR_threshold,
                title=feature)

    def _get_volcano_sub_data(self, df, mode, target, effect_threshold=0,
            FDR_threshold=20):
        # using data related to the given drug

        if mode == 'Drug id':
            other = 'FEATURE'
        elif mode == 'FEATURE':
            other = 'Drug id'
        else:
            raise ValueError("mode parameter must be 'FEATURE' or 'Drug id'")

        subdf = df[df[mode] == target].copy()

        varname_pval = 'FEATURE_ANOVA_pval'
        varname_qval = 'ANOVA FEATURE FDR %'
        delta = subdf['FEATURE_deltaMEAN_IC50']
        effects = subdf['FEATURE_IC50_effect_size']
        signed_effects = np.sign(delta) * effects

        qvals = list(subdf[varname_qval])
        pvals = list(subdf[varname_pval])
        features = subdf[other]

        colors = []
        annotations = []
        labels = features
        if self.settings.analysisType == 'PANCAN':
            for sign, qval, pval, label in zip(signed_effects, qvals,
                    pvals, labels):
                if sign <= -effect_threshold and qval <= FDR_threshold:
                    colors.append('green')
                    annotations.append((sign,pval,label))
                elif sign >= effect_threshold and qval <= FDR_threshold:
                    colors.append('red')
                    annotations.append((sign,pval,label))
                else:
                    colors.append('grey')
        else:
            raise NotImplementedError
            #COL[which(qvals<=fdrth &
            #             pval<=gdscANOVA.settings.pval_TH & delta>0)]<-redcol
            #     COL[which(qvals<=fdrth &
            #             pval<=gdscANOVA.settings.pval_TH & delta<0)]<-greencol

        # here we normalise wrt the drug. In R code, normalised
        # my max across all data (minN, maxN)
        markersize = subdf['N_FEATURE_pos'] / subdf['N_FEATURE_pos'].max()
        markersize = list(markersize*800)
        markersize = [x if x>50 else 50 for x in markersize]

        dd = {'colors':colors, 'annotations':annotations,
                'signed_effects':signed_effects,
                'markersize':markersize, 'qvals':qvals, 'pvals':pvals}
        dd = AttrDict(**dd)
        return dd


    def volcano_plot_one_drug(self, df, drug_id, FDR_threshold=20,
            effect_threshold=0):
        """Volcano plot for one drug

        :param df: output of :meth:`anova_all`

        """
        # needs to run :meth:`anova_all` first

        # using all data, get the FDR limits
        stats = self._get_volcano_global_data(df, FDR_threshold)
        data = self._get_volcano_sub_data(df, 'Drug id', drug_id,
                effect_threshold=effect_threshold, FDR_threshold=FDR_threshold)

        self.volcano_plot(data.signed_effects, data.pvals, data.markersize,
                data.colors, data.annotations, stats, FDR_threshold,
                title=drug_id)

    def volcano_plot(self, signed_effects, pvals, markersize,
            colors, annotations, stats, FDR_threshold, title=''):

        Y = -np.log10(list(pvals)) # somehow should be cast to list ?

        pylab.clf()

        pylab.scatter(list(signed_effects), Y, s=markersize, alpha=0.4, c=colors,
                linewidth=0)

        m = abs(signed_effects.min())
        M = abs(signed_effects.max())
        l = max([m, M]) * 1.1
        pylab.xlim([-l, l])
        pylab.xlabel("Signed effect size", fontsize=self.settings.fontsize)
        pylab.ylabel('-log10(pvalues)', fontsize=self.settings.fontsize)
        pylab.ylim([0, pylab.ylim()[1]])

        #print(fdrlim, fdrlim1, fdrlim2, fdrlim3)
        fdrlim = stats['fdrs'][FDR_threshold]
        fdrlim1 = stats['fdrs'][10]
        fdrlim2 = stats['fdrs'][1]
        fdrlim3 = stats['fdrs'][0.01]
        pylab.axhline(-np.log10(fdrlim), linestyle='--',
            color='gray', alpha=1, label="FDR %s pct" % FDR_threshold)
        pylab.axhline(-np.log10(fdrlim1), linestyle='-.',
            color='gray', alpha=1, label="FDR 10 pct")
        pylab.axhline(-np.log10(fdrlim2), linestyle=':',
            color='gray', alpha=1, label="FDR 1 pct")
        pylab.axhline(-np.log10(fdrlim3), linestyle='--',
            color='black', alpha=1, label="FDR 0.01 pct")

        pylab.axvline(0, color='gray', alpha=0.5)
        ax = pylab.legend(loc='best')
        ax.set_zorder(-1) # in case there is a circle behind the legend.
        pylab.title("%s" % title.replace("_","\_"))
        for this in annotations:
            x,y,text = this
            pylab.text(x,-pylab.log10(y),text.replace("_", "\_"))


    def _get_boxplot_data(self, mode='tissue'):
        # should be called by anova_one_drug_one_feature
        # since masked_tissue, masked_ic50 attributes must
        # be populated.
        assert mode in ['tissue', 'msi']
        df = pd.DataFrame(
            {'tissue':self.masked_tissue.values,
             'ic50':self.masked_ic50,
             'feature':self.masked_features,
             'msi':self.masked_msi.values})
        if mode == 'tissue':
            df.drop('msi', inplace=True, axis=1)
        elif mode == 'msi':
            df.drop('tissue', inplace=True, axis=1)

        groups = df.groupby(['feature', mode])
        # counts items in each category and fill with NA
        counts = groups.count().unstack().fillna(0)

        cc = (counts>=2).all()
        tissues = list(cc.unstack().columns[cc])

        print(tissues)
        groups = df.query(mode + ' in @tissues').groupby([mode, 'feature'])

        # TODO; move all this if block into a method
        # figure out the delta between pos and neg
        means = groups.mean().unstack(mode)
        if len(means):
            delta = means.ix[0] - means.ix[1]
            delta.sort()
            significance = {}
            data = []
            names = []
            for tissue in delta.ix['ic50'].index:
                prefix_query = mode+"==@tissue"
                neg = df.query(prefix_query+' and feature==0')['ic50']
                pos = df.query(prefix_query+' and feature==1')['ic50']
                # HERE in the original code, equal_var is False. why ?
                res = scipy.stats.ttest_ind(neg, pos, equal_var=False)
                significance[tissue] = res[1] # p-values
                data.append(neg.values)
                data.append(pos.values)
                if mode == 'tissue':
                    name = tissue
                elif mode == 'msi':
                    if tissue == 0:
                        name = 'MSI-stable'
                    elif tissue == 1:
                        name = 'MSI-unstable'

                for this in [0.05, 0.01, 0.001]:
                    if significance[tissue] < this:
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




