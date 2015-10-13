import pandas as pd
import scipy
import pylab
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import fdrcorrection
from easydev import Progress, AttrDict
from gdsc import boxswarm
from gdsc import reader
from cno.misc.profiler import do_profile
# See reader module to get the format. The file IC50_input.txt was provided
# by Howard as a test case
#data = reader.IC50()
# Matrix with all IC50 values
#ic50 = data.ic50
# Structure containing the input features to be correlated with drug response
#features = data.features


class Features(object):
    def __init__(self):
        pass

        # read features
        # selection / adding
        # additional_features ?
        # flag for additional features only
        # exclude hyper methylation
        # one gene only

class GDSCReader(object):
    pass


class GDSC_ANOVA(object):
    def __init__(self, IC50, features):
        """


        g = GDSC_ANOVA()
        g.read_ic50('ic50.tsv')
        g.read_features('features.tsv')
        g.ic50  # get df
        g.features # get df
        g.cosmicIDs # get cosmic ID used in IC50 and features
        g.featureIDs # get features identifiers for the cosmicIDs
        g.drugIds  # get drugIds used in IC50

        # For each drug, scan all featuers and search for associations.
        g.analysis()

        """
        self.ic50 = IC50.copy()
        self.features = features.copy()

        # save the tissues
        self.tissue_factor = self.features['Tissue Factor Value']

        # and MSI (Microsatellite instability) status of the samples.
        self.msi_factor = self.features['MS-instability Factor Value']

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
            'equal_var_ttest': True
            }
        # makes this dict keys accessible as attributes
        self.settings = AttrDict(**self.settings)

        self.ols_tissue_msi_feature = ols('Y ~ C(tissue) + C(msi) + feature',
                data=pd.DataFrame({'tissue':[1,2], 'Y':[1,2], 'msi':[0,1],
                    'feature':[1,2]}, 
                    columns=['Y', 'feature', 'msi', 'tissue']))
        self.init()

    #@do_profile()
    def anova_one_drug_one_feature(self, drug_id='Drug_1_IC50',
            feature_name='ABCB1_mut', show_boxplot=False):

        """


        """
        # select IC50 of a given drug
        data = self.ic50[drug_id]
        mask = data.isnull() == False

        # and the respective features
        self.masked_ic50 = data[mask]
        self.masked_features = self.features[feature_name][mask]

        # select only relevant tissues
        self.masked_tissue = self.tissue_factor[self.masked_ic50.index]
        self.masked_msi = self.msi_factor[self.masked_ic50.index]

        # Let us create an alias to indicate the main Y variable to be regressed
        self.Y = self.masked_ic50


        positive_feature = (self.masked_features == 1).sum()
        negative_feature = (self.masked_features == 0).sum()
        positive_msi = (self.masked_msi == 1).sum()
        negative_msi = (self.masked_msi == 0).sum()

        A = self.settings.includeMSI_factor and\
            positive_feature >= self.settings.featFactorPopulationTh and\
            negative_feature >= self.settings.featFactorPopulationTh and\
            positive_msi >= self.settings.MSIfactorPopulationTh and\
            positive_msi >= self.settings.MSIfactorPopulationTh
        B = (not self.settings.includeMSI_factor) and\
            positive_feature >= self.settings.featFactorPopulationTh and\
            negative_feature >= self.settings.featFactorPopulationTh

        if (A is False) and (B is False):
            results = {'FEATURE': feature_name,
                'Drug id': drug_id,
                'Drug name': 'NA',
                'Drug target': 'NA',
                'N_FEATURE_pos': Npos,
                'N_FEATURE_neg': Nneg,
                'log_max.Conc.tested': None,
                'log_max.Conc.tested2': None,
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

        # else, we do the real anova analysis.

        # First, let us create a data frame to hold the relevant data sets
        self.data = pd.DataFrame({'Y':self.Y, 'feature': self.masked_features,
            'msi':self.masked_msi, 'tissue':self.masked_tissue},
            columns=['Y', 'feature', 'msi', 'tissue'])

        # Order is important... Does not change total sum of square
        # but may change individual effects of the categorical
        # components. Not sure how important this is and will need
        # to consider other cases for robustness testing maybe.
        if self.settings.analysisType == 'PANCAN':
            self.data_lm = ols('Y ~ C(tissue) + C(msi) + feature',
                data=self.data, missing='raise').fit() #Specify C for Categorical

            #from statsmodels.formula import handle_formula_data

            #formula = self.ols_tissue_msi_feature.formula
            #tmp = handle_formula_data(data, None, formula, depth=2,
            #                                    missing='drop')
            #((endog, exog), missing_idx, design_info) = tmp

            #kwargs.update({'missing_idx': missing_idx,
            #               'missing': missing,
            #              'design_info': design_info})


            #self.ols_tissue_msi_feature.data = self.data
            #self.data_lm2 = self.ols_tissue_msi_feature.fit()
            #self.data_lm2.model.data.ynames = 'Y'

        elif self.settings.includeMSI_factor is True:
            self.data_lm = ols('Y ~ C(msi) + feature',
                data=self.data).fit() #Specify C for Categorical
        else:
            self.data_lm = ols('Y ~ feature',
                data=self.data).fit() #Specify C for Categorical

        self.stats = sm.stats.anova_lm(self.data_lm, typ=1)

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
                bb.xlabel = '%s log(IC50)' % drug_id
                bb.title = 'FEATURE/Cancer-type interactions'
                bb.plot(vert=False)
                ax = pylab.twinx()
                ax.set_yticks([i+0.5 for i in range(0, len(names))])
                ax.set_yticklabels([len(this) for this in data])
                pylab.tight_layout()
                self.plot2_tuning = bb

            if self.settings.includeMSI_factor:
                pylab.figure(3)
                results = self._get_boxplot_data('msi')
                if results is None:
                    print("INFO: MSI with at least 2 pos and 2 neg found. " +
                        "No image created.")
                else:
                    data, names, significance = results
                    bb = boxswarm.BoxSwarm(data, names)
                    bb.xlabel = '%s log(IC50)' % drug_id
                    bb.title = 'FEATURE/MS-instability interactions'
                    bb.plot(vert=False)
                    ax = pylab.twinx()
                    ax.set_yticks([i+0.5 for i in range(0, len(names))])
                    ax.set_yticklabels([len(this) for this in data])
                    pylab.tight_layout()
                    self.plot3_tuning = bb


        pvalues = self.stats['PR(>F)']

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

        FEATURE_IC50_WTT_pvalue = self.tfit[1]

        positives = self.masked_ic50[self.masked_features==1]
        negatives = self.masked_ic50[self.masked_features==0]

        pos_IC50_mean = positives.mean()
        neg_IC50_mean = negatives.mean()
        delta_mean_IC50 = pos_IC50_mean - neg_IC50_mean

        pos_IC50_std = positives.std(ddof=1)
        neg_IC50_std = negatives.std(ddof=1)

        Npos = len(positives)
        Nneg = len(negatives)

        import cohens, glass
        EFFECTSIZE_IC50 = cohens.cohens(positives, negatives)
        GLASS_d = glass.glass(positives, negatives)
        # compute cohens between IC50 where feature is pos and IC50
        # where feature is
        # negative. same for GLASS

        if drug_id.startswith("Drug"):
            drug_id = int(drug_id.split("_")[1])

        results = {'FEATURE': feature_name,
                'Drug id': drug_id,
                'Drug name': None,
                'Drug target': None,
                'N_FEATURE_pos': Npos,
                'N_FEATURE_neg': Nneg,
                'log_max.Conc.tested': None,
                'log_max.Conc.tested2': None,
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
                'FEATURE_IC50_T_pval': FEATURE_IC50_WTT_pvalue
                }

        df = pd.DataFrame(results, index=[1])
        return df

    #98% of time in  method anova_one_drug_one_feature
    def anova_one_drug(self, drug_id, animate=True):

        # some features can be dropped
        # TODO: parameters for settings here
        # FIXME : not sure if sum should be across row or columns
        mask = self.features.sum(axis=0) >= 3
        #TODO: MSI, tissues, name must always be kept
        selected_features = self.features[self.features.columns[mask]]
        print(len(selected_features.columns))

        # scan all features for a given drug
        assert drug_id in self.ic50.columns
        N = len(selected_features.columns)-3
        pb = Progress(N, 10)
        res = {}
        # note that we start at idnex 4 to drop sample name, tissue and MSI
        for i,feature in enumerate(selected_features.columns[3:]):
            res[feature] = self.anova_one_drug_one_feature(drug_id, feature)
            if animate is True:
                pb.animate(i+1)
        df = pd.concat(res, ignore_index=True)

        # TODO: drop rows where FEATURE_ANOVA_PVAL is None
        return df

    def anova_all(self, animate=True, drugs=None, features=None):
        """


        :param drugs: not used yet but maybe used to select a subset of
            features
        :param features: not used yet but maybe used to select a subset of
            features
        """
        # drop DRUG where number of IC50 (non-null) is below 5
        # axis=0 is default but we emphasize that sum is over column (i.e. drug
        vv = (self.ic50.isnull() == False).sum(axis=0)
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
        df.insert(len(df.columns), 'FDR', FDR)

        # insert a unique identifier as first column
        N = len(df)
        df.insert(0, 'assoc_id', range(1,N+1))

        # save as attribute
        self.anova_df = df
        return df

    def volcano_plot_one_drug(self, drug_id, FDR_TH=20):
        # needs to run :meth:`anova_all` first
        df = self.anova_all()
        df = df[df['Drug id'] == drug_id]

        minN = df['N_FEATURE_pos'].min()
        maxN = df['N_FEATURE_pos'].max()
        qvals = df['FDR']
        pvals = df['FEATURE_ANOVA_pval']
        delta = df['FEATURE_deltaMEAN_IC50']
        effects = df['FEATURE_IC50_effect_size']
        signed_effects = np.sign(delta) * effects
        signed_effects = signed_effects.fillna(0)

        markersize = df['N_FEATURE_pos']
        Y = -log10(pvals)
        plot(signed_effects, Y, markersize=markersize)


    def init(self):
        self.individual_anova = {}

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



    
