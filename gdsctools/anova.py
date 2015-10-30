# -*- python -*-
# -*- coding utf-8 -*-

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
and genomic features



"""
import os
import pandas as pd
import scipy
import pylab
import numpy as np

#import statsmodels.api as sm
#from statsmodels.formula.api import ols
from statsmodels.formula.api import OLS
from statsmodels.stats.multitest import fdrcorrection

from easydev import Progress, AttrDict, Logging
import easydev
from gdsctools import boxswarm, readers
try:
    from cno.misc.profiler import do_profile
except:
    pass
#from gdsctools import cohens, glass
# See reader module to get the format. The file IC50_input.txt was provided
# by Howard as a test case
#data = reader.IC50()
# Matrix with all IC50 values
#ic50 = data.ic50
# Structure containing the input features to be correlated with drug response
#features = data.features

from gdsctools.report import Report, HTMLTable
from gdsctools.tools import Savefig
from colormap import cmap_builder
from gdsctools.volcano import VolcanoANOVA
from gdsctools.settings import ANOVASettings


__all__ = [ 'ANOVA', 'ANOVAReport']


class ColumnTypes(object):

    def __init__(self):

        self.mapping = {
             'Drug Target': np.dtype('O'),
             'Drug id': np.dtype('O'),
             'Drug name': np.dtype('O'),
             'FEATURE': np.dtype('O'),
             'FEATURE_ANOVA_pval': np.dtype('float64'),
             'FEATURE_IC50_T_pval': np.dtype('float64'),
             'FEATURE_IC50_effect_size': np.dtype('float64'),
             'FEATURE_deltaMEAN_IC50': np.dtype('float64'),
             'FEATUREneg_Glass_delta': np.dtype('float64'),
             'FEATUREneg_IC50_sd': np.dtype('float64'),
             'FEATUREneg_logIC50_MEAN':  np.dtype('float64'),
             'FEATUREpos_Glass_delta': np.dtype('float64'),
             'FEATUREpos_IC50_sd': np.dtype('float64'),
             'FEATUREpos_logIC50_MEAN': np.dtype('float64'),
             'MSI_ANOVA_pval': np.dtype('O'),
             'N_FEATURE_neg': np.dtype('int64'),
             'N_FEATURE_pos': np.dtype('int64'),
             'Tissue_ANOVA_pval': np.dtype('O'),
             'log max.Conc.tested': np.dtype('O'),
             'log max.Conc.tested2': np.dtype('O')}

    def astype(self, df):
        try:
            # does not work in python3.3 on travis but should work
            # we newer pandas version.
            df = df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
        except:
            pass
        return df



class ANOVAReport(Savefig):
    """Class used to create final HTML reports on an :class:`ANOVA` analysis


    ::

        from gdsctools import *

        # Perform the analysis itself to get a set of results (dataframe)
        an = ANOVA(ic50_test)
        results = an.anova_all()

        # now, we can create the report. Note that we must currently provide
        # a file that contains the concentrations used for each drug
        r = ANOVAReport(gdsc=an, results=results, 
            concentrations='concentrations.csv')
        r.settings.pvalue_threshold = 0.001
        r.settings.FDR_threshold = 28
        r.settings.directory = 'testing'
        r.report()

    """
    def __init__(self, gdsc, results, concentrations=None, sep="\t"):
        """.. rubric:: Constructor


        :param gdsc: the instance with which you created the results to report
        :param results: the results returned by :meth:`ANOVA.anova_all`
        :param concentrations: todo

        """
        super(ANOVAReport, self).__init__()

        data = results
        # data can be a file with all results as exported
        # by ANOVA analysis
        if data is not None and isinstance(data, str):
            print("Reading the results from a file")
            self.df = self.read_csv(data, sep=sep)
        elif data is not None:
            try:
                self.df = data.df.copy()
            except:
                # or an instance of a dataframe
                self.df = data.copy()

        self.settings = ANOVASettings()
        for k, v in gdsc.settings.items():
            self.settings[k] = v

        self._colname_drug_id = 'Drug id'
        self.varname_pval = 'FEATURE_ANOVA_pval'
        self.varname_qval = 'ANOVA FEATURE FDR %'

        self.ic50 = readers.IC50(gdsc.ic50)
        self.input_features = readers.GenomicFeatures(gdsc.features).df

        if concentrations:
            # input may not have the concentrations columns right now.
            # This should be fixed in input data set
            if "log max.Conc.tested" in self.df.columns:
                print("your dataframe already contains concentration. replace them")
            self.df.drop('log max.Conc.tested', axis=1, inplace=True)
            self.df.drop('log max.Conc.tested2', axis=1, inplace=True)
            drugid = self._colname_drug_id
            self.conc = pd.read_csv(concentrations, sep='\t')
            newdata = self.conc[drugid].apply(lambda x: "Drug_"+str(x)+"_IC50")
            self.conc[drugid] = newdata
            self.conc.set_index(drugid, inplace=True)
            df = self.df.join(self.conc, on=drugid, how='left')
            self.df = df
            del self.conc

        # create some data
        self._set_sensible_df()

    def _get_ndrugs(self):
        return len(self.df[self._colname_drug_id].unique())
    n_drugs = property(_get_ndrugs, doc="return number of drugs")

    def _get_nfeatures(self):
        # !! -3 to remove sample name, tissue, msi columns
        return len(self.input_features.columns) - 3
    n_features = property(_get_nfeatures, 
            doc="return number of features ignoring MSI, sample and tissue")

    def _get_ntests(self):
        return len(self.df.index)
    n_tests = property(_get_ntests)

    def _get_ncelllines(self):
        return len(self.input_features.index)
    n_celllines = property(_get_ncelllines, 
            doc="return number of cell lines")

    def _df_append(self, df, data):
        count = len(df)
        df.ix[count] = data
        return df

    def diagnostics(self):
        """Return summary of the analysis (dataframe)"""
        df = pd.DataFrame({'text':[], 'value':[]})
        txt = []

        N = float(self.n_drugs*self.n_features)
        ratio = float(self.n_tests)/(N) * 100
        ratio = easydev.precision(ratio, digit=2)

        msg = "Type of analysis"
        df = self._df_append(df, [msg, self.settings.analysis_type])

        msg = "Total number of possible drug/feature associations"
        df = self._df_append(df, [msg, N])
        msg = "Total number of ANOVA tests performed"
        df = self._df_append(df, [msg, self.n_tests])
        msg = "Percentage of tests performed"
        df = self._df_append(df, [msg, ratio])
        
        # trick to have an empty line
        df = self._df_append(df, ["", ""])

        msg = "Total number of tested drugs"
        df = self._df_append(df, [msg, self.n_drugs])
        msg = "Total number of genomic features"
        df = self._df_append(df, [msg, self.n_features])

        msg = "Total number of screened cell lines"
        df = self._df_append(df, [msg, self.n_celllines])

        msg = "MicroSatellite instability included as factor"
        msi = self.settings.includeMSI_factor
        df = self._df_append(df, [msg, msi])

        # trick to have an empty line
        df = self._df_append(df, ["", ""])
        nsens = len(self.sensible_df)
        nres = len(self.resistant_df)
        msg = "Total number of significant associations"
        df = self._df_append(df, [msg, nsens+nres])
        msg = " - sensitive"
        df = self._df_append(df, [msg, nsens])
        msg = " - resistant"
        df = self._df_append(df, [msg, nres])

        msg = "p-value significance threshold"
        df = self._df_append(df, [msg, self.settings.FDR_threshold])
        
        p1, p2 = self._get_pval_range()
        msg = 'Range of significant p-values'
        value = "[{:.4}, {:.4}]".format(p1,p2)
        df = self._df_append(df, [msg, value])

        f1, f2 = self._get_fdr_range()
        f1 = easydev.precision(f1, 3)
        f2 = easydev.precision(f2, 3)
        msg = "range of significant % FDRs"
        value = '[{} {}]'.format(f1,f2)
        df = self._df_append(df, [msg, value])
        #return "\n".join(txt)
        return df

    def _get_pval_range(self):
        """Get pvalues range of the significant hits"""
        nsens = len(self.sensible_df)
        nres = len(self.resistant_df)
        N = nsens + nres
        if N == 0:
            return 0,0
        name = self.varname_pval
        data = self.df[name].ix[0:N-1]
        m, M = data.min(), data.max()
        return m,M

    def _get_fdr_range(self):
        """Get FDR range of the significant hits"""
        name = self.varname_qval
        data = self.df[name][(self.df[name]< self.settings.FDR_threshold)]
        if len(data) == 0:
            return 0,0
        m, M = data.min(), data.max()
        return m,M

    def _get_pvalue_from_fdr(self):
        qvals = df[self.varname_qval]
        pvals = df[self.varname_pval]
        pvalue = pvals[qvals < self.settings.FDR_threshold].max()
        return pvalue

    def read_csv(self, filename, sep="\t"):
        """Reads a result file and store it in :attr:`df`"""
        self.df = pd.read_csv(filename, sep=sep,
                comment="#")
        return self.df

    def _set_sensible_df(self):
        # just an alias
        logand = np.logical_and

        # select sensible data set
        mask1 = self.df['ANOVA FEATURE FDR %'] < self.settings.FDR_threshold
        mask2 = self.df['FEATURE_ANOVA_pval'] < self.settings.pvalue_threshold
        mask3 = self.df['FEATURE_deltaMEAN_IC50'] < 0
        self.sensible_df = self.df[logand(logand(mask1, mask2), mask3)]

        # select resistant data set
        mask3 = self.df['FEATURE_deltaMEAN_IC50'] >= 0
        self.resistant_df = self.df[logand(logand(mask1, mask2), mask3)]

    def get_significant_set(self):
        """Return signiificant hits (resistant and sensible)"""
        # a property that is long to compute
        # and may change if FDR changes.
        self._set_sensible_df()
        df = pd.concat([self.sensible_df, self.resistant_df])
        try:
            df.sort_values('assoc_id', inplace=True)
        except:
            df.sort('assoc_id', inplace=True)
        return df

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
        try:
            df_count.sort_values(by=['total', 'name'], ascending=False,
                inplace=True)
        except:
            df_count.sort(columns=['total', 'name'], ascending=False,
                inplace=True)

        df_count.drop('name', axis=1, inplace=True)
        return df_count

    def drug_summary(self,  top=50, fontsize=10):
        """Return dataframe with significant drugs"""
        # get sensible and resistant sub dataframes
        self._set_sensible_df()

        # group by drug
        df_count_sensible = self.sensible_df.groupby('Drug id').count()
        df_count_resistant = self.resistant_df.groupby('Drug id').count()

        df_count = self._get_data(df_count_sensible, df_count_resistant)

        if len(df_count):
            self._plot(df_count, 'drug', top)
            if self.settings.savefig is True:
                self.savefig('drug_summary.png')

        return df_count

    def feature_summary(self, top=50, fontsize=10):
        """Return dataframe with significant features"""
        # get sensible and resistant sub dataframes
        self._set_sensible_df()

        df_count_sensible = self.sensible_df.groupby('FEATURE').count()
        df_count_resistant = self.resistant_df.groupby('FEATURE').count()

        df_count = self._get_data(df_count_sensible, df_count_resistant)

        if len(df_count)>0:
            self._plot(df_count, 'feature', top, fontsize=fontsize)
            if self.settings.savefig is True:
                self.savefig('feature_summary.png')
        return df_count

    def _plot(self, df_count, title_tag, top, fontsize=10):
        """Used by drug_summary and feature_summary to plot the
        bar plot"""
        if top > len(df_count):
            top = len(df_count)

        df = df_count.ix[0:top][[u'sens assoc', u'res assoc']]
        labels = list(df.index)
        labels = [x.replace('_', '\_') for x in labels]
        ind = range(0, len(labels))
        # reverse does not exist with python3
        try:
            ind.reverse()
        except:
            ind = list(ind)
            ind.reverse()
        data1 = df['sens assoc'].values
        data2 = df['res assoc'].values
        pylab.figure(1)
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
                    % ("$>$", self.settings.FDR_threshold, "$\%$"),
                    fontsize=15)
        pylab.legend(loc='lower right')
        pylab.tight_layout()

    def get_significant_hits(self, concentrations='concentrations.tsv', 
            show=True):
        """Return a summary of significant hits
        
        :param show: show a plot with the distribution of significant hits

        """
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

        if show is True:
            pylab.clf()
            ax = pylab.gca()
            df.plot(kind='bar', width=.8, 
                    color=['r', 'gray', 'orange', 'black'],
                    rot=0, ax=ax)
            pylab.grid()
        # original is 'aquamarine4','cyan2','cornflowerblue    ','aquamarine'),
        return df

    def family_based_drug_summary(self):
        raise NotImplementedError

    def __str__(self):
        self.df.info()
        return ""

    def check(self):
        """for developers (check local file agreement). will be removed."""
        rold = ANOVAReport("anova_all.tsv",
                concentrations='concentrations.tsv')
        rold.df = rold.df[self.df.columns]
        for x in self.df.columns:
            try:
                print(x, max(self.df[x] - rold.df[x]))
            except:
                print(x, all(self.df[x] == rold.df[x]))

    def create_html_associations(self):
        """Create an HTML page for each significant association"""
        print("\n\nCreating individual HTML pages for each association")
        df = self.get_significant_set()
        drugs = df['Drug id'].values
        features = df['FEATURE'].values
        assocs = df['assoc_id'].values
        fdrs = df['ANOVA FEATURE FDR %'].values
        N = len(drugs)
        pb = Progress(N)
        html = OneDrugOneFeature(self.ic50, self.input_features,
                drug='dummy', feature='dummy', fdr='dummy',
                directory=self.settings.directory)
        html.settings = self.settings
        html.factory.settings.includeMSI_factor = self.settings.includeMSI_factor
        html.factory.settings.analysis_type = self.settings.analysis_type

        self.html = html
        for i in range(N):
            html.drug = drugs[i]
            html.feature = features[i]
            html._filename = str(assocs[i]) + '.html'
            html.fdr = fdrs[i]
            html.assoc_id = assocs[i]
            html._init_report() # since we have one shared instance
            html.report(onweb=False)
            pb.animate(i+1)

    def create_html_features(self):
        """Create an HTML page for each significant feature"""
        df = self.get_significant_set()
        groups = df.groupby('FEATURE')
        print("\n\nCreating individual HTML pages for each feature")
        N = len(groups.indices.keys())
        pb = Progress(N)
        for i, feature in enumerate(groups.indices.keys()):
            # get the indices and therefore subgroup
            #if feature not in ['BRAF_mut', 'MLL2_mut', 'KRAS_mut']:
            #    continue
            subdf = groups.get_group(feature)
            metadata = {}
            metadata['n_cell_lines'] = self.input_features[feature].sum()
            # get concentration range
            metadata['feature'] = feature
            html = HTMLOneFeature(self.df, subdf, metadata,
                    directory=self.settings.directory)
            html.settings = self.settings
            self.subdf = subdf
            html.report(onweb=False)
            pb.animate(i+1)

    def create_html_drugs(self):
        """Create an HTML page for each drug that has at least one significant association"""
        # group by driugs
        df = self.get_significant_set()
        groups = df.groupby('Drug id')
        print("\n\nCreating individual HTML pages for each drug")
        N = len(groups.indices.keys())
        pb = Progress(N)
        for i, drug in enumerate(groups.indices.keys()):
            #if drug not in ['Drug_330_IC50', 'Drug_1047_IC50']:
            #    continue
            # get the indices and therefore subgroup
            subdf = groups.get_group(drug)
            metadata = {}
            metadata['n_cell_lines'] = len(self.ic50.df[drug].dropna())

            # get concentration range
            conc = subdf['log max.Conc.tested'].unique()[0]
            if pd.isnull(conc) is False:
                metadata['conc_min'] = conc / 4.**4
                metadata['conc_max'] = conc
            else:
                metadata['conc_min'] = "?"
                metadata['conc_max'] = "?"
            metadata['drug'] = drug

            html = HTMLOneDrug(self.df, subdf, metadata,
                    directory=self.settings.directory)
            html.settings = self.settings
            html.report(onweb=False)
            pb.animate(i+1)

    def create_html_main(self, onweb=False):
        """Create HTML main document (summary)"""
        print("\n\nCreating main HTML page in directory %s" %
                (self.settings.directory))
        buffer = self.settings.savefig
        self.settings.savefig = True
        html = HTML_main(self, 'index.html', directory=self.settings.directory)
        html._init_report() # created the directory
        html.settings = self.settings
        html.report(onweb=onweb)
        self.settings.savefig = buffer

    def create_html_manova(self):
        """Create summary table with all significant hits"""
        df = self.get_significant_set()
        html = HTMLManova(df, directory=self.settings.directory)
        html.report(onweb=False)

    def create_html_pages(self, onweb=False):
        """Create all HTML pages"""
        self._set_sensible_df()
        self.create_html_main(onweb=onweb)
        self.create_html_drugs()
        self.create_html_features()
        self.create_html_associations()
        self.create_html_manova()


class ANOVA(Logging):
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
    
    .. see also:: data format
    .. todo:: data format

    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import IC50, ANOVA, ic50_test
        ic = IC50(ic50_test)
        an = ANOVA(ic)
        # This is to select a specific tissue
        an.set_cancer_type('breast')
        df = an.anova_one_drug_one_feature('Drug_1047_IC50',
            'TP53_mut', show_boxplot=True)

    :Details about the anova analysis: In the example above, we perform a
        regression/anova test based on OLS regression. This is done for
        one feature one drug across all cell lines (tissue) in the method
        :meth:`anova_one_drug`. The regression
        takes into account the following factors: tissue, MSI and features. 
        The order matters. If there is only one tissue, this factor is 
        dropped. If the number of MSI values is less than a pre-defined 
        parameter (see :class:`ANOVASettings`), it is dropped. The other
        methods :meth:`anova_one_drug` and :meth:`anova_all` are wrappers 
        around :meth:`anova_one_drug_one_feature` to loop over all drugs, and
        loop over all drugs and all features, respectively.



    """
    def __init__(self, ic50, features=None, verbose='INFO'):
        """.. rubric:: Constructor

        :param DataFrame IC50: a dataframe with the IC50. Rows should be
            the COSMIC identifiers and columns should be the Drug names
            (or identifiers)
        :param features: another dataframe with rows as in the IC50 matrix
            and columns as features.  The first 3 columns must be named
            specifically to hold tissues, MSI (see format).
        :param verbose: verbosity in "WARNING", "ERROR", "DEBUG", "INFO"
            

        The attribute :attr:`settings` contains specific settings related
        to the analysis or visulation.
        """
        import scipy # placed here to speed up initialisation of gdsctools
        super(ANOVA, self).__init__(level=verbose)
        # Reads IC50
        self.logging.info('Reading data and building data structures')
        self.ic50 = readers.IC50(ic50)

        # Create a dictionary version of the data
        # to be accessed per drug where NA have already been
        # removed. Each drug is a dictionary with 2 keys:
        # Y for the data and indices for the cosmicID where
        # there is an IC50 measured.
        ic50_parse = self.ic50.df.copy().unstack().dropna()
        self.ic50_dict = dict([(d, {'indices': ic50_parse.ix[d].index,
            'Y':ic50_parse.ix[d].values}) for d in self.ic50.drugIds])

        # Reads features
        if features is None:
            # Reads default version provided with the package
            self.features = readers.GenomicFeatures()
        else:
            self.features = readers.GenomicFeatures(features)

        # settings
        self.settings = ANOVASettings()

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
        self._odof_dict = dict([(name, None)
            for name in self.column_names[1:]])

        # a cache to compute ANOVA
        self.individual_anova = {}

        self._init()

    def _autoset_msi(self):
        # if the number of positives(or negatives) factors is not large enough then
        # the MSI factor is not used
        self.msi_factor = self.features.df['MS-instability Factor Value']
        total = len(self.msi_factor)
        positives = self.msi_factor.sum()
        negatives = total - positives

        # we must have at least 2 positives or 2 negative
        # This is therefore a < comparison here below. See in
        # _get_one_drug_one_feature_data that we use >= which
        # is consistent
        if positives <  self.settings.MSIfactorPopulationTh:
            self.settings.includeMSI_factor = False
        # and 1 negative
        if negatives <  self.settings.MSIfactorPopulationTh:
            self.settings.includeMSI_factor = False

    def _autoset_tissue(self):
        # select tissue based on the features
        self.tissue_factor = self.features.df['Tissue Factor Value']
        if len(self.tissue_factor.unique()) == 1:
            # there is only one tissue
            tissue = self.tissue_factor.unique()[0]
            self.settings.analysis_type = tissue
        else:
            # this is a PANCAN analysis
            self.settings.analysis_type = 'PANCAN'

    def set_cancer_type(self, ctype=None):
        """Select only a set of tissues.
        
        Input IC50 may be PANCAN (several cancer tissues).
        This  function can be used to select a subset of tissues. 
        This function changes the :attr:`ic50` dataframe and possibly
        the feature as well is some are not relevant anymore (sum of the
        column is zero for instance).

        """
        if ctype is None:
            return
        ctype = easydev.to_list(ctype)
        for this in ctype:
           assert this in self.features.tissues


        # keep only features that correspond to the tissue
        # and have at least featFactorPopulationTh positives
        self.features.keep_tissue_in(ctype)

        self.ic50.df = self.ic50.df.ix[self.features.df.index]
        self._init()

    def _init(self):
        # Some preprocessing to speed up data access
        ic50_parse = self.ic50.df.copy().unstack().dropna()
        self.ic50_dict = dict([(d,
            {'indices': ic50_parse.ix[d].index,
             'Y':ic50_parse.ix[d].values}) for d in self.ic50.drugIds])

        # save the tissues
        self._autoset_tissue()

        # and MSI (Microsatellite instability) status of the samples.
        self._autoset_msi()

        # alias to speed up some code. Those are dictionary version of the
        # 3 dataframes above.
        self.features_dict = {}
        self.msi_dict = {}
        self.tissue_dict = {}

        # FIXME not sure we need a copy here. could be a reference if not
        # changed.
        for drug_name in self.ic50.drugIds:
            indices = self.ic50_dict[drug_name]['indices']
            self.features_dict[drug_name] = self.features.df.ix[indices].copy()
            self.msi_dict[drug_name] = self.msi_factor.ix[indices].copy()
            self.tissue_dict[drug_name] = self.tissue_factor.ix[indices].copy()

        # some preprocessing for the OLS compuation.
        # We create the dummies for the tissue factor once for all
        self._tissue_dummies = pd.get_dummies(self.tissue_factor)
        columns = self._tissue_dummies.columns
        columns = ['C(tissue)[T.' + x + ']' for x in columns]
        self._tissue_dummies.columns = columns
        N = len(self._tissue_dummies)
        self._tissue_dummies['C(msi)[T.1]'] = [1]*N
        self._tissue_dummies['feature'] = [1] * N
        self._tissue_dummies.insert(0, 'Intercept', [1] * N)

    def _get_cosmics(self):
        return self.ic50.cosmicIds
    def _set_cosmics(self, cosmics):
        self.ic50.cosmicIds = cosmics
        self.features.cosmicIds = cosmics
        self._init()
    cosmicIds = property(_get_cosmics, _set_cosmics,
        doc="get/set the cosmic identifiers in the IC50 and feature matrices")

    def _get_drug_names(self):
        return self.ic50.drugIds
    def _set_drug_names(self, drugs):
        self.ic50.drugIds = drugs
        self._init()
    drugIds = property(_get_drug_names, _set_drug_names)

    def _get_feature_names(self):
        return self.features.features
    def _set_features_names(self, features):
        self.features.features = features
        self._init()
    feature_names = property(_get_feature_names, _set_features_names)

    def _get_analysis_mode(self):
        modes = []
        if self.settings.analysis_type == 'PANCAN':
            modes.append('tissue')

        if self.settings.includeMSI_factor is True:
            modes.append('msi')

        modes.append('feature')
        return modes

    def diagnostics(self):
        """Return dataframe with information about the analysis

        173390 feasible tests in v17 (96.65)
        """
        n_drugs = len(self.ic50.drugIds)
        n_features = len(self.features.features) - 3
        n_combos = n_drugs * n_features
        feasible = 0
        pb = Progress(n_drugs, 1)
        counter = 0
        for drug in self.ic50.drugIds:
            for feature in self.features.features[3:]:
                status = self._get_one_drug_one_feature_data(drug, feature,
                        diagnostic_only=True)
                if status is True:
                    feasible += 1
            counter += 1
            pb.animate(counter)

        results = {
                'n_drug': n_drugs,
                'n_combos': n_combos,
                'feasible_tests': feasible,
                'percentage_feasible_tests': float(feasible)/n_combos*100}
        return results

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
        self.indices = indices
        # select only relevant tissues/msi/features
        # This line takes 50% of the time
        dd.masked_features = self.features_dict[drug_name][feature_name]
        dd.masked_tissue = self.tissue_dict[drug_name]
        dd.masked_msi = self.msi_dict[drug_name]

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

        # We could of course use the mean() and std() functions from pandas or
        # numpy. We could also use the glass and cohens modules but the
        # following code is now optimised to speed up this function 
        # call by 5/10 times.
        #dd.positives = dd.Y.values[dd.masked_features.values==1]
        #dd.negatives = dd.Y.values[dd.masked_features.values==0]

        dd.positives = dd.Y[dd.masked_features.values==1]
        dd.negatives = dd.Y[dd.masked_features.values==0]
        dd.Npos = len(dd.positives)
        dd.Nneg = len(dd.negatives)

        dd.A = A
        dd.B = B
        if (A == False) and (B == False):
            dd.status = False
            return dd
        else:
            dd.status = True

        if diagnostic_only is True:
            return dd.status

        # compute mean and std of pos and neg sets; using mean() takes 15us and
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

        csd = (dd.Npos -1.) * dd.pos_IC50_std**2 + \
                (dd.Nneg - 1.)*dd.neg_IC50_std**2
        csd /= dd.Npos + dd.Nneg -2.  # make sure this is float
        dd.effectsize_ic50 = md / np.sqrt(csd)

        # additional information
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
        if drug_name not in self.drugIds:
            raise ValueError('Unknown drug name %s. Use e.g., %s'
                    % (drug_name, self.drugIds[0]))
        if feature_name not in self.feature_names:
            raise ValueError('Unknown feature name %s. Use e.g., %s'
                    % (feature_name, self.feature_names[0]))
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
            #self._mydata = pd.DataFrame({'Y': odof.Y,
            #    'tissue':odof.masked_tissue,
            #    'msi':  odof.masked_msi, 'feature': odof.masked_features})
            #self.data_lm = ols('Y ~ C(msi) + feature',
            #    data=self._mydata).fit() #Specify C for Categorical

            # FIXME: 40% of the time is used to create this data structure
            # We could use pd.get_dummies but pretty slow
            # instead we create the full matrix of dummies and
            # then set relevant data. The issue is that some
            # columns end up with sum == 0 and needs to be dropped
            df = self._tissue_dummies.ix[odof.masked_tissue.index]
            todrop = df.columns[df.values.sum(axis=0)==0]
            df = df.drop(todrop, axis=1)
            #df = pd.get_dummies(odof.masked_tissue)
            #df.columns = ['C(tissue)[T.'+x +']' for x in
            #        odof.masked_tissue.unique()]
            # Ntissue = len(df.columns) - 3

            Ntissue = len(df.columns) - 3
            # Here we set other variables with dataframe columns' names as
            # expected by OLS
            df['C(msi)[T.1]'] = odof.masked_msi.values
            df['feature'] = odof.masked_features.values
            #df.insert(0, 'Intercept', [1] * (odof.Npos + odof.Nneg))

            # Here, we need to get rid of some of the cases to agree with R....
            # ??? why ?? TODO FIXME Could be that aov in R drops
            # some columns if number of positive is not large enough ?
            # dropping bladder across all give correct answers
            # across all drugs and features ?!
            # ?? if we remove skin instead of bladder, same results...
            # ?? if we remove 2 tissues, then we starts to have different
            # thinkgs???
            try:
                df = df.drop('C(tissue)[T.Bladder]', axis=1)
                Ntissue -= 1
            except:
                pass

            self.data_lm = OLS(odof.Y, df.values).fit()
            

            # SKLearn is also a possiblity
            # works for msi+feature but not if we include tissues ?
            # compared to R lm (not aov), we get the same
            # self.dff = df
            # ols = sklearn.linear_model.LinearRegression(fit_intercept=False)
            # ols.fit(an.dff, an.Y).coef_
            self.odof = odof
        elif self.settings.includeMSI_factor is True:
            #self._mydata = pd.DataFrame({'Y': odof.Y,
            #    'msi':  odof.masked_msi, 'feature': odof.masked_features})
            #self.data_lm = ols('Y ~ C(msi) + feature',
            #    data=self._mydata).fit() #Specify C for Categorical
            df = pd.DataFrame()
            df['C(msi)[T.1]'] = odof.masked_msi.values
            df['feature'] = odof.masked_features.values
            df.insert(0, 'Intercept', [1] * (odof.Npos + odof.Nneg))
            self.data_lm = OLS(odof.Y, df.values).fit()
            Ntissue = 0
        else:
            df = pd.DataFrame()
            df['feature'] = odof.masked_features.values
            df.insert(0, 'Intercept', [1] * (odof.Npos + odof.Nneg))
            self.data_lm = OLS(odof.Y, df.values).fit()
            Ntissue = 0
            #self._mydata = pd.DataFrame({'Y': odof.Y,
            #    'feature': odof.masked_features})
            #self.data_lm = ols('Y ~ feature',
            #    data=self._mydata).fit() #Specify C for Categorical
            #Ntissue = 0

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
                self._boxplot_pancan(odof, directory=directory,
                    savefig=savefig, fignum=3, mode='msi')


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
        #feature_name = odof.feature_name.replace("_", "\_")

        results = self._get_boxplot_data(odof, mode)
        if results is None:
            self.logging.info("INFO: no tissue with at least 2 pos and 2 neg found. " + "No image created.")
            return

        pylab.figure(fignum)
        pylab.clf()
        data, names, significance = results
        bb = boxswarm.BoxSwarm(data, names)
        bb.xlabel = r'%s log(IC50)' % drug_name
        if mode == 'tissue':
            bb.title = 'FEATURE/Cancer-type interactions'
        else:
            bb.title = 'FEATURE/MS-instability interactions'
        ax = bb.plot(vert=False)
        # get info from left axis
        common_ylim = ax.get_ylim()
        common_ticks = ax.get_yticks()

        self.ax = ax.twinx()
        self.ax.set_ylim(common_ylim)
        self.ax.set_yticks(common_ticks)
        self.ax.set_yticklabels([len(this) for this in data])
        #pylab.tight_layout()
        if savefig is True:
            filename = directory + os.sep
            filename += 'ODOF_{}_{}____{}'.format(mode,
                    odof.drug_name, odof.feature_name)
            pylab.savefig(filename + '.png', bbox_inches='tight')
            #pylab.savefig(filename + '.svg')

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

        #pylab.tight_layout()
        if savefig is True:
            filename = directory + os.sep
            filename += 'ODOF_all_{}____{}'.format(data.drug_name,
                    data.feature_name)
            pylab.savefig(filename + '.png', bbox_inches='tight')
            #pylab.savefig(filename + '.svg')

    def _get_anova_summary(self, data_lm, Ntissue, output='dict'):
        # could use this with statsmodels
        # The only values we want is PR(>F)
        # self.stats = sm.stats.anova_lm(self.data_lm, typ=1)
        # note that typ=1 does not work and typ=2,3 are differnt
        # from R version, that uses type I surely.
        q, r = np.linalg.qr(data_lm.model.data.exog)
        effects = np.dot(q.T, data_lm.model.data.endog)

        # create the W matrix using tissue and MSI if requested
        # default is that the 3 features are used
        modes = self._get_analysis_mode()
        Ncolumns = data_lm.model.data.exog.shape[1]

        if 'tissue' in modes and 'msi' in modes:
            dof = [Ntissue, 1, 1]
            indices = ['tissue', 'msi', 'feature', 'Residuals']
            # 4 stands for intercept + tissue + msi +feature
            arr = np.zeros((4, Ncolumns))
            arr[1, slice(1, Ntissue+1)] = 1
            arr[2, Ntissue + 1] = 1
            arr[3, Ntissue + 2] = 1
        elif 'tissue' not in modes and 'msi' in modes:
            dof = [1, 1]
            indices = ['msi', 'feature', 'Residuals']
            # 3 stands for intercept + msi +feature
            arr = np.zeros((3, Ncolumns))
            arr[1, Ntissue + 1] = 1
            arr[2, Ntissue + 2] = 1
        elif 'tissue' not in modes and 'msi' not in modes:
            dof = [1]
            indices = ['feature', 'Residuals']
            # 3 stands for intercept + msi +feature
            arr = np.zeros((2, Ncolumns))
            arr[1, Ntissue + 1] = 1
        arr[0, 0] = 1                   # intercept

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

        features = features[features.columns[3:]]
        mask = features.sum(axis=0) >= 3

        # TODO: MSI, tissues, name must always be kept
        #
        selected_features = features[features.columns[mask]]

        # scan all features for a given drug
        assert drug_id in self.ic50.df.columns
        N = len(selected_features.columns)
        pb = Progress(N, 10)
        res = {}
        # note that we start at idnex 4 to drop sample name, tissue and MSI
        for i, feature in enumerate(selected_features.columns):
            # production True, means we do not want to create a DataFrame
            # for each call to the anova_one_drug_one_feature function
            # Instead, we require dictionaries
            this = self.anova_one_drug_one_feature(drug_id, feature,
                    production=True)
            if this['FEATURE_ANOVA_pval'] is not None:
                res[feature] = this
            if animate is True:
                pb.animate(i+1)

        # if production is False:
        # df = pid.concat(res, ignore_index=True)
        df = pd.DataFrame.from_records(res)
        df = df.T

        df = ColumnTypes().astype(df)
        # TODO: drop rows where FEATURE_ANOVA_PVAL is None
        return df

    def anova_all(self, animate=True, drugs=None):
        """Run all ANOVA tests for all drugs and all features.

        :param drugs: you may select a subset of drugs

        Loops over all drugs calling :meth:`anova_one_drug` for each 
        drug and concatenating all results together. Note that once all
        data are gathered, an extra column containing the FDR corrections
        is added to the dataframe using :meth:`add_fdr_column` method.

        .. note:: comparison with version contained in this package
            gives same results. FDR (~1e-6) and FEATURE_IC50_T_pval differs
            slighlty (1e-14) especially for FDR variable with large FDR
            close to 1 but nothing to worry about.
        """
        # drop DRUG where number of IC50 (non-null) is below 5
        # axis=0 is default but we emphasize that sum is over 
        # column (i.e. drug
        vv = (self.ic50.df.isnull() == False).sum(axis=0)
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
            # TODO: try/except
            if drug_name in self.individual_anova.keys():
                pass
            else:
                res = self.anova_one_drug(drug_name, animate=False)
                self.individual_anova[drug_name] = res
            if animate is True:
                pb.animate(i+1)

        df = pd.concat(self.individual_anova, ignore_index=True)

        if len(df) == 0:
            return df
        # sort all data by ANOVA p-values
        try:
            df.sort_values('FEATURE_ANOVA_pval', inplace=True)
        except:
            df.sort('FEATURE_ANOVA_pval', inplace=True)

        # all ANOVA have been compute individually for each drug and each
        # feature. Now, we need to compute the FDR correction
        df = self.add_fdr_column(df)

        # insert a unique identifier as first column
        N = len(df)
        df.insert(0, 'assoc_id', range(1, N+1))

        # ordered as defined in the __init__ method
        df = df[self.column_names]
        df.reset_index(inplace=True)

        return df

    def add_fdr_column(self, df):
        """Add the FDR column in a datafream based on pvalues

        .. seealso:: :meth:`anova_all`
        """
        if len(df) == 0:
            return
        fdr = self._compute_fdr(df)
        # insert FDR as last column.
        try:
            df.insert(len(df.columns), 'ANOVA FEATURE FDR %', fdr)
        except:
            # replaces it
            df['ANOVA FEATURE FDR %'] = fdr
        return df

    def _compute_fdr(self, df):
        if self.settings.pval_correction_method == 'fdr':
            data = df['FEATURE_ANOVA_pval'].values
            fdr = fdrcorrection(data)[1]  * 100 # percentage ??
        else:
            raise NotImplementedError
            # should be qvalue correction (see qvalue library in R)
            fdr = [None] * len(df)
        return fdr

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
        # used in the query
        categories = list(cc.unstack().columns[cc])

        groups = df.query(mode + ' in @categories',
                engine='python').groupby([mode, 'feature'])

        # TODO; move all this if block into a method
        # figure out the delta between pos and neg
        means = groups.mean().unstack(mode)
        if len(means):
            delta = means.ix[0] - means.ix[1]
            try:
                # new pandas v0.17
                delta.sort_values(inplace=True)
            except:
                # sort_values not in anaconda for py3.3
                delta.sort(inplace=True)

            significance = {}
            data = []
            names = []
            for category in delta.ix['ic50'].index:
                prefix_query = mode+"==@category"
                neg = df.query(prefix_query+' and feature==0',
                        engine='python')['ic50']
                pos = df.query(prefix_query+' and feature==1',
                        engine='python')['ic50']
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

    def __str__(self):
        txt = self.ic50.__str__()
        txt += "\n" + self.features.__str__()
        return txt



def multicore(ic50, maxcpu=4):
    """Using 4 cores, the entire analysis took 15 minutes using
    4 CPUs (16 Oct 2015).

    :param ic50: a filename or :class:`IC50` instance.

    """

    import time
    t1 = time.time()
    master = ANOVA(ic50)

    drugs = master.ic50.drugIds

    from easydev import MultiProcessing
    t = MultiProcessing(maxcpu=maxcpu)
    # add all jobs (one per drug)
    for i, drug in enumerate(drugs):
        t.add_job(_analyse_one_drug, master, drug)
    t.run()

    # populate the ANOVA instance with the results
    for this in t.results:
        drug = this[0]
        result = this[1]
        master.individual_anova[drug] = result

    print("\nTook " + str(time.time() - t1) +  "seconds.")

    return master

def _analyse_one_drug(master, drug):
    res = master.anova_one_drug(drug_id=drug, animate=False)
    return (drug, res)


class HTMLManova(Report):
    """Creates an HTML page dedicated to significant hits

    ::

        # analyse the data across all drugs and features
        df = anova.anova_all()
        # Create a table for the first 10 significant hits (table is sorted
        # by ascending FDR)
        h = HTMLManova(df.ix[0:10, directory='.')
        # Create the HTML page (pops up by default)
        h.report()

    """
    def __init__(self, df, directory='gdsc'):
        """

        :param : a dataframe as output by :meth:`ANOVA.anova_all`
        :param directory: where to save the file

        The HTML filename is stored in the :attr:`filename`, which can
        be changes (default is manova.html)
        """
        self.df = df
        self._filename = 'manova.html'
        super(HTMLManova, self).__init__(directory=directory,
                filename=self.filename)

    def _create_report(self):
        hits = SignificantHits(self.df, 'all hits')
        html = hits.to_html()
        self.add_section(html, 'MANOVA results (significant only) summary')


class OneDrugOneFeature(Report):
    def __init__(self, ic50, features=None, drug=None, feature=None,
            directory='gdsc', fdr='?', assoc_id='?'):
        # FIXME here we lose the setttings since we create a new instance
        self.factory = ANOVA(ic50, features=features)
        self.factory.settings.directory = directory
        self.assoc_id = assoc_id

        self.drug = drug
        self.feature = feature
        self.fdr = fdr
        self.add_settings = False
        filename = "{0}____{1}.html".format(self.drug,
                self.feature.replace(" ", "_"))

        super(OneDrugOneFeature, self).__init__(directory=directory,
                filename=filename)

    def run(self):
        #pylab.ioff()
        df = self.factory.anova_one_drug_one_feature(self.drug,
                self.feature, savefig=True, show_boxplot=True,
                directory=self.directory)
        #pylab.ion()
        # FIXME assoc id
        df.insert(0, 'assoc_id', self.assoc_id)
        df['ANOVA FEATURE FDR %'] = self.fdr
        return df

    def to_html(self, df, precision=6):
        sign = SignificantHits(df, 'features')
        html = sign.to_html(escape=False, header=True, index=False)
        return html

    def _create_report(self, onweb=True):
        # generated pictures and results
        #print('Generating data, images and HTML')
        df = self.run()

        # Create the table and add it
        html_table = self.to_html(df, precision=2)
        self.add_section(html_table, 'Individual association analysis')

        section = ""
        for prefix in ['ODOF_all', 'ODOF_msi', 'ODOF_tissue']:
            tag = "{0}_{1}____{2}.png".format(prefix, self.drug, self.feature)
            section += '<img alt="association {0}" src="{0}">\n'.format(tag)
        self.add_section(section, "Boxplots")

        if self.add_settings is True:
            table = ANOVASettings(**self.factory.settings)
            self.add_section(table.to_html(), 'Settings')



class HTMLOneFeature(Report):
    def __init__(self, data, subdata, metadata, directory='gdsc'):
        self.df = data
        self.subdf = subdata
        self.feature = metadata['feature']
        self.metadata = metadata
        filename = "{0}.html".format(self.feature)
        super(HTMLOneFeature, self).__init__(directory=directory,
                filename=filename)

    def run(self, N=20):
        v = VolcanoANOVA(self.df)
        v.settings = self.settings # get fdr, pval
        v.settings.savefig = True
        v.settings.directory = self.directory
        # FIXME: insiced volvano_plot, we add tooltip
        # when we call volcano again, the tooltip form the
        # previous call are still there. The only solution
        # so far is to close the figure.
        pylab.close(1)
        v.volcano_plot_one_feature(self.feature)
        v.savefig('volcano_{}.png'.format(self.feature))

    def to_html(self, df, precision=6):
        sign = SignificantHits(df, 'features')
        html = sign.to_html(escape=False, header=True, index=False)
        return html

    def _create_report(self, onweb=True):
        # generated pictures and results
        #print('Generating data, images and HTML')
        self.run()
        # could be a table with no border ?
        summary = """
        Binary feature equal to %(binary)s for samples harboring mutations in MLL2
        Number of cell lines positive for this feature: %(n_cell_lines)s
        <br>
        """
        self.metadata['binary'] = '?'

        N = len(self.subdf)

        # Just a paragraph to give the drug name
        self.add_pretoc("<p><br>Individual Feature analysis :<br> {0}</p>".format(self.feature))

        # Intro section
        section = summary % self.metadata
        self.add_section(section, "Feature and screening numbers")

        # Table section
        section = "Statistically significant associations identified = {0}".format(N)
        if N>0:
            section += self.to_html(self.subdf, precision=2)
        self.add_section(section, 'Individual association analysis')

        # image section
        section = '<img alt="volcano plot {0}" src="volcano_{0}.png">\n'.format(self.feature)

        self.add_section(section, "All-tests volcano plot")



class HTMLOneDrug(Report):
    def __init__(self, data, subdata, metadata,
            directory='gdsc'):
        """
        data is a dataframe with all results from anova_all
        subdata is a dataframe with significant associations. Can be empty
        metadata is a dict with those keys:
            n_cell_lines
            min_conc
            max_conc
            drug
        """
        self.df = data
        self.subdf = subdata
        self.drug = metadata['drug']
        self.metadata = metadata
        filename = "{0}.html".format(self.drug)
        super(HTMLOneDrug, self).__init__(directory=directory,
                filename=filename)

    def create_pictures(self):
        v = VolcanoANOVA(self.df)
        v.settings = self.settings # get fdr, pval
        v.settings.savefig = True
        v.settings.directory = self.directory
        v.volcano_plot_one_drug(self.drug)
        v.savefig('volcano_{}.png'.format(self.drug))

    def to_html(self, df, precision=6):
        # why ?
        #try:df.drop('log max.Conc.tested')
        #except: pass
        #try: df.drop('log max.Conc.tested2')
        #except:pass

        sign = SignificantHits(df, 'drugs')
        html = sign.to_html(escape=False, header=True, index=False)
        return html

    def _create_report(self, onweb=True):
        # generated pictures and results
        #print('Generating data, images and HTML')
        self.create_pictures()
        # could be a table with no border ?
        summary = """
        Drug Name: %(drug_name)s<br>
        Drug ID: %(drug_id)s<br>
        Synonyms: %(synonyms)s<br>
        Brand name: %(brand_name)s<br>
        Target: %(target)s<br>
<br>
        Number of cell lines screeened: %(n_cell_lines)s<br>
        Screening concentration range(2): %(conc_min)s to %(conc_max)s uM<br>
        """

        N = len(self.subdf)
        if N >= 0:
            # add the table
            self.metadata['drug_name'] = self.subdf['Drug name'].unique()[0]
            self.metadata['drug_id'] = self.subdf['Drug id'].unique()[0]
            self.metadata['synonyms'] = ''
            self.metadata['brand_name'] = ''
            self.metadata['target'] = self.subdf['Drug Target'].unique()[0]
        else:
            pass

        # Just a paragraph to give the drug name
        self.add_pretoc("<p><br>Drug id :<br> {0}</p>".format(self.drug))

        # Intro section
        section = summary % self.metadata
        self.add_section(section, "Drug details and screening numbers")

        # Table section
        section = "Statistically significant associations identified = {0}".format(N)
        if N>0:
            section += self.to_html(self.subdf, precision=2)
        self.add_section(section, 'Individual association analysis')

        # image section
        section = '<img alt="volcano plot" src="volcano_{0}.png">\n'.format(self.drug)

        self.add_section(section, "All-tests volcano plot")


class HTML_main(Report):
    def __init__(self, results, filename='index.html', directory='gdsc'):
        super(HTML_main, self).__init__(directory=directory,
                filename=filename)
        self.results = results
        self.add_dependencies = True
        self.settings = self.results.settings

    def _create_report(self, onweb=True):
        # A summary table
        diag = self.results.diagnostics()
        table = HTMLTable(diag, 'summary')
        txt = ''
        for index, row in diag.iterrows():
            txt += row.text + ": " +  str(row.value) + "<br/>"
        self.add_section(txt, 'Summary')

        print('Create summary plots')
        # this can be pretty slow. so keep only 1000 most relevant
        # values and 1000 random ones to get an idea of the distribution
        v = VolcanoANOVA(self.results.df)
        v.selector(v.df, 1000,1000, inplace=True)

        v.settings = self.settings # get fdr, pval
        v.settings.savefig = True
        v.settings.directory = self.directory
        v.volcano_plot_all()
        # volcano plot
        html = """
<h3></h3>
<img alt="volcanot plot for all associations" src="volcano_all.png">

<br>
Possibly, a javascript version is available
<a href="volcano_all_js.html">here</a>

        """
        try:
            import mpld3
            htmljs = mpld3.fig_to_html(v.current_fig)
        except:
            htmljs = ""
        fh = open(self.directory + os.sep + "volcano_all_js.html","w")
        fh.write(htmljs)
        fh.close()

        self.add_section(html, 'Volcano plot')

        # feature summary
        df_features = self.results.feature_summary()
        filename = 'features_summary.tsv'
        df_features.to_csv(self.directory + os.sep + filename, sep='\t')
        html = """
<h3>Features most frequently associated with a drug response</h3>
<img alt="summary of significant features" src="feature_summary.png">
<br>
You can <a href="{}">download the significant-features table</a> in tsv format.
""".format(filename)
        self.add_section(html, 'Feature summary')

        # MANOVA link
        self.add_section('Explore all significant associations following this link <a href="manova.html">manova</a>', "MANOVA")

        # drug summary
        df_drugs = self.results.drug_summary()
        filename = 'drugs_summary.tsv'
        df_drugs.to_csv(self.directory + os.sep + filename, sep='\t')
        html = """
<h3>Drug whose response is frequently associated with afeature</h3>
<img alt="summary of significant drugs" src="drug_summary.png">
<br>You can <a href="{}">download the significant-features table</a> in tsv format.
""".format( filename)

        self.add_section(html, 'Drug summary')

        # Create table with links to all drugs
        groups = self.results.df.groupby('Drug id')
        try:
            try:
                df = groups.mean()['ANOVA FEATURE FDR %'].sort_values()
            except:
                df = groups.mean()['ANOVA FEATURE FDR %'].sort()
            df = df.reset_index() # get back the Drug id in the dataframe columns
            # add another set of drug_id but sorted in alpha numerical order
            drugs = list(df['Drug id'].values)
            alphanum = lambda x: int(x.split("_")[1])
            drugs.sort(key=alphanum)
            df['Drug id (alpha order)'] = drugs
            table = HTMLTable(df, 'drugs')
            table.add_href('Drug id')
            table.add_href('Drug id (alpha order)')
        except:
            # FIXME 
            table = HTMLTable(self.results.df, 'drugs')
            table.add_href('Drug id')
        # rename one columns
        table.df.columns = [x.replace('ANOVA FEATURE FDR',
            'mean ANOVA FEATURE FDR') for x in table.df.columns]


        html = "The following table provides links to dedicated pages for each drug (sorted by ascending FDR)"
        html += table.to_html(escape=False, header=True, index=False)
        self.add_section(html, 'Drug wise associations browse')

        # Create full table with links to all features
        df = pd.DataFrame({'FEATURE': self.results.df['FEATURE'].unique()})
        try:
            df.sort_values(by='FEATURE', inplace=True)
        except:
            df.sort('FEATURE', inplace=True)

        groups = self.results.get_significant_set().groupby('FEATURE').groups

        count = []
        for feature in df['FEATURE'].values:
            if feature in groups.keys():
                count.append(len(groups[feature]))
            else:
                count.append(0)
        df['hits'] = count

        table = HTMLTable(df, 'features')
        table.add_href('FEATURE')
        html = table.to_html(escape=False, header=True, index=False)
        self.add_section(html, 'Feature wise associations browse')

        # Save the settings
        self.add_section(self.settings.to_html(), 'Settings')


class SignificantHits(object):
    def __init__(self, df, name):
        self.df = df.copy() # to not alter original version
        self.cmap_clip = cmap_builder('#ffffff', '#0070FF')
        self.cmap_absmax = cmap_builder('green', 'white', 'red')
        self.name = name
        self.clip_threshold = 2

        columns = [u'assoc_id', 'FEATURE',
            'Drug id', u'Drug name', 'Drug Target',
            'N_FEATURE_neg', 'N_FEATURE_pos',
            'FEATUREpos_logIC50_MEAN',
            'FEATUREneg_logIC50_MEAN',
            'FEATURE_deltaMEAN_IC50',
            'FEATURE_IC50_effect_size',
            'FEATUREneg_Glass_delta',
            'FEATUREpos_Glass_delta',
            'FEATURE_ANOVA_pval',
            'Tissue_ANOVA_pval',
            'MSI_ANOVA_pval',
            'ANOVA FEATURE FDR %']
        self.df = self.df[columns]

    def to_html(self, escape=False, header=True, index=False):
        # If there is a value below 0.01, a scientific notation is
        # use. we prefer to use 2 digits and write <0.01 
        colname = 'ANOVA FEATURE FDR %'
        self.df.loc[self.df[colname] <0.01, colname] = '<0.01'
        html = HTMLTable(self.df, self.name)
        # Those columns should be links
        for this in ['FEATURE', 'Drug id', 'assoc_id']:
            html.add_href(this)


        for this in ['FEATURE_IC50_effect_size', 'FEATUREneg_Glass_delta',
                'FEATUREpos_Glass_delta']:
            html.add_bgcolor(this, self.cmap_clip, mode='clip',
                    threshold=self.clip_threshold)

        # normalise data and annotate with color
        html.add_bgcolor('FEATURE_deltaMEAN_IC50', self.cmap_absmax,
                mode='absmax')

        html.df.columns = [x.replace("_", " ") for x in html.df.columns]
        return html.to_html(escape=escape, header=header, index=index,
                justify='center')
