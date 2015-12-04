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
import os
import shutil
import pandas as pd
import scipy
import pylab
import numpy as np

from statsmodels.formula.api import OLS

from easydev import Progress, AttrDict, Logging
import easydev
from colormap import cmap_builder

from gdsctools.stats import MultipleTesting
from gdsctools import readers
from gdsctools.boxplots import BoxPlots
from gdsctools.report import HTMLTable, ReportMAIN
from gdsctools.tools import Savefig
from gdsctools.volcano import VolcanoANOVA
from gdsctools.settings import ANOVASettings
from gdsctools.anova_results import ANOVAResults

__all__ = ['ANOVA', 'ANOVAReport']



class ANOVAReport(object):
    """Class used to interpret the results and create final HTML report

    Results is a data structure returned by :meth:`ANOVA.anova_all`.

    ::

        from gdsctools import *

        # Perform the analysis itself to get a set of results (dataframe)
        an = ANOVA(ic50_test)
        results = an.anova_all()

        # now, we can create the report.
        r = ANOVAReport(gdsc=an, results=results)

        # we can tune some settings
        r.settings.pvalue_threshold = 0.001
        r.settings.FDR_threshold = 28
        r.settings.directory = 'testing'
        r.create_html_pages()

    :Significant association: an association is significant if one of those
        conditions is fulfilled:
         - The field *ANOVA_FEATURE_FDR_%* must be < FDR_threshold
         - The field *FEATURE_ANOVA_pval* must be < pvalue_threshold
         - The field *FEATURE_deltaMEAN_IC50* must be < 0 (sensible) or
           >= 0 (resistant)
        It can be **resistant** or **sensitive**.

    """
    def __init__(self, gdsc, results, sep="\t", drug_decoder=None):
        """.. rubric:: Constructor

        :param gdsc: the instance with which you created the results to report
        :param results: the results returned by :meth:`ANOVA.anova_all`

        """
        self.figtools = Savefig()
        self.gdsc = gdsc
        self.df = ANOVAResults(results).df # this does a copy and sanity check

        self.settings = ANOVASettings()
        for k, v in gdsc.settings.items():
            self.settings[k] = v

        self._colname_drug_id = 'DRUG_ID'
        self.varname_pval = 'FEATURE_ANOVA_pval'
        self.varname_qval = 'ANOVA_FEATURE_FDR_%'

        # maybe there was not drug_decoder in the gdsc parameter,
        # so a user may have provide a file, in which case, we need
        # to update the content of the dur_decoder.
        if len(gdsc.drug_decoder) == 0 and drug_decoder is None:
            print('\nWARNING no drug name or target will be populated')
            print('You can read one if you wish using read_drug_decoder')
        elif drug_decoder is not None:
            self.read_drug_decoder(drug_decoder)
        else: # should be in gdsc.drug_decoder
            pass

        """if concentrations:
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
        """
        # create some data
        self._set_sensible_df()

        # just to create the directory
        ReportMAIN(directory=self.settings.directory)

    def _get_ndrugs(self):
        return len(self.df[self._colname_drug_id].unique())
    n_drugs = property(_get_ndrugs, doc="return number of drugs")

    def read_drug_decoder(self, filename):
        """Read file with the DRUG information

        .. seealso:: :class:`gdsctools.readers.DrugDecoder`
        """
        if filename is not None:
            self.gdsc.read_drug_decoder(filename)
            self.df = self.gdsc.drug_annotations(self.df)

    def _get_ntests(self):
        return len(self.df.index)
    n_tests = property(_get_ntests)

    def _get_ncelllines(self):
        return len(self.gdsc.features.df.index)
    n_celllines = property(_get_ncelllines,
            doc="return number of cell lines")

    def _df_append(self, df, data):
        count = len(df)
        df.ix[count] = data
        return df

    def diagnostics(self):
        """Return summary of the analysis (dataframe)"""
        df = pd.DataFrame({'text': [], 'value': []})

        n_features = len(self.gdsc.features.df.columns)
        n_features -= self.gdsc.features.shift
        n_drugs = len(self.df[self._colname_drug_id].unique())

        N = float(n_drugs * n_features)
        ratio = float(self.n_tests)/(N) * 100
        ratio = easydev.precision(ratio, digit=2)

        msg = "Type of analysis"
        df = self._df_append(df, [msg, self.settings.analysis_type])

        msg = "Total number of possible drug/feature associations"
        df = self._df_append(df, [msg, int(N)])
        msg = "Total number of ANOVA tests performed"
        df = self._df_append(df, [msg, self.n_tests])
        msg = "Percentage of tests performed"
        df = self._df_append(df, [msg, ratio])

        # trick to have an empty line
        df = self._df_append(df, ["", ""])

        msg = "Total number of tested drugs"
        df = self._df_append(df, [msg, n_drugs])
        msg = "Total number of genomic features used"
        df = self._df_append(df, [msg, n_features])

        msg = "Total number of screened cell lines"
        df = self._df_append(df, [msg, self.n_celllines])

        msg = "MicroSatellite instability included as factor"
        msi = self.settings.include_MSI_factor
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
        value = "[{:.4}, {:.4}]".format(p1, p2)
        df = self._df_append(df, [msg, value])

        f1, f2 = self._get_fdr_range()
        msg = "Range of significant % FDRs"
        value = '[{:.4} {:.4}]'.format(f1, f2)
        df = self._df_append(df, [msg, value])
        return df

    def _get_pval_range(self):
        """Get pvalues range of the significant hits"""
        nsens = len(self.sensible_df)
        nres = len(self.resistant_df)
        N = nsens + nres
        if N == 0:
            return 0., 0.
        name = self.varname_pval
        data = self.df[name].ix[0:N-1]
        m, M = data.min(), data.max()
        return m, M

    def _get_fdr_range(self):
        """Get FDR range of the significant hits"""
        name = self.varname_qval
        data = self.df[name][(self.df[name] < self.settings.FDR_threshold)]
        if len(data) == 0:
            return 0., 0.
        m, M = data.min(), data.max()
        return m, M

    def _set_sensible_df(self):
        # just an alias
        logand = np.logical_and

        # select sensible data set
        mask1 = self.df['ANOVA_FEATURE_FDR_%'] < self.settings.FDR_threshold
        mask2 = self.df['FEATURE_ANOVA_pval'] < self.settings.pvalue_threshold
        mask3 = self.df['FEATURE_deltaMEAN_IC50'] < 0
        self.sensible_df = self.df[logand(logand(mask1, mask2), mask3)]

        # select resistant data set
        mask3 = self.df['FEATURE_deltaMEAN_IC50'] >= 0
        self.resistant_df = self.df[logand(logand(mask1, mask2), mask3)]

    def get_significant_set(self):
        """Return significant hits (resistant and sensible)"""
        # a property that is long to compute
        # and may change if FDR changes.
        self._set_sensible_df()
        df = pd.concat([self.sensible_df, self.resistant_df])
        try:
            df.sort_values('ASSOC_ID', inplace=True)
        except:
            df.sort('ASSOC_ID', inplace=True)
        return df

    def _get_data(self, df_count_sensible, df_count_resistant):
        # we can drop all columns except one, which is renamed as count
        df1 = df_count_sensible['ASSOC_ID']
        df1.name = 'sens assoc'
        df2 = df_count_resistant['ASSOC_ID']
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

    def drug_summary(self,  top=50, fontsize=15, filename=None):
        """Return dataframe with significant drugs

        :param fontsize:
        :param top: max number of significant associations to show
        :param filename: if provided, save the file in the directory
        """
        # get sensible and resistant sub dataframes
        self._set_sensible_df()

        # group by drug
        colname = self._colname_drug_id
        df_count_sensible = self.sensible_df.groupby(colname).count()
        df_count_resistant = self.resistant_df.groupby(colname).count()

        df_count = self._get_data(df_count_sensible, df_count_resistant)

        if len(df_count):
            self._plot(df_count, 'drug', top)
            fig = pylab.gcf()
            self.figtools.directory = self.settings.directory
            self.figtools.savefig(filename, size_inches=(12, 14), 
                    bbox_inches='tight')
        return df_count

    def feature_summary(self, filename=None, top=50, fontsize=15):
        """Return dataframe with significant features

        :param fontsize:
        :param top: max number of significant associations to show
        :param filename: if provided, save the file in the directory
        """
        # get sensible and resistant sub dataframes
        self._set_sensible_df()

        df_count_sensible = self.sensible_df.groupby('FEATURE').count()
        df_count_resistant = self.resistant_df.groupby('FEATURE').count()

        df_count = self._get_data(df_count_sensible, df_count_resistant)

        if len(df_count) > 0:
            self._plot(df_count, 'feature', top)
            #fig = pylab.gcf()
            self.figtools.directory = self.settings.directory
            self.figtools.savefig(filename, set_inches=(12, 14),
                    bbox_inches='tight')
        return df_count

    def _plot(self, df_count, title_tag, top):
        """Used by drug_summary and feature_summary to plot the
        bar plot"""
        if top > len(df_count):
            top = len(df_count)

        df = df_count.ix[0:top][[u'sens assoc', u'res assoc']]
        labels = list(df.index)
        # add drug name
        if len(self.gdsc.drug_decoder) > 0:
            for i, label in enumerate(labels):
                name = self.gdsc.drug_decoder.get_name(label)
                if name is not None:
                    labels[i] = labels[i] + " - " + name
                else:
                    pass

        labels = [x.replace('_', ' ') for x in labels]
        # restrict size to first 30 characters
        labels = [x[0:30] for x in labels]
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
        ax.set_yticks([x + 0.5 for x in ind])
        ax.set_yticklabels(labels, fontsize=12)
        pylab.grid()
        pylab.title(r"Top %s %s most frequently " % (top, title_tag) + \
                    "\nassociated with drug  response", 
                    fontsize=self.settings.fontsize/1.2)
        pylab.xlabel(r'Number of significant associations (FDR %s %s %s) '
                     % ("$>$", self.settings.FDR_threshold, "$\%$"),
                     fontsize=18)

        M = max(data1+data2)
        #ax.set_xticks()
        #ax.set_xticklabels(labels, fontsize=fontsize)
        ax.set_xlim([0, M+1])
        pylab.legend(loc='lower right')
        try:pylab.tight_layout()
        except:pass

    def get_significant_hits(self,  show=True):
        """Return a summary of significant hits

        :param show: show a plot with the distribution of significant hits

        .. todo:: to finalise
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
            res = self.df['ANOVA_FEATURE_FDR_%']<fdr
            significants.append(res.sum())

            # meaningful hits
            indices = np.logical_and(self.df['ANOVA_FEATURE_FDR_%']<fdr,
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

    def __str__(self):
        self.df.info()
        return ""

    def create_html_associations(self):
        """Create an HTML page for each significant association

        The name of the output HTML file is **<association id>.html**
        where association id is stored in :attr:`df`.

        """
        print("\n\nCreating individual HTML pages for each association")
        df = self.get_significant_set()

        drugs = df['DRUG_ID'].values
        features = df['FEATURE'].values
        assocs = df['ASSOC_ID'].values
        fdrs = df['ANOVA_FEATURE_FDR_%'].values

        N = len(df)
        pb = Progress(N)

        html = Association(self, drug='dummy', feature='dummy', 
                fdr='dummy')

        for i in range(N):
            html.drug = drugs[i]
            html.feature = features[i]
            html._filename = str(assocs[i]) + '.html'
            html.fdr = fdrs[i]
            html.assoc_id = assocs[i]
            html._init_report() # since we have one shared instance
            html.create_report(onweb=False)
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
            subdf = groups.get_group(feature)
            html = HTMLOneFeature(self, self.df, subdf, feature)
            html.create_report(onweb=False)
            pb.animate(i+1)

    def create_html_drugs(self):
        """Create an HTML page for each drug that has at
        least one significant association

        Actually, we are interested in each drug, could be a flag

        """
        # group by drugs
        all_drugs = list(self.df['DRUG_ID'].unique())

        df = self.get_significant_set()
        groups = df.groupby('DRUG_ID')
        print("\n\nCreating individual HTML pages for each drug")
        N = len(groups.indices.keys())
        N = len(all_drugs)
        pb = Progress(N)
        #all_drugs = list(df.DRUG_ID) + ['Drug_1050_IC50']
        for i, drug in enumerate(all_drugs):
            # enumerate(groups.indices.keys()):
            # get the indices and therefore subgroup
            if drug in groups.groups.keys():
                subdf = groups.get_group(drug)
            else:
                subdf = {}

            html = HTMLOneDrug(self, self.df, subdf, drug)
            html.create_report(onweb=False)
            pb.animate(i+1)

    def create_html_main(self, onweb=False):
        """Create HTML main document (summary)"""
        print("\n\nCreating main HTML page in directory %s" %
                (self.settings.directory))
        buffer_ = self.settings.savefig
        self.settings.savefig = True
        html = HTMLPageMain(self, 'index.html')
        html._init_report() # created the directory
        html.create_report(onweb=onweb)
        self.settings.savefig = buffer_

    def create_html_manova(self, onweb=True):
        """Create summary table with all significant hits
        
        :param onweb: open browser with the created HTML page.

        """
        df = self.get_significant_set()
        page = HTMLPageMANOVA(self.gdsc, df)
        page.create_report(onweb)

    def create_html_pages(self, onweb=True):
        """Create all HTML pages"""
        self._set_sensible_df()
        self.create_html_main(onweb=onweb)
        self.create_html_drugs()
        self.create_html_features()
        self.create_html_associations()
        self.create_html_manova(onweb=False)

    def onweb(self):
        from easydev import onweb
        onweb(self.settings.directory + os.sep + 'index.html')


# Not that Logging is not used: it is not pickable and prevent
# multicore analysis.
class ANOVA(object): #Logging):
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

    """
    def __init__(self, ic50, genomic_features=None, 
            drug_decoder=None, verbose='INFO', low_memory=False):
        """.. rubric:: Constructor

        :param DataFrame IC50: a dataframe with the IC50. Rows should be
            the COSMIC identifiers and columns should be the Drug names
            (or identifiers)
        :param features: another dataframe with rows as in the IC50 matrix
            and columns as features.  The first 3 columns must be named
            specifically to hold tissues, MSI (see format).
        :param drug_decoder: a 3 column CSV file with drug's name and targets
            see :mod:`readers` for more information.
        :param verbose: verbosity in "WARNING", "ERROR", "DEBUG", "INFO"

        The attribute :attr:`settings` contains specific settings related
        to the analysis or visulation.
        """
        # We first need to read the IC50 using a dedicated reader
        self.ic50 = readers.IC50(ic50)

        # Create a dictionary version of the data
        # to be accessed per drug where NA have already been
        # removed. Each drug is a dictionary with 2 keys:
        # Y for the data and indices for the cosmicID where
        # there is an IC50 measured.
        ic50_parse = self.ic50.df.copy().unstack().dropna()
        self.ic50_dict = dict([(d, {'indices': ic50_parse.ix[d].index,
            'Y':ic50_parse.ix[d].values}) for d in self.ic50.drugIds])

        # Reads features if provided, otherwise use a default data set
        if genomic_features is None:
            # Reads default version provided with the package
            self.features = readers.GenomicFeatures()
        else:
            self.features = readers.GenomicFeatures(genomic_features)

        #: a CSV with 3 columns used in the report
        self.read_drug_decoder(drug_decoder)

        # create the multiple testing factory used in anova_all()
        self.multiple_testing = MultipleTesting()

        # We prune the genomic features by settings the cosmic ids of
        # the features to be those of the cosmic ids of the IC50. See
        # readers module. This affectation, prune the features dataframe
        # automatically. This fails if a cosmic identifier is not
        # found in the features' cosmic ids, so let us catch the error
        # before hand to give a
        unknowns = set(self.ic50.cosmicIds).difference(
                set(self.features.cosmicIds))
        if len(unknowns) > 0:
            print("WARNING: " +
                "%s cosmic identifiers in your IC50 " % len(unknowns) +
                "could not be found in the genomic feature matrix. "+
                "They will be dropped. Consider using a user-defined " +
                "genomic features matrix")

        self.ic50.drop_cosmic(list(unknowns))
        self.features.cosmicIds  = self.ic50.cosmicIds
        #self.cosmicIds = self.ic50.cosmicIds

        #: an instance of :class:`~gdsctools.settings.ANOVASettings`
        self.settings = ANOVASettings()
        self.settings.low_memory = low_memory

        # alias to all column names to store results (unordered)
        # cast to list because of  python3 error.
        self.column_names = list(ANOVAResults().mapping.keys())

        # skip assoc_id for now
        self._odof_dict = dict([(name, None)
            for name in self.column_names])

        # a cache to store ANOVA results for each drug
        self.individual_anova = {}

        # must be called if ic50 or features are changed.
        self._init()

    def _autoset_msi(self):
        if self.features.found_msi:
            # if the number of pos. (or neg.) factors is not large enough then
            # the MSI factor is not used
            msi_name = self.features.colnames.msi
            self.msi_factor = self.features.df[msi_name]
            total = len(self.msi_factor)
            positives = self.msi_factor.sum()
            negatives = total - positives

            # we must have at least 2 positives or 2 negative
            # This is therefore a < comparison here below. See in
            # _get_one_drug_one_feature_data that we use >= which
            # is consistent.
            if positives < self.settings.MSI_factor_threshold:
                self.settings.include_MSI_factor = False
            if negatives < self.settings.MSI_factor_threshold:
                self.settings.include_MSI_factor = False
        else:
            self.settings.include_MSI_factor = False

    def _autoset_tissue(self):
        # select tissue based on the features
        tissue_name = self.features.colnames.tissue
        self.tissue_factor = self.features.df[tissue_name]
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
        the feature as well if some are not relevant anymore (sum of the
        column is zero for instance).

        """
        if ctype is None:
            return
        # input may be unicode. transform to str otherwise, error pandas.
        ctype = str(ctype)

        if ctype == 'PANCAN':
            return


        if isinstance(ctype, str):
            ctype = [ctype]

        for this in ctype:
            assert this in self.features.tissues, "%s not found" % ctype

        # keep only features that correspond to the tissue
        self.features.keep_tissue_in(ctype)

        self.ic50.df = self.ic50.df.ix[self.features.df.index]
        self._init()

    def read_settings(self, settings):
        """Read settings and update cancer type if set"""
        self.settings.from_json(settings)
        self.set_cancer_type(self.settings.analysis_type)

    def _init(self):
        # Some preprocessing to speed up data access
        ic50_parse = self.ic50.df.copy().unstack().dropna()
        # for each drug, we store the IC50s (Y) and corresponding indices
        # of cosmic identifiers
        self.ic50_dict = dict([
            (d, {'indices': ic50_parse.ix[d].index,
             'Y': ic50_parse.ix[d].values}) for d in self.ic50.drugIds])

        # save the tissues
        self._autoset_tissue()

        # and MSI (Microsatellite instability) status of the samples.
        self._autoset_msi()

        # dictionaries to speed up code.
        self.features_dict = {}
        self.msi_dict = {}
        self.tissue_dict = {}
        # fill the dictionaries for each drug once for all
        for drug_name in self.ic50.drugIds:
            indices = self.ic50_dict[drug_name]['indices']
            # if we were to store all drugs /features, this takes
            # 1Gb of memory for 265 drugs and 680 features. This is
            # therefore not scalable, especially for multiprocessing.
            if self.settings.low_memory is True:
                pass
            else:
                self.features_dict[drug_name] = self.features.df.ix[indices]

            # MSI and tissue can be store
            if self.features.found_msi:
                self.msi_dict[drug_name] = self.msi_factor.ix[indices]
            self.tissue_dict[drug_name] = self.tissue_factor.ix[indices]

        # some preprocessing for the OLS computation.
        # We create the dummies for the tissue factor once for all
        # to agree with R convention, we sort the column
        # by alpha order where a<B==b<c unlike in R where A<B<C<a<b<c
        self._tissue_dummies = pd.get_dummies(self.tissue_factor)
        columns = self._tissue_dummies.columns
        columns = sorted(columns, key=lambda s: s.lower())

        columns = ['C(tissue)[T.' + x + ']' for x in columns]
        self._tissue_dummies.columns = columns
        N = len(self._tissue_dummies)
        self._tissue_dummies['C(msi)[T.1]'] = [1]*N
        self._tissue_dummies['feature'] = [1] * N
        self._tissue_dummies.insert(0, 'Intercept', [1] * N)

        # drop first feature in the tissues used a reference in the regression
        # # should be if there are at least 2 tissues ?
        tissues = [x for x in self._tissue_dummies.columns if 'tissue' in x]
        self._tissue_dummies.drop(tissues[0], axis=1, inplace=True)

        # reset the buffer.
        self.individual_anova = {}

    def _get_cosmics(self):
        return self.ic50.cosmicIds
    def _set_cosmics(self, cosmics):
        self.ic50.cosmicIds = cosmics
        self.features.cosmicIds = cosmics
        self._init()
        self.individual_anova = {}
    cosmicIds = property(_get_cosmics, _set_cosmics,
        doc="get/set the cosmic identifiers in the IC50 and feature matrices")

    def _get_drug_names(self):
        return self.ic50.drugIds
    def _set_drug_names(self, drugs):
        self.ic50.drugIds = drugs
        self._init()
        # not need to init this again ? self.individual_anova = {}
    drugIds = property(_get_drug_names, _set_drug_names,
            doc="Get/Set drug identifers")

    def _get_feature_names(self):
        shift = self.features.shift
        return self.features.features[shift:]
    def _set_features_names(self, features):
        self.features.features = features
        self._init()
        self.individual_anova = {}
    feature_names = property(_get_feature_names, _set_features_names,
            doc="Get/Set feature names")

    def _get_analysis_mode(self):
        modes = []
        if self.settings.analysis_type == 'PANCAN':
            modes.append('tissue')

        if self.settings.include_MSI_factor is True:
            modes.append('msi')

        modes.append('feature')
        return modes

    def diagnostics(self):
        """Return dataframe with information about the analysis

        """
        n_drugs = len(self.ic50.drugIds)
        n_features = len(self.features.features) - self.features.shift
        n_combos = n_drugs * n_features
        feasible = 0
        pb = Progress(n_drugs, 1)
        counter = 0
        for drug in self.ic50.drugIds:
            for feature in self.features.features[3:]:
                dd = self._get_one_drug_one_feature_data(drug, feature,
                        diagnostic_only=True)
                if dd.status is True:
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

    def read_drug_decoder(self, filename=None):
        """Read file with the DRUG information

        .. seealso:: :class:`gdsctools.readers.DrugDecoder`
        """
        # Read the DRUG decoder file into a DrugDecoder/Reader instance
        self.drug_decoder = readers.DrugDecoder(filename)

    def drug_annotations(self, df):
        """Populate the drug_name and drug_target field if possible

        :param df: input dataframe as given by e.g., :meth:`anova_one_drug`
        :return df: same as input but with the FDR column populated
        """
        if len(self.drug_decoder.df) == 0:
            print("Nothing done. DrugDecoder file not provided.")

        # aliases
        decoder = self.drug_decoder.df
        drugs = df.DRUG_ID.values


        drug_names = [decoder.ix[x].DRUG_NAME if x in decoder.index else None
                 for x in drugs]
        drug_target = [decoder.ix[x].DRUG_TARGET if x in decoder.index
                else None for x in drugs]

        # this is not clean. It works but could be simpler surely.
        df['DRUG_NAME'] = drug_names
        df['DRUG_TARGET'] =  drug_target
        return df

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
        drug_name = self.drug_decoder.get_name(drug_id)
        drug_target = self.drug_decoder.get_target(drug_id)

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

            # We could use pd.get_dummies but pretty slow
            # instead we create the full matrix in _init() method.
            # One issue is that some columns end up with sum == 0
            # and needs to be dropped.
            df = self._tissue_dummies.ix[odof.masked_tissue.index]
            todrop = df.columns[df.values.sum(axis=0) == 0]
            if len(todrop) > 0: # use if since drop() is slow
                df = df.drop(todrop, axis=1)

            # Here we set other variables with dataframe columns' names as
            # expected by OLS.
            if self.settings.include_media_factor is True:
                raise NotImplementedError

            df['C(msi)[T.1]'] = odof.masked_msi.values
            df['feature'] = odof.masked_features.values

            self.Y = odof.Y
            self.EV = df.values
            # The regression and anova summary are done here
            #
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
            #self.anova_pvalues = self._get_anova_summary(self.data_lm,
            #    output='dict')

            # example of computing null model ?
            # Example of computing pvalues ourself
            """self.samples1 = []
            self.samples2 = []
            self.samples3 = []
            Y = odof.Y.copy()
            pb = Progress(10000,20)
            for i in range(0,10000):
                #pylab.shuffle(Y)
                #data_lm = OLS(Y, df.values).fit()
                data_lm = OLS(Y+0.3*pylab.randn(len(Y)), df.values).fit()
                anova_pvalues = self._get_anova_summary(data_lm,
                    output='dict')
                self.samples1.append(anova_pvalues['msi'])
                self.samples2.append(anova_pvalues['feature'])
                self.samples3.append(anova_pvalues['tissue'])
                pb.animate(i)
            """

        elif self.settings.include_MSI_factor is True:
            #self._mydata = pd.DataFrame({'Y': odof.Y,
            #    'msi':  odof.masked_msi, 'feature': odof.masked_features})
            #self.data_lm = ols('Y ~ C(msi) + feature',
            #    data=self._mydata).fit() #Specify C for Categorical
            df = pd.DataFrame()
            df['C(msi)[T.1]'] = odof.masked_msi.values
            df['feature'] = odof.masked_features.values
            df.insert(0, 'Intercept', [1] * (odof.Npos + odof.Nneg))
            self.data_lm = OLS(odof.Y, df.values).fit()
        else:
            df = pd.DataFrame()
            df['feature'] = odof.masked_features.values
            df.insert(0, 'Intercept', [1] * (odof.Npos + odof.Nneg))
            self.data_lm = OLS(odof.Y, df.values).fit()
            #self._mydata = pd.DataFrame({'Y': odof.Y,
            #    'feature': odof.masked_features})
            #self.data_lm = ols('Y ~ feature',
            #    data=self._mydata).fit() #Specify C for Categorical

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
                'TISSUE_ANOVA_pval': tissue_PVAL,
                'MSI_ANOVA_pval': MSI_PVAL,
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
        alphas = pylab.linspace(0,5, Nalpha)

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
            alphas = pylab.linspace(0,1, Nalpha)

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

        # create the W matrix using tissue and MSI if requested
        # default is that the 3 features are used
        if 'tissue' in modes and 'msi' in modes:
            dof = [Ntissue, 1, 1]
            indices = ['tissue', 'msi', 'feature', 'Residuals']
            # 4 stands for intercept + tissue + msi +feature
            arr = np.zeros((4, Ncolumns))
            arr[1, slice(1, Ntissue)] = 1
            arr[2, Ntissue + 1] = 1
            arr[3, Ntissue + 2] = 1
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
            return {'tissue': F_pvalues[0], 'msi':F_pvalues[1],
                    'feature':F_pvalues[2]}
        elif self.settings.include_MSI_factor is True:
            return {'msi': F_pvalues[0], 'feature':F_pvalues[1]}
        else:
            return {'feature': F_pvalues[0]}

        #return anova

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
        print("""
        Analysis of Variance Table

        Response: Y
        Df  Sum Sq Mean Sq F value  Pr(>F)
        TISSUEpattern  26  352.35 13.5517  9.2685 < 2e-16 ***
        MSIpattern      1    5.31  5.3094  3.6313 0.05705 .
        FEATpattern     1    3.19  3.1861  2.1791 0.14028
        Residuals     817 1194.55  1.4621
        """)

    #98% of time in  method anova_one_drug_one_feature
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

        df = ANOVAResults().astype(df)
        if len(df) == 0:
            return df

        if len(self.drug_decoder) > 0:
            df = self.drug_annotations(df)
        # TODO: drop rows where FEATURE_ANOVA_PVAL is None
        if output != 'object':
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
        :return: an :class:`ANOVAResults` instance with the dataframe
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

        if len(self.individual_anova) == 0:
            return ANOVAResults()

        df = pd.concat(self.individual_anova, ignore_index=True)

        if len(df) == 0:
            return df
        # sort all data by ANOVA p-values
        try:
            df.sort_values('FEATURE_ANOVA_pval', inplace=True)
        except:
            df.sort('FEATURE_ANOVA_pval', inplace=True)

        # all ANOVA have been computed individually for each drug and each
        # feature. Now, we need to compute the multiple testing corrections
        df = self.add_pvalues_correction(df)

        # insert a unique identifier as first column
        df.insert(0, 'ASSOC_ID', range(1, len(df) + 1))

        # order the column names as defined in the __init__ method
        df = df[self.column_names]
        df.reset_index(inplace=True, drop=True)

        results = ANOVAResults()
        results.df = df
        results.settings = ANOVASettings(**self.settings)
        return results

    def add_pvalues_correction(self, df, colname='FEATURE_ANOVA_pval'):
        """Add the corrected pvalues column in a dataframe based on pvalues

        The default method (FDR correction) is stored in
        :attr:`settings.pval_correction_method` and can be changed to other
        methods (e.g., *qvalue*)

        .. seealso:: :meth:`anova_all`,
            :class:`~gdsctools.stats.MultipleTesting`
        """
        if len(df) == 0:
            return

        # extract pvalues
        data = df[colname].values

        # set the method and compute new pvalues
        self.multiple_testing.method = self.settings.pval_correction_method
        new_pvalues = self.multiple_testing.get_corrected_pvalues(data)
        new_pvalues *= 100
        # insert new columns.
        colname = 'ANOVA_FEATURE_FDR_%'

        try:
            df.insert(len(df.columns), colname, new_pvalues)
        except:
            # replaces it otherwise
            df[colname] = new_pvalues
        return df

    def reset_buffer(self):
        self.individual_anova = {}

    def __str__(self):
        txt = self.ic50.__str__()
        txt += "\n" + self.features.__str__()
        return txt

    def __repr__(self):
        txt = self.__str__()
        return txt


class HTMLPageMANOVA(ReportMAIN):
    """Creates an HTML page dedicated to significant hits

    ::

        # analyse the data across all drugs and features
        results = gdsc.anova_all()
        # Create a table for the first 10 significant hits (table is sorted
        # by ascending FDR)
        h = HTMLPageMANOVA(gdsc, results.df[0:10])
        # Create the HTML page (pops up by default)
        h.report()

    """
    def __init__(self, gdsc, df):
        """.. rubric:: constructor

        :param : an ANOVA instance.
        :param directory: where to save the file

        The HTML filename is stored in the :attr:`filename`, which can
        be changes (default is manova.html)
        """
        super(HTMLPageMANOVA, self).__init__(filename='manova.html',
                directory=gdsc.settings.directory, 
                template_filename='manova.html')
        #self.template = self.env.get_template('manova.html')

        html = SignificantHits(df, 'all hits').to_html()
        self.jinja['manova'] = html
        self.jinja['analysis_domain'] = gdsc.settings.analysis_type


##############################################################################
#                                                                            #
#                                                                            #
#                    HTML REPORTS RELATED                                    #
#                                                                            #
#                                                                            #
##############################################################################


class Association(ReportMAIN):
    def __init__(self, report, drug=None, feature=None,
            fdr=-1, assoc_id=-1):
        try:
            # here report is expected to be ANOVAReport
            # gdsctools interface from ipython
            self.factory = report.gdsc
        except:
            # here report is actually ANOVA itself
            # pipeline/standalone version
            self.factory = report
        # Does that changes the main settings ?? YES
        self.factory.settings.savefig = True
        self.assoc_id = assoc_id

        self.drug = drug
        self.feature = feature
        self.fdr = fdr
        filename = "{0}____{1}.html".format(self.drug,
                self.feature.replace(" ", "_"))

        super(Association, self).__init__(
                directory=report.settings.directory,
                filename=filename, template_filename='association.html')
        self.jinja['analysis_domain'] = report.settings.analysis_type

    def run(self):
        # to keep . Used in the standalone version
        df = self.factory.anova_one_drug_one_feature(self.drug,
                self.feature, show=True,
                directory=self.directory + os.sep + 'images')
        df['ASSOC_ID'] = self.assoc_id
        df['ANOVA_FEATURE_FDR_%'] = self.fdr
        return df

    def _create_report(self, onweb=True):
        # generated pictures and results
        df = self.run()

        # Create the table and add it
        sign = SignificantHits(df, 'features')
        html_table = sign.to_html(escape=False, header=True, index=False)
        self.jinja['association_table'] = html_table

        # Main boxplot always included
        prefix = 'images/ODOF_all'
        tag = "{0}_{1}____{2}.png".format(prefix, self.drug, self.feature)
        section = '<img alt="association {0}" src="{0}">\n'.format(tag)

        if self.factory.settings.include_MSI_factor:
            prefix = 'ODOF_msi'
            tag = "{0}_{1}____{2}.png".format(prefix, self.drug, self.feature)
            section += '<img alt="association {0}" src="{0}">\n'.format(tag)
        if self.factory.settings.analysis_type == 'PANCAN':
            prefix = 'ODOF_tissue'
            tag = "{0}_{1}____{2}.png".format(prefix, self.drug, self.feature)
            section += '<img alt="association {0}" src="{0}">\n'.format(tag)
        self.jinja['boxplots'] = section


class HTMLOneFeature(ReportMAIN):
    def __init__(self, report, data, subdata, feature):
        self.df = data
        self.subdf = subdata
        self.feature = feature
        self.settings = report.settings

        filename = "{0}.html".format(self.feature)
        super(HTMLOneFeature, self).__init__(
                directory=report.settings.directory,
                filename=filename, template_filename='feature.html')
        self.title = 'Single Feature analysis (%s)' % self.feature

        self.jinja['n_cell_lines'] = report.gdsc.features.df[feature].sum()
        self.jinja['feature_name'] = feature
        self.jinja['analysis_domain'] = report.settings.analysis_type

    def create_pictures(self):
        v = VolcanoANOVA(self.df, settings=self.settings)
        v.volcano_plot_one_feature(self.feature)
        v.savefig_and_js('images/volcano_{}'.format(self.feature))

        # See https://github.com/CancerRxGene/gdsctools/issues/79
        v.current_fig.canvas.mpl_disconnect(v.cid)
        try:
            import mpld3
            mpld3.plugins.clear(v.current_fig)
        except:
            pass

    def _create_report(self, onweb=True):
        self.create_pictures()

        self.jinja['N_hits'] = len(self.subdf)
        if len(self.subdf) > 0:
            sign = SignificantHits(self.subdf, 'feature')
            html = sign.to_html(escape=False, header=True, index=False)
            self.jinja['association_table'] = html

        # image section
        self.jinja['image_filename'] = "images/volcano_{0}".format(self.feature)


class HTMLOneDrug(ReportMAIN):
    def __init__(self, report, data, subdata, drug):
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
        self.drug = drug
        self.settings = report.settings

        filename = "{0}.html".format(self.drug)
        super(HTMLOneDrug, self).__init__(
                directory=report.settings.directory,
                filename=filename, template_filename='drug.html')
        self.title = 'Single Drug analysis (%s)' % self.drug
           
        self.jinja['n_cell_lines'] = len(report.gdsc.ic50.df[drug].dropna())
        self.jinja['drug_id'] = drug
        self.jinja['drug_name'] = report.gdsc.drug_decoder.get_name(drug)
        self.jinja['drug_target'] = report.gdsc.drug_decoder.get_target(drug)
        self.jinja['analysis_domain'] = report.settings.analysis_type

    def create_pictures(self):
        v = VolcanoANOVA(self.df, settings=self.settings)
        v.volcano_plot_one_drug(self.drug)
        v.savefig_and_js('images/volcano_{}'.format(self.drug), size_inches=(10,10))

        # See https://github.com/CancerRxGene/gdsctools/issues/79
        v.current_fig.canvas.mpl_disconnect(v.cid)
        try:
            import mpld3
            mpld3.plugins.clear(v.current_fig)
        except:
            pass

    def _create_report(self, onweb=True):
        self.create_pictures()

        # add the table
        self.jinja['synonyms'] = ''
        self.jinja['brand_name'] = ''
        self.jinja['conc_min'] = '?'
        self.jinja['conc_max'] = '?'

        # Table section
        self.jinja['N_hits'] = len(self.subdf)
        if len(self.subdf)>0:
            sign = SignificantHits(self.subdf, 'drugs')
            html = sign.to_html(escape=False, header=True, index=False)
            self.jinja['association_table'] = html

        # image section
        self.jinja['image_filename'] = "images/volcano_{0}".format(self.drug)


class HTMLPageMain(ReportMAIN):
    def __init__(self, report, filename='index.html'):
        super(HTMLPageMain, self).__init__(
                directory=report.settings.directory,
                filename=filename)
        self.report = report
        self.settings = report.settings
        self.jinja['analysis_domain'] = report.settings.analysis_type
        self.jinja['settings_table'] = self.settings.to_html()

    def _create_report(self, onweb=True):
        # A summary table
        diag = self.report.diagnostics()
        table = HTMLTable(diag, 'summary')
        txt = ''
        for index, row in diag.iterrows():
            if len(row.text) == 0 and len(row.value) == 0:
                txt += '----<br/>'
            else:
                txt += row.text + ": " +  str(row.value) + "<br/>"
        self.jinja['summary'] = txt

        print('Creating volcano plots')
        # this can be pretty slow. so keep only 1000 most relevant
        # values and 1000 random ones to get an idea of the distribution
        v = VolcanoANOVA(self.report.df, settings=self.settings)
        v.selector(v.df, 1500, 1500, inplace=True)
        v.volcano_plot_all()
        v.savefig_and_js("volcano_all_js")

        self.jinja['volcano'] = """
            <h3></h3>
            <a href="volcano_all_js.html">
                <img alt="volcano plot for all associations"
                    src="volcano_all_js.png">
            </a>
            <br/>
            <p> A javascript version is available 
                <a href="volcano_all_js.html">here</a> (
                or click on the image).</p>
        """

        # MANOVA link
        N = len(self.report.get_significant_set())
        self.jinja['manova'] = """
        There were %(N)s significant associations found.
        All significant associations have been gatherered
        in the following link: <br/><a href="manova.html">manova results</a>.
        """ % {'N': N}

        # feature summary
        df_features = self.report.feature_summary("feature_summary.png")
        filename = 'OUTPUT' + os.sep + 'features_summary.csv'
        df_features.to_csv(self.directory + os.sep + filename, sep=',')

        # drug summary
        not_tested = [x for x in self.report.gdsc.drugIds if x not in 
                self.report.df.DRUG_ID.unique()]
        if len(not_tested) > 0:
            not_tested = """Those drugs have not been analysed due to 
            lack of valid data points: """ + ", ".join(not_tested)
        else:
            not_tested = ""
        self.jinja['drug_not_tested'] = not_tested

        df_drugs = self.report.drug_summary(filename="drug_summary.png")
        get_name = self.report.gdsc.drug_decoder.get_name
        if len(self.report.gdsc.drug_decoder.df) > 0:
            df_drugs.index = [x + "-" + get_name(x) for x in df_drugs.index]
        filename = 'OUTPUT' + os.sep + 'drugs_summary.csv'
        df_drugs.to_csv(self.directory + os.sep + filename, sep=',')

        # --------------------------- Create table with links to all drugs
        groups = self.report.df.groupby('DRUG_ID')
        try:
            df = groups.mean()['ANOVA_FEATURE_FDR_%'].sort_values()
        except:
            # note double brackets for pythonn3.3
            df = groups.mean()[['ANOVA_FEATURE_FDR_%']].sort()
        df = df.reset_index() # get back the Drug id in the dframe columns

        # let us add also the drug name
        df = self.report.gdsc.drug_annotations(df)

        # let us also add number of associations computed
        counts = [len(groups.groups[k]) for k in df.DRUG_ID]
        df['Number of associations computed'] = counts
        groups = self.report.get_significant_set().groupby('DRUG_ID').groups
        count = []
        for drug in df['DRUG_ID'].values:
            if drug in groups.keys():
                count.append(len(groups[drug]))
            else:
                count.append(0)
        df['hits'] = count

        # add another set of drug_id but sorted in alpha numerical order
        table = HTMLTable(df, 'drugs')
        table.add_href('DRUG_ID')
        table.df.columns = [x.replace('ANOVA_FEATURE_FDR',
            'mean ANOVA FEATURE FDR') for x in table.df.columns]
        table.add_bgcolor('hits', mode='max',
                cmap=cmap_builder('white', 'orange', 'red'))

        self.jinja['drug_table'] = table.to_html(escape=False,
                header=True, index=False)

        # ---------------------- Create full table with links to all features
        df = pd.DataFrame({'FEATURE': self.report.df['FEATURE'].unique()})
        try:
            df.sort_values(by='FEATURE', inplace=True)
        except:
            df.sort('FEATURE', inplace=True)

        groups = self.report.get_significant_set().groupby('FEATURE').groups

        count = []
        for feature in df['FEATURE'].values:
            if feature in groups.keys():
                count.append(len(groups[feature]))
            else:
                count.append(0)
        df['hits'] = count

        table = HTMLTable(df, 'features')
        table.sort('hits', ascending=False)
        table.add_href('FEATURE')
        table.add_bgcolor('hits', mode='max',
                cmap=cmap_builder('white', 'orange', 'red'))
        self.jinja['feature_table'] = table.to_html(escape=False,
                header=True, index=False)

        # -------------------------------------- COSMIC table for completeness
        colnames = self.report.gdsc.features._special_names
        df = self.report.gdsc.features.df[colnames]

        df = df.reset_index()
        table = HTMLTable(df)
        url = "http://cancer.sanger.ac.uk/cell_lines/sample/overview?id="
        table.add_href('COSMIC_ID', url=url, newtab=True)
        self.jinja['cosmic_table'] = table.to_html()

        # -------------------------------------- settings and INPUT files
        input_dir = self.directory + os.sep + 'INPUT'
        filename = self.report.gdsc.ic50._filename
        html = ''
        if filename is not None:
            shutil.copy(filename, input_dir)
            filename = os.path.basename(filename)
            ic50_filename = filename
            filename = os.sep.join(['INPUT', filename])
            self.jinja['ic50_file'] = filename
        else:
            ic50_filename = 'unknown'

        # the genomic features, which may be the default version or without
        # location
        gf_filename = self.report.gdsc.features._filename
        if gf_filename is None: 
            gf_filename = os.sep.join([input_dir, 'genomic_features.csv'])
            self.report.gdsc.features.to_csv(gf_filename)
            html = """Saved <a href="INPUT/genomic_features.csv">Genomic 
                      Features</a> file<br/> (possibly the default 
                      version)."""
            self.jinja['gf_file'] = html
        else:
            shutil.copy(gf_filename, input_dir)
            gf_filename = os.path.basename(gf_filename)
            txt = """Get <a href="INPUT/%s">Genomic Features</a> 
                     file.<br/>"""
            self.jinja['gf_file'] = txt % gf_filename

        # the drug decode file
        filename = self.report.gdsc.drug_decoder._filename
        if filename is not None:
            output_filename = 'DRUG_DECODE.csv'
            shutil.copy(filename, input_dir + os.sep + output_filename)
            filename = os.path.basename(filename)
            html = 'Get <a href="INPUT/%s">Drug DECODE file</a>' % output_filename
        else:
            html = 'No Drug DECODE file was provided'
        self.jinja['drug_decode'] = html

        # Save settings as json file
        filename = os.sep.join([input_dir, 'settings.json'])
        self.settings.to_json(filename)
        filename = os.path.basename(filename)
        self.jinja['settings'] = \
                """Get the settings as a <a href="INPUT/%s">
                json file</a>.""" % filename

        # Save all Results dataframe
        filename = os.sep.join([self.settings.directory, 'OUTPUT', 
            'results.csv'])
        ANOVAResults(self.report.df).to_csv(filename)

        code = """from gdsctools import *
import os

def getfile(filename, where='../INPUT'):
    return os.sep.join([where, filename])

# reback the IC50 and genomic features matrices
gdsc = ANOVA(getfile('%(ic50)s'), getfile('%(gf_filename)s'))
gdsc.settings.from_json(getfile('settings.json'))

# Analyse the data
results = gdsc.anova_all()

# Create the HTML report
r = ANOVAReport(gdsc, results)
r.create_html_pages(onweb=False)"""
        code = code % {
                'ic50': ic50_filename, 
                'gf_filename': gf_filename}

        filename = os.sep.join([self.settings.directory, 'code','rerun.py'])
        fh = open(filename, 'w')
        fh.write(code)
        fh.close()



class SignificantHits(object):
    def __init__(self, df, name):
        self.df = df.copy() # to not alter original version
        self.cmap_clip = cmap_builder('#ffffff', '#0070FF')
        self.cmap_absmax = cmap_builder('green', 'white', 'red')
        self.name = name
        self.clip_threshold = 2

        columns = [u'ASSOC_ID', 'FEATURE',
            'DRUG_ID', u'DRUG_NAME', 'DRUG_TARGET',
            'N_FEATURE_neg', 'N_FEATURE_pos',
            'FEATUREpos_logIC50_MEAN',
            'FEATUREneg_logIC50_MEAN',
            'FEATURE_deltaMEAN_IC50',
            'FEATURE_IC50_effect_size',
            'FEATUREneg_Glass_delta',
            'FEATUREpos_Glass_delta',
            'FEATURE_ANOVA_pval',
            'TISSUE_ANOVA_pval',
            'MSI_ANOVA_pval',
            'ANOVA_FEATURE_FDR_%']
        self.df = self.df[columns]

    def to_html(self, escape=False, header=True, index=False):
        # If there is a value below 0.01, a scientific notation is
        # use. we prefer to use 2 digits and write <0.01
        colname = 'ANOVA_FEATURE_FDR_%'
        self.df.loc[self.df[colname] <0.01, colname] = '<0.01'
        html = HTMLTable(self.df, self.name)
        # Those columns should be links
        for this in ['FEATURE', 'DRUG_ID', 'ASSOC_ID']:
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
