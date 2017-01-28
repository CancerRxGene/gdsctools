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
import warnings

from gdsctools.report import HTMLTable, ReportMain
from gdsctools.tools import Savefig
from gdsctools.volcano import VolcanoANOVAJS
from gdsctools.boxplots import BoxPlotsJS
from gdsctools.settings import ANOVASettings
from gdsctools.anova_results import ANOVAResults
from gdsctools.readers import DrugDecode

import pandas as pd
import pylab
import numpy as np

from easydev import Progress
import easydev
from colormap import cmap_builder


__all__ = ['ANOVAReport']


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


    .. rubric:: Significant association

    An association is significant if

      - The field *ANOVA_FEATURE_FDR* must be < FDR_threshold
      - The field *ANOVA_FEATURE_pval* must be < pvalue_threshold

    It is then labelled **sensible** if *FEATURE_delta_MEAN_IC50* is 
    below 0, otherwise it is **resistant**.


    """
    def __init__(self, gdsc, results=None, sep="\t", drug_decode=None, 
            verbose=True):
        """.. rubric:: Constructor

        :param gdsc: the instance with which you created the results to report
        :param results: the results returned by :meth:`ANOVA.anova_all`. If
            not provided, the ANOVA is run on the fly.

        """
        self.verbose = verbose
        self.figtools = Savefig(verbose=False)
        self.gdsc = gdsc

        if results is None:
            results = gdsc.anova_all()
        self.df = ANOVAResults(results).df # this does a copy and sanity check


        # Make sure the DRUG are integers
        self.df.DRUG_ID = self.df.DRUG_ID.astype(int)


        self.settings = ANOVASettings()
        for k, v in gdsc.settings.items():
            self.settings[k] = v

        self._colname_drug_id = 'DRUG_ID'
        self.varname_pval = 'ANOVA_FEATURE_pval'
        self.varname_qval = 'ANOVA_FEATURE_FDR'

        # maybe there was not drug_decode in the gdsc parameter,
        # so a user may have provide a file, in which case, we need
        # to update the content of the dur_decoder.
        if len(gdsc.drug_decode) == 0 and drug_decode is None:
            warnings.warn("No drug name or target will be populated."
                          "You may want to provide a DRUG_DECODE file.")
            self.drug_decode = DrugDecode()
        elif drug_decode is not None:
            # Read a file
            self.drug_decode = DrugDecode(drug_decode)
        else:
            # Copy from gdsc instance
            self.drug_decode = DrugDecode(gdsc.drug_decode)
        self.df = self.drug_decode.drug_annotations(self.df)
        #if sum(self.df == np.inf).sum()>0:
        #    print("WARNING: infinite values were found in your results... Set to zero")
        try:
            self.df = self.df.replace({np.inf:0, -np.inf:0})
        except:
            pass

        # create some data
        self._set_sensible_df()

        self.company = None

        # just to create the directory
        # ReportMain(directory=self.settings.directory, verbose=self.verbose)

    def _get_ndrugs(self):
        return len(self.df[self._colname_drug_id].unique())
    n_drugs = property(_get_ndrugs, doc="return number of drugs")

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
        self._set_sensible_df()

        df = pd.DataFrame({'text': [], 'value': []})

        n_features = len(self.gdsc.features.df.columns)
        n_features -= self.gdsc.features.shift
        n_drugs = len(self.df[self._colname_drug_id].unique())

        N = float(n_drugs * n_features)

        if N == 0:
            ratio = 0
        else:
            ratio = float(self.n_tests)/(N) * 100

        try:
            ratio = easydev.precision(ratio, digit=2)
        except:
            # Fixme: this is a hack for the infinite values but should not
            # happen...
            ratio = 0

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
        df = self._df_append(df, [msg, self.settings.pvalue_threshold])

        msg = "FDR significance threshold"
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
        mask1 = self.df['ANOVA_FEATURE_FDR'] < self.settings.FDR_threshold
        mask2 = self.df['ANOVA_FEATURE_pval'] < self.settings.pvalue_threshold
        mask3 = self.df['FEATURE_delta_MEAN_IC50'] < 0
        self.sensible_df = self.df[logand(logand(mask1, mask2), mask3)]

        # select resistant data set
        mask3 = self.df['FEATURE_delta_MEAN_IC50'] >= 0
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
            df_count.sort_values(by=['total', 'name'], 
                    ascending=[False, True], inplace=True)
        except:
            df_count.sort(columns=['total', 'name'], 
                    ascending=[False, True], inplace=True)

        df_count.drop('name', axis=1, inplace=True)
        return df_count

    def get_drug_summary_data(self):
        """Return dataframe with drug summary"""
        # get sensible and resistant sub dataframes
        self._set_sensible_df()

        # group by drug
        colname = self._colname_drug_id
        df_count_sensible = self.sensible_df.groupby(colname).count()
        df_count_resistant = self.resistant_df.groupby(colname).count()

        df_count = self._get_data(df_count_sensible, df_count_resistant)
        return df_count

    def drug_summary(self,  top=50, fontsize=15, filename=None):
        """Return dataframe with significant drugs and plot figure

        :param fontsize:
        :param top: max number of significant associations to show
        :param filename: if provided, save the file in the directory
        """

        df_count = self.get_drug_summary_data()

        if len(df_count):
            self._plot(df_count, 'drug', top)
            fig = pylab.gcf()
            self.figtools.directory = self.settings.directory
            self.figtools.savefig(filename, size_inches=(12, 14),
                    bbox_inches='tight')
        return df_count

    def get_feature_summary_data(self):
        """Return dataframe with feature summary"""
        # get sensible and resistant sub dataframes
        self._set_sensible_df()
        df_count_sensible = self.sensible_df.groupby('FEATURE').count()
        df_count_resistant = self.resistant_df.groupby('FEATURE').count()
        df_count = self._get_data(df_count_sensible, df_count_resistant)
        return df_count

    def feature_summary(self, filename=None, top=50, fontsize=15):
        """Return dataframe with significant features and plot figure

        :param fontsize:
        :param top: max number of significant associations to show
        :param filename: if provided, save the file in the directory
        """
        df_count = self.get_feature_summary_data()

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

        df = df_count.iloc[0:top][[u'sens assoc', u'res assoc']]
        labels = list(df.index)
        # add drug name
        if len(self.drug_decode) > 0:
            for i, label in enumerate(labels):
                if title_tag == 'drug':
                    name = self.drug_decode.get_name(label)
                    if name is not None:
                        labels[i] = "{}-{}".format(labels[i], name)
                else:
                    pass

        labels = [str(x).replace('_', ' ') for x in labels]
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

        xticks = ax.get_xticks()
        ax.set_xticklabels(
                [int(x) if divmod(x,1)[1] == 0 else "" for x in xticks])


        pylab.grid()
        pylab.title(r"Top %s %s(s) most frequently " % (top, title_tag) + \
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

        MC1 = 1
        MC2 = 2
        mask2 = self.df['FEATURE_pos_logIC50_MEAN'] < MC1
        mask3 = self.df['FEATURE_pos_logIC50_MEAN'] < MC2
        mask4 = self.df['FEATURE_neg_logIC50_MEAN'] < MC1
        mask5 = self.df['FEATURE_neg_logIC50_MEAN'] < MC2
        maskMC = mask2 + mask3 + mask4 + mask5

        for fdr in fdrs:
            # significant hits
            res = self.df['ANOVA_FEATURE_FDR'] < fdr
            significants.append(res.sum())

            # meaningful hits
            indices = np.logical_and(self.df['ANOVA_FEATURE_FDR'] < fdr,
                    maskMC)
            significant_meaningful.append(indices.sum())

            # meaningful strong hits
            mask1 = self.df.ix[indices]['FEATURE_pos_Glass_delta'] >= 1
            mask2 = self.df.ix[indices]['FEATURE_neg_Glass_delta'] >= 1
            strong_hits.append(np.logical_or(mask1, mask2).sum())

            # meaningful full strong hits
            mask1 = self.df.ix[indices]['FEATURE_pos_Glass_delta'] >= 1
            mask2 = self.df.ix[indices]['FEATURE_neg_Glass_delta'] >= 1
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
        if self.verbose:
            print("Creating individual HTML pages for each significant association")
        df = self.get_significant_set()

        drugs = df['DRUG_ID'].values
        features = df['FEATURE'].values
        assocs = df['ASSOC_ID'].values
        fdrs = df['ANOVA_FEATURE_FDR'].values

        N = len(df)
        pb = Progress(N)

        html = Association(self, drug='dummy', feature='dummy', fdr='dummy',
                company=self.company)

        for i in range(N):
            html.drug = drugs[i]
            html.feature = features[i]
            if str(assocs[i]).startswith("a"):
                html._filename = str(assocs[i]) + '.html'
            else:
                html._filename = "a" + str(assocs[i]) + '.html'
            html.fdr = fdrs[i]
            html.assoc_id = assocs[i]
            #html._init_report() # since we have one shared instance
            html.create_report(onweb=False)
            if self.settings.animate:
                pb.animate(i+1)
        if self.settings.animate: print("\n")

    def create_html_features(self):
        """Create an HTML page for each significant feature"""
        df = self.get_significant_set()
        groups = df.groupby('FEATURE')
        if self.verbose:
            print("Creating individual HTML pages for each feature")
        N = len(groups.indices.keys())
        pb = Progress(N)
        for i, feature in enumerate(groups.indices.keys()):
            # get the indices and therefore subgroup
            subdf = groups.get_group(feature)
            html = HTMLOneFeature(self, self.df, subdf, feature)
            html.create_report(onweb=False)
            if self.settings.animate:
                pb.animate(i+1)
        if self.settings.animate: print("\n")

    def create_html_drugs(self):
        """Create an HTML page for each drug"""
        # group by drugs
        all_drugs = list(self.df['DRUG_ID'].unique())

        df = self.get_significant_set()
        groups = df.groupby('DRUG_ID')
        if self.verbose:
            print("Creating individual HTML pages for each drug")
        N = len(groups.indices.keys())
        N = len(all_drugs)
        pb = Progress(N)
        for i, drug in enumerate(all_drugs):
            # enumerate(groups.indices.keys()):
            # get the indices and therefore subgroup
            if drug in groups.groups.keys():
                subdf = groups.get_group(drug)
            else:
                subdf = {}

            html = HTMLOneDrug(self, self.df, subdf, drug)
            html.create_report(onweb=False)
            if self.settings.animate:
                pb.animate(i+1)
        if self.settings.animate: print("\n")

    def create_html_main(self, onweb=False):
        """Create HTML main document (summary)"""
        self._set_sensible_df()

        if self.verbose:
            print("Creating main HTML page in directory %s" %
                (self.settings.directory))
        ReportMain(directory=self.settings.directory, verbose=self.verbose)

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
        page = HTMLPageMANOVA(self, df, self.company)
        page.create_report(onweb)

    def create_html_pages(self, onweb=True):
        """Create all HTML pages"""
        self.create_html_main(onweb=onweb)
        self.create_html_manova(onweb=False)
        self.create_html_drugs()
        self.create_html_features()
        self.create_html_associations()

    def onweb(self):
        from easydev import onweb
        onweb(self.settings.directory + os.sep + 'index.html')


class HTMLPageMANOVA(ReportMain):
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
    def __init__(self, report, df, company):
        """.. rubric:: constructor

        :param : an ANOVA instance.
        :param directory: where to save the file

        The HTML filename is stored in the :attr:`filename`, which can
        be changes (default is manova.html)
        """
        super(HTMLPageMANOVA, self).__init__(filename='manova.html',
                directory=report.settings.directory+os.sep+"associations",
                template_filename='manova.html', init_report=False)

        html = ANOVAResults(df).get_html_table(collapse_table=False)

        self.jinja['manova'] = html
        self.jinja['analysis_domain'] = report.settings.analysis_type
        self.jinja['resource_path'] = ".."
        self.jinja["collaborator"] = company

    def _create_report(self):
        pass # used to avoid warning

##############################################################################
#                                                                            #
#                                                                            #
#                    HTML REPORTS RELATED                                    #
#                                                                            #
#                                                                            #
##############################################################################


class Association(ReportMain):
    def __init__(self, report, drug=None, feature=None,
            fdr=-1, assoc_id=-1, company=None):

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
                directory=report.settings.directory + os.sep + "associations",
                filename=filename, template_filename='association.html',
                init_report=False)
        self.jinja['analysis_domain'] = report.settings.analysis_type
        self.jinja['resource_path'] = ".."
        self.jinja["collaborator"] = company

    def run(self):
        # to keep . Used in the standalone version
        df = self.factory.anova_one_drug_one_feature(self.drug,
                self.feature, show=False,
                directory=self.directory)
        df['ASSOC_ID'] = self.assoc_id
        df['ANOVA_FEATURE_FDR'] = self.fdr
        return df

    def _create_report(self, onweb=True):
        # generated pictures and results
        df = self.run()
        odof = self.factory._get_one_drug_one_feature_data(self.drug, 
                self.feature)

        # Create the table and add it
        sign = ANOVAResults(df)

        html_table = sign.get_html_table(escape=False, header=True,
                index=False)

        self.jinja['association_table'] = html_table

        # Javascript version
        bx = BoxPlotsJS(odof)
        self.jinja["boxplot_all_jsdata"] = bx.get_html_association()

        section = ""
        if self.factory.settings.include_MSI_factor:
            self.jinja["boxplot_msi_jsdata"] = bx.get_html_msi()
        if self.factory.settings.include_media_factor:
            self.jinja["boxplot_media_jsdata"] = bx.get_html_media()
        if self.factory.settings.analysis_type == 'PANCAN':
            self.jinja["boxplot_tissue_jsdata"] = bx.get_html_tissue()

        self.jinja['boxplots'] = section

class HTMLOneFeature(ReportMain):
    def __init__(self, report, data, subdata, feature):
        self.df = data
        self.subdf = subdata
        self.feature = feature
        self.settings = report.settings

        self.n_hits = 0
        if len(report.resistant_df):
            self.n_hits += sum([True for x in report.resistant_df.FEATURE
                if x == feature])
        if len(report.sensible_df):
            self.n_hits += sum([True for x in report.sensible_df.FEATURE
                if x == feature])

        filename = "{0}.html".format(self.feature)
        super(HTMLOneFeature, self).__init__(
                directory=report.settings.directory + os.sep + "associations",
                filename=filename, template_filename='feature.html',
                init_report=False)
        self.title = 'Single Feature analysis (%s)' % self.feature

        self.jinja['n_cell_lines'] = report.gdsc.features.df[feature].sum()
        self.jinja['feature_name'] = feature
        self.jinja['analysis_domaiin'] = report.settings.analysis_type
        self.jinja['n_hits'] = self.n_hits
        self.jinja['resource_path'] = ".."
        self.jinja["collaborator"] = report.company

    def create_pictures(self):
        v = VolcanoANOVAJS(self.df, settings=self.settings)
        html = v.render_feature(self.feature)
        return html

    def _create_report(self, onweb=True):
        self.jinja['N_hits'] = len(self.subdf)
        if len(self.subdf) > 0:
            sign = ANOVAResults(self.subdf)
            html = sign.get_html_table(escape=False, 
                    header=True, index=False)
            self.jinja['association_table'] = html
        self.jinja['volcano_jsdata'] =  self.create_pictures()


class HTMLOneDrug(ReportMain):
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
        self.drug = int(drug)
        self.settings = report.settings

        filename = "drug_{0}.html".format(self.drug)
        super(HTMLOneDrug, self).__init__(
                directory=report.settings.directory + os.sep + "associations",
                filename=filename, template_filename='drug.html',
                init_report=False)
        self.title = 'Single Drug analysis (%s)' % self.drug

        self.nhits = 0
        if len(report.resistant_df):
            self.nhits += sum([True for x in report.resistant_df.DRUG_ID
                if x == drug])
        if len(report.sensible_df):
            self.nhits += sum([True for x in report.sensible_df.DRUG_ID
                if x == drug])

        self.jinja['n_cell_lines'] = len(report.gdsc.ic50.df[self.drug].dropna())
        self.jinja['drug_id'] = self.drug


        self.jinja['drug_name'] = report.drug_decode.get_name(self.drug)
        self.jinja['drug_target'] = report.drug_decode.get_target(self.drug)
        #self.jinja['drug_synonyms'] = report.drug_decode._get_row(self.drug, 'SYNONYMS')
        self.jinja['drug_owner'] = report.drug_decode._get_row(self.drug, 'OWNED_BY')

        self.jinja['analysis_domain'] = report.settings.analysis_type
        self.jinja['n_hits'] = self.nhits
        self.jinja['resource_path'] = ".."
        self.jinja["collaborator"] = report.company

    def create_pictures(self):
        v = VolcanoANOVAJS(self.df, settings=self.settings)
        html = v.render_drug(self.drug)
        return html

    def _create_report(self, onweb=True):
        #self.create_pictures()

        # add the table
        self.jinja['synonyms'] = ''
        self.jinja['brand_name'] = ''
        self.jinja['conc_min'] = '?'
        self.jinja['conc_max'] = '?'

        # Table section
        self.jinja['N_hits'] = len(self.subdf)
        if len(self.subdf) > 0:
            sign = ANOVAResults(self.subdf)
            html = sign.get_html_table(escape=False, 
                    header=True, index=False)
            self.jinja['association_table'] = html

        # image section
        self.jinja['volcano_jsdata'] =  self.create_pictures()


class HTMLPageMain(ReportMain):
    def __init__(self, report, filename='index.html'):
        super(HTMLPageMain, self).__init__(
                directory=report.settings.directory,
                filename=filename)
        self.report = report
        self.settings = report.settings
        self.jinja['analysis_domain'] = report.settings.analysis_type
        self.jinja['settings_table'] = self.settings.to_html()
        self.jinja["collaborator"] = report.company

    def _create_report(self, onweb=True):
        self.add_summary()
        if len(self.report.df):
            self.add_volcano()
            self.add_manova()
            self.add_features()
        else:
            self.jinja["__section_skip_volcano"] = True
            self.jinja["__section_skip_manova"] = True
            self.jinja["__section_skip_features"] = True
            self.jinja["__section_skip_feature_summary"] = True
            self.jinja["__section_skip_drug_summary"] = True
            self.jinja["__section_skip_drug_associations"] = True
            self.jinja["__section_skip_feature_associations"] = True
        self.add_cosmic()
        self.add_settings()

    def add_summary(self):
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

    def add_volcano(self):
        """
        As of version 0.13, we will use Javascript directly. 

        # this can be pretty slow. so keep only 1000 most relevant
        # values and 1000 random ones to get an idea of the distribution
        v = VolcanoANOVA(self.report.df, settings=self.settings)
        #v.selector(v.df, 1500, 1500, inplace=True)
        v.volcano_plot_all()
        v.savefig("volcano_all_js")

        self.jinja['volcano'] = '''
            <h3></h3>
            <a href="volcano_all_js.html">
                <img alt="volcano plot for all associations"
                    src="volcano_all_js.png">
            </a>
            <br/>
            <p> A javascript version is available
                <a href="volcano_all_js.html">here</a> (
                or click on the image).</p>
        '''
        """
        v = VolcanoANOVAJS(self.report.df, settings=self.settings)
        self.jinja['volcano_jsdata'] = v.render_all()

    def add_manova(self):
        # MANOVA link
        N = len(self.report.get_significant_set())
        self.jinja['manova'] = str(N)

    def add_features(self):
        
        # feature summary
        df_features = self.report.feature_summary("feature_summary.png")
        filename = 'OUTPUT' + os.sep + 'features_summary.csv'
        df_features.to_csv(self.directory + os.sep + filename, sep=',')

        not_tested = ""
        self.jinja['drug_not_tested'] = not_tested

        df_drugs = self.report.drug_summary(filename="drug_summary.png")
        get_name = self.report.drug_decode.get_name
        if len(self.report.drug_decode.df) > 0:
            df_drugs.index = ["{}-{}".format(x, get_name(x)) for x in df_drugs.index]
        filename = 'OUTPUT' + os.sep + 'drugs_summary.csv'
        df_drugs.to_csv(self.directory + os.sep + filename, sep=',')

        if len(self.report.df) == 0:
            return

        # --------------------------- Create table with links to all drugs
        groups = self.report.df.groupby('DRUG_ID')
        try:
            df = groups.mean()['ANOVA_FEATURE_FDR'].sort_values()
        except:
            # note double brackets for pythonn3.3
            df = groups.mean()[['ANOVA_FEATURE_FDR']].sort()

        df = df.reset_index() # get back the Drug id in the dframe columns
        # let us add also the drug name
        df = self.report.drug_decode.drug_annotations(df)

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
        table.add_href('DRUG_ID', url="associations/drug_", suffix=".html")
        table.df.columns = [x.replace('ANOVA_FEATURE_FDR',
            'mean FEATURE ANOVA FDR') for x in table.df.columns]
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
        table.add_href('FEATURE', url="associations/", suffix=".html")
        table.add_bgcolor('hits', mode='max',
                cmap=cmap_builder('white', 'orange', 'red'))
        self.jinja['feature_table'] = table.to_html(escape=False,
                header=True, index=False)

    def add_cosmic(self):
        # -------------------------------------- COSMIC table for completeness
        colnames = self.report.gdsc.features._special_names
        df = self.report.gdsc.features.df[colnames]

        # TODO
        # add other columns if possible e.g., GDSC1, GDSC2, TCGA

        df = df.reset_index()
        table = HTMLTable(df)
        url = "http://cancer.sanger.ac.uk/cell_lines/sample/overview?id="
        table.add_href('COSMIC_ID', url=url, newtab=True)
        self.jinja['cosmic_table'] = table.to_html()

    def add_settings(self):
        # -------------------------------------- settings and INPUT files
        input_dir = self.directory + os.sep + 'INPUT'
        filename = 'ANOVA_input.csv'
        filename = os.sep.join([input_dir, filename])
        self.report.gdsc.ic50.to_csv(filename)
        filename = os.sep.join(['INPUT', 'ANOVA_input.csv'])
        self.jinja['ic50_file'] = filename

        # the genomic features, which may be the default version
        # one provided by the user. It may have been changed
        gf_filename = os.sep.join([input_dir, 'genomic_features.csv'])
        self.report.gdsc.features.to_csv(gf_filename)
        html = """Saved <a href="INPUT/genomic_features.csv">Genomic
                  Features</a> file<br/> (possibly the default
                  version)."""
        self.jinja['gf_file'] = html

        # Always save DRUG_DECODE file even if empty
        # It may be be interpreted in other pipeline or for reproducibility
        output_filename = input_dir + os.sep + 'DRUG_DECODE.csv'
        self.report.drug_decode.to_csv(output_filename)
        html = 'Get <a href="INPUT/DRUG_DECODE.csv">Drug DECODE file</a>'
        if len(self.report.drug_decode) == 0:
            html += 'Note that DRUG_DECODE file was not provided (empty?).'
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
gdsc = ANOVA(getfile('%(ic50)s'), getfile('%(gf_filename)s'),
        getfile('DRUG_DECODE.csv'))
gdsc.settings.from_json(getfile('settings.json'))
gdsc.init()

# Analyse the data
results = gdsc.anova_all()

# Create the HTML report
r = ANOVAReport(gdsc, results)
r.create_html_pages(onweb=False)"""
        code = code % {
                'ic50': 'ANOVA_input.csv',
                'gf_filename': 'genomic_features.csv'}

        filename = os.sep.join([self.settings.directory, 'code','rerun.py'])
        fh = open(filename, 'w')
        fh.write(code)
        fh.close()

