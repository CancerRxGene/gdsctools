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
"""ANOVAResults data structure to store the output of the ANOVA analysis"""
from collections import OrderedDict

import pandas as pd
import numpy as np
import pylab

from colormap import cmap_builder

from gdsctools import readers
from gdsctools.volcano import VolcanoANOVA
from gdsctools.report import HTMLTable

__all__ = ['ANOVAResults']


class ANOVAResults(object):
    """Class to handle results of the ANOVA analysis

    The :class:`ANOVA` class and in particular its method
    :meth:`~gdsctools.anova.ANOVA.anova_all` returns the results of
    the ANOVA analysis for each drug and genomic feature. The results
    are stored in a data structure defined in this class, which is 
    just a dataframe stored in :attr:`df` attribute with the following
    header:

    =========================== ===============================================
         Column name                                  Description
    =========================== ===============================================
    ASSOC_ID                    Alphanumeric identifier of the interaction
    FEATURE                     The CFE involved in the interaction, it
                                can be a mutated cancer driver gene
                                (CG) [suffix _mut], an abberrantly
                                fused protein [suffix fusion], a copy
                                number altered chromosomal region (RACS)
                                [prefix gain for amplifications or loss
                                for deletions];
    DRUG_ID                     Numerical id of the drug involved in the
                                interaction;
    DRUG_TARGET                 Putative target of the drug involved in the
                                interaction;
    N_FEATURE_pos               Number of cell lines harbouring the CFE
                                indicated in column E and that have been
                                screened with the drug indicated in columns
                                F and G, therefore have been included in
                                the test;
    N_FEATURE_neg               Number of cell lines not harbouring the
                                CFE indicated in column E and that have been
                                screened with the drug indicated in columns
                                F and G, therefore have been included in the
                                test;
    FEATURE_pos_logIC50_MEAN    Average log IC50 of the population of cell
                                lines accounted in colum i;
    FEATURE_neg_logIC50_MEAN    Average log IC50 of the population of cell
                                lines accounted in colum j;
    FEATURE_delta_MEAN_IC50     Difference between the two average natural
                                log IC50 values in the previous two columns
                                (j - i). A negative value indicates an
                                interaction for sensitivity, whereas a
                                positive value indicates an interaction
                                for resistance;
    FEATURE_pos_IC50_sd         Log IC50 Standard deviation for the
                                population of cell lines accounted in column
                                i;
    FEATURE_neg_IC50_sd         Log IC50 Standard deviation for the
                                population of cell lines accounted in column
                                j;
    FEATURE_IC50_effect_size    Cohen's d, quantifying the effect size of
                                the interaction. A value >=0.5 indicates
                                a moderate effect size. A value >=1
                                indicates a large effect size (i.e. difference
                                in mean log IC50 values greater than their
                                pooled standard deviations). A value >= 2
                                indicates a very large effect size (i.e.
                                difference in mean log IC50 is at least two
                                times their pooled standard deviation);
    FEATURE_pos_Glass_delta     Glass delta, quantifying the effect size of
                                the interaction as the ratio between the
                                difference of the mean log IC50 values and
                                the standard deviation of the log IC50 values
                                of the population of cell lines accounted in
                                column i;
    FEATURE_neg_Glass_delta     Glass delta Same as above for the negative set.
    ANOVA_FEATURE_pval          ANOVA test p-value quantyfing the interaction
                                significance;
    ANOVA_TISSUE_pval           ANOVA test p-value quantifying the
                                significance of the interaction between drug
                                response and the tissue of origin of the
                                cell lines; for the cancer-specific
                                interactions this value is NA;
    ANOVA_MEDIA_pval            ANOVA test p-value quantifying the
                                significance of the interaction between drug
                                response and the screening medium of the cell
                                lines; for the cancer-specific interactions
                                this value is NA;
    ANOVA_MSI_pval              ANOVA test p-value quantifying the
                                significance of the interaction between drug
                                response and the micro-satellite instability
                                status of the cell lines; for the cancer type
                                with no micro-satellite instable cell line
                                samples this value is NA;
    ANOVA_FEATURE_FDR           False discovery rate obtained by correcting
                                the p-values in column u, on an individual
                                analysis basis, for multiple hypothesis
                                testing with the q-value correction method
                                (Storey & TIbshirani, 2003)
    =========================== ===============================================

    Note that those column names are renamed internally (and if the data is
    saved in a new file):

    ======================= ================================
    ======================= ================================
    assoc_id                ASSOC_ID
    Drug id                 DRUG_ID
    Owned_by                OWNED_BY
    FEATUREpos_IC50_sd      FEATURE_pos_IC50_sd
    FEATUREneg_IC50_sd      FEATURE_neg_IC50_sd
    FEATUREpos_Glass_delta  FEATURE_pos_Glass_delta
    FEATUREneg_Glass_delta  FEATURE_neg_Glass_delta
    FEATUREpos_logIC50_MEAN FEATURE_pos_logIC50_MEAN
    FEATUREneg_logIC50_MEAN FEATURE_neg_logIC50_MEAN
    Drug Target             DRUG_TARGET
    FEATURE_deltaMEAN_IC50  FEATURE_delta_MEAN_IC50
    FEATURE_ANOVA_pval      ANOVA_FEATURE_pval
    ANOVA FEATURE FDR %     ANOVA_FEATURE_FDR
    MSI_ANOVA_pval          ANOVA_MSI_pval
    Tissue_ANOVA_pval       ANOVA_TISSUE_pval
    MEDIA_ANOVA_pval        ANOVA_MEDIA_pval
    TISSUE_ANOVA_pval       ANOVA_TISSUE_pval
    Drug name               DRUG_NAME
    ======================= ================================



    """
    _colname_drug_id = 'DRUG_ID'

    def __init__(self, filename=None, settings=None):
        """.. rubric:: Constructor

        :param str filename: Another ANOVAResults instance or a 
            compatible CSV file with the correct header.
            The filename may also be set to None (default) and 
            populated later.

        """
        if filename is not None and isinstance(filename, str):
            self.read_csv(filename)
        elif filename is None:
            self._df = pd.DataFrame()
        else:
            try:
                self._df = filename.df.copy()
            except:
                self._df = filename.copy()
            assert isinstance(self._df, pd.core.frame.DataFrame), \
                "excepts a dataframe or filename"

        #: dictionary with the relevant column names and their expected types
        self.mapping = OrderedDict()
        self.mapping['ASSOC_ID'] = np.dtype('int64')
        self.mapping['FEATURE'] = np.dtype('O')
        self.mapping['DRUG_ID'] = np.dtype('int64')
        self.mapping['DRUG_NAME'] = np.dtype('O')
        self.mapping['DRUG_TARGET'] = np.dtype('O')
        self.mapping['N_FEATURE_neg'] = np.dtype('int64')
        self.mapping['N_FEATURE_pos'] = np.dtype('int64')
        self.mapping['FEATURE_pos_logIC50_MEAN'] = np.dtype('float64')
        self.mapping['FEATURE_neg_logIC50_MEAN'] = np.dtype('float64')
        self.mapping['FEATURE_delta_MEAN_IC50'] = np.dtype('float64')
        self.mapping['FEATURE_IC50_effect_size'] = np.dtype('float64')
        self.mapping['FEATURE_neg_Glass_delta'] = np.dtype('float64')
        self.mapping['FEATURE_pos_Glass_delta'] = np.dtype('float64')
        self.mapping['FEATURE_neg_IC50_sd'] = np.dtype('float64')
        self.mapping['FEATURE_pos_IC50_sd'] = np.dtype('float64')
        self.mapping['FEATURE_IC50_T_pval'] = np.dtype('float64')
        self.mapping['ANOVA_FEATURE_pval'] = np.dtype('float64')
        self.mapping['ANOVA_TISSUE_pval'] = np.dtype('float64')
        self.mapping['ANOVA_MSI_pval'] = np.dtype('float64')
        self.mapping['ANOVA_MEDIA_pval'] = np.dtype('float64')
        self.mapping['ANOVA_FEATURE_FDR'] = np.dtype('float64')


        # If the dataframe is empty, we still fill the columns so that reports 
        # and other code will find the column names.
        if len(self.df) == 0:
            self.df = pd.DataFrame(columns=self.mapping.keys())
            self.df = self.astype(self.df)


        # before gdsctools, columns names were a bit different.
        # We need to rename some column names
        self.df.rename(columns={
            'assoc_id': 'ASSOC_ID',
            'Drug id': 'DRUG_ID',
            'Owned_by': 'OWNED_BY',
            'FEATUREpos_IC50_sd': 'FEATURE_pos_IC50_sd',
            'FEATUREneg_IC50_sd': 'FEATURE_neg_IC50_sd',
            'FEATUREpos_Glass_delta': 'FEATURE_pos_Glass_delta',
            'FEATUREneg_Glass_delta': 'FEATURE_neg_Glass_delta',
            'FEATUREpos_logIC50_MEAN': 'FEATURE_pos_logIC50_MEAN',
            'FEATUREneg_logIC50_MEAN': 'FEATURE_neg_logIC50_MEAN',
            'Drug Target': 'DRUG_TARGET',
            'FEATURE_deltaMEAN_IC50': 'FEATURE_delta_MEAN_IC50',
            'FEATURE_ANOVA_pval': 'ANOVA_FEATURE_pval',
            'ANOVA FEATURE FDR %': 'ANOVA_FEATURE_FDR',
            'MSI_ANOVA_pval': 'ANOVA_MSI_pval',
            'Tissue_ANOVA_pval': 'ANOVA_TISSUE_pval',
            'MEDIA_ANOVA_pval': 'ANOVA_MEDIA_pval',
            'TISSUE_ANOVA_pval': 'ANOVA_TISSUE_pval',
            'Drug name': 'DRUG_NAME', 'A': 'B'}, inplace=True)

        self.colnames_subset = [
            'ASSOC_ID', 'FEATURE',
            'DRUG_ID', 'DRUG_NAME', 'DRUG_TARGET',
            'N_FEATURE_neg', 'N_FEATURE_pos',
            'FEATURE_pos_logIC50_MEAN', 'FEATURE_neg_logIC50_MEAN',
            'FEATURE_delta_MEAN_IC50',
            'FEATURE_IC50_effect_size',
            'FEATURE_neg_Glass_delta', 'FEATURE_pos_Glass_delta',
            'ANOVA_FEATURE_pval',
            'ANOVA_TISSUE_pval',
            'ANOVA_MSI_pval',
            'ANOVA_MEDIA_pval',
            'ANOVA_FEATURE_FDR']

        self._df.reset_index(drop=True)
        self.settings = settings


    def astype(self, df):
        try:
            # does not work in python3.3 on travis but should work
            # we newer pandas version.
            df = df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
        except:
            for col in df.columns:
                if col in self.mapping.keys():
                    df[col] = df[col].astype(self.mapping[col])
        return df

    def _get_df(self):
        return self._df
    def _set_df(self, df):
        # TODO check that all columns are found and with correct type.
        self._df = df
    df = property(_get_df, _set_df, doc="dataframe with all results")

    def to_csv(self, filename):
        """Save the ANOVAResults dataframe into a CSV file"""
        assert filename.endswith('.csv'), "filename should end in .csv"
        self.df.to_csv(filename, sep=',', index=False)

    def read_csv(self, filename):
        """Read an ANOVAResults file from a CSV file

        .. todo:: check validity of the header
        """
        self.reader = readers.Reader(filename)
        self._df = self.reader.df

    def __len__(self):
        return len(self.df)

    def _get_drugIds(self):
        if len(self) == 0:
            return []
        else:
            return self.df[self._colname_drug_id].unique()
    drugIds = property(_get_drugIds,
            doc="Returns the list of drug identifiers")

    def volcano(self, settings=None):
        """Calls :class:`VolcanoANOVA` on the results

        x-value is sign(FEATURE_delta_MEAN_IC50) times FEATURE_IC50_effect_size
        y-value is the FDR correction

        See the online documentation for details on gdsctools.readthedocs.io.

        """
        if settings is None:
            settings = self.settings
        if len(self.df) == 0:
            print("No data to plot")
            return
        self.handle_volcano = VolcanoANOVA(self.df, settings=settings)
        self.handle_volcano.volcano_plot_all()

    def __str__(self):
        txt = 'Total number of ANOVA tests performed: %s ' % len(self.df)
        return txt

    def __repr__(self):
        txt = 'ANOVAResults (%s tests): ' % len(self.df)
        return txt

    def copy(self):
        """Returns a copy """
        return ANOVAResults(self.df.copy())

    def get_html_table(self, collapse_table=False, clip_threshold=2,
            index=False, header=True, escape=False):
        """Return an HTML table for the reports"""
        cmap_clip = cmap_builder('#ffffff', '#0070FF')
        cmap_absmax = cmap_builder('green', 'white', 'red')

        # The copy is used because we'll change it afterwards
        df = self.df[self.colnames_subset].copy()

        colname = 'ANOVA_FEATURE_FDR'

        df.loc[df[colname] < 0.01, colname] = '<0.01'
        # In the assoc column, we remove the first "a" letter so that
        # the column is properly sorted by Id but the link should be with the
        # "a" as prefix
        df.ASSOC_ID = df.ASSOC_ID.apply(lambda x: int(str(x).replace("a", "")))

        html = HTMLTable(df, 'notused')
        # Those columns should be links
        html.add_href("FEATURE")
        html.add_href("ASSOC_ID", url="a", suffix=".html") # here url works like a prefix
        html.add_href("DRUG_ID", url="drug_", suffix=".html") # here url works like a prefix

        for this in ['FEATURE_IC50_effect_size', 'FEATURE_neg_Glass_delta',
                'FEATURE_pos_Glass_delta']:
            html.add_bgcolor(this, cmap_clip, mode='clip',
                    threshold=clip_threshold)

        # normalise data and annotate with color
        html.add_bgcolor('FEATURE_delta_MEAN_IC50', cmap_absmax,
            mode='absmax')

        html.df.columns = [x.replace("_", " ") for x in html.df.columns]
        return html.to_html(escape=escape, header=header, index=index,
                collapse_table=collapse_table, justify='center')

    def barplot_effect_size(self):
        """Dev not for production"""
        # barplot of the IC50 effect size
        if len(self.df) == 0:
            print("No data to plot")
        data = np.sign(self.df.FEATURE_delta_MEAN_IC50) * self.df.FEATURE_IC50_effect_size
        try:
            data = data.sort_values()
        except:
            data.sort(inplace=True)

        n_green = len(data[data<0])
        n_red = len(data[data>=0])

        print(n_green, n_red)
        data.plot(kind='barh', width=1, alpha=0.5,
            color=['green']*n_green + ['red'] * n_red)

        pylab.xlabel("Effect size")
        pylab.ylabel("Drug name")



