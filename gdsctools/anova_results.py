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
"""ANOVAResults stored the output of the ANOVA analysis"""
from collections import OrderedDict

import pandas as pd
import numpy as np

from gdsctools import readers
from gdsctools.volcano import VolcanoANOVA


__all__ = ['ANOVAResults']


class ANOVAResults(object):
    """Class to handle results of the ANOVA analysis

     Used to wrap the results returned by the
    :meth:`gdsctools.anova.ANOVA.anova_all` method.


    =========================== ===============================================
         column name                                  description
    =========================== ===============================================
    Domain                      The analysis in which the interaction
                                has been identified (can be PANCAN
                                for Pan-Cancer or one of the TCGA labels
                                indicating one of the cancer types
                                included in the study, see Table S1A);
    assoc_id                    Alphanumerical identifier of the
                                interaction;
    Feature                     The CFE involved in the interaction, it
                                can be a mutated cancer driver gene
                                (CG) [suffix _mut], an abberrantly
                                fused protein [suffix fusion], a copy
                                number altered chromosomal region (RACS)
                                [prefix gain for amplifications or loss
                                for deletions];
    Drug id                     Numerical id of the drug involved in the
                                interaction;
    Drug Target                 Putative target of the drug involved in the
                                interaction;
    N_FEATURE_pos:              Number of cell lines harbouring the CFE
                                indicated in column E and that have been
                                screened with the drug indicated in columns
                                F and G, therefore have been included in
                                the test;
    N_FEATURE_neg               Number of cell lines not harbouring the
                                CFE indicated in column E and that have been
                                screened with the drug indicated in columns
                                F and G, therefore have been included in the
                                test;
    log max.Conc.tested         Maximal concentration tested (values are in
                                natural log);
    log max.Conc.tested2        2nd Maximal concentration tested (values are
                                in natural log); some compounds have been
                                screened at two different ranges of
                                concentrations. When this is not the case,
                                this column contains NA;
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
    ANOVA_FEATURE_FDR_%         False discovery rate obtained by correcting
                                the p-values in column u, on an individual
                                analysis basis, for multiple hypothesis
                                testing with the q-value correction method
                                (Storey & TIbshirani, 2003)
    =========================== ===============================================

    """
    def __init__(self, filename=None):
        """.. rubric:: Constructor

        :param str filename: Another ANOVAResults instance of a saved
            dataframe that can be read by this class, that is a CSV
            with the official header. This parameter can also be set
            to None (default) and populated later.

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

        # TODO: check the header of the dataframe

        self.drug_target = 'DRUG_TARGET'
        self.drug_name = 'DRUG_NAME'

        #: dictionary with the relevant column names and their expected types
        self.mapping = OrderedDict()
        self.mapping['ASSOC_ID'] =  np.dtype('int64')
        self.mapping['FEATURE'] = np.dtype('O')
        self.mapping['DRUG_ID'] = np.dtype('O')
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
       
        # FIXME will be removed in the future
        self.mapping['log max.Conc.tested'] = np.dtype('O'),
        self.mapping['log max.Conc.tested2'] = np.dtype('O'),
        
        self.mapping['ANOVA_FEATURE_pval'] = np.dtype('float64')
        self.mapping['ANOVA_TISSUE_pval'] = np.dtype('float64')
        self.mapping['ANOVA_MSI_pval'] = np.dtype('float64')
        self.mapping['ANOVA_MEDIA_pval'] = np.dtype('float64')
  
        self.mapping['ANOVA_FEATURE_FDR_%'] = np.dtype('float64')

        self.colnames_subset = ['ASSOC_ID', 'FEATURE',
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
            'ANOVA_FEATURE_FDR_%']

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
        """Save dataframe into a file using comma separated values"""
        assert filename.endswith('.csv'), "filename should end in .csv"
        self.df.to_csv(filename, sep=',', index=False)

    def read_csv(self, filename):
        """Read a CSV file

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
            return self.df[self.drug_id].unique()
    drugIds = property(_get_drugIds,
            doc="Returns the list of drug identifiers")

    def volcano(self):
        """Calls :class:`VolcanoANOVA` on the results"""
        v = VolcanoANOVA(self.df, settings=self.settings)
        v.volcano_plot_all()



