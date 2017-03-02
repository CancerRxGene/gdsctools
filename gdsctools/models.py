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
import collections

import pandas as pd

from easydev import Progress

from gdsctools.stats import MultipleTesting
from gdsctools import readers
from gdsctools.settings import ANOVASettings
from gdsctools.anova_results import ANOVAResults
from gdsctools.errors import GDSCToolsDuplicatedDrugError

import colorlog as logger

__all__ = ['BaseModels']


class BaseModels(object): 
    """A Base class for ANOVA / ElaticNet models


    """
    def __init__(self, ic50, genomic_features=None,
            drug_decode=None, verbose=True, 
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
        self.verbose = verbose
        self._init_called = False

        # We first need to read the IC50 using a dedicated reader
        try:
            # Simple one without duplicated
            self.ic50 = readers.IC50(ic50)
        except GDSCToolsDuplicatedDrugError:  
            print("duplicated error")
            try:
                from gdsctools.gdsc import IC50Cluster
                self.ic50 = IC50Cluster(ic50)
            except Exception as err:
                raise(err)

        # Reads features if provided, otherwise use a default data set
        if genomic_features is None:
            # Reads default version provided with the package
            self.features = readers.GenomicFeatures()
        else:
            self.features = readers.GenomicFeatures(genomic_features)

        if self.features.found_media is False and \
                set_media_factor is True:
            if self.verbose:
                print('Populating MEDIA Factor in the Genomic Feature matrix')
            self.features.fill_media_factor()

        #: a CSV with 3 columns used in the report
        self.read_drug_decode(drug_decode)

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

        if len(unknowns) > 0 and self.verbose:
            print("WARNING: " +
                "%s cosmic identifiers in your IC50 " % len(unknowns) +
                "could not be found in the genomic feature matrix. " +
                "They will be dropped. Consider using a user-defined " +
                "genomic features matrix")

        self.ic50.drop_cosmic(list(unknowns))
        self.features.cosmicIds = self.ic50.cosmicIds
        #self.cosmicIds = self.ic50.cosmicIds

        #: an instance of :class:`~gdsctools.settings.ANOVASettings`
        self.settings = ANOVASettings()

        # alias to all column names to store results
        # cast to list (Python3).
        self.column_names = list(ANOVAResults().mapping.keys())

        # skip assoc_id for now
        self._odof_dict = dict([(name, None)
            for name in self.column_names])

        # a cache to store ANOVA results for each drug
        self.individual_anova = {}

        # must be called if ic50 or features are changed.
        self.init()

    def _autoset_msi_factor(self):
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
            self.settings.analysis_type = 'feature_only'

    def _autoset_tissue_factor(self):
        # select tissue based on the features
        tissue_name = self.features.colnames.tissue
        self.tissue_factor = self.features.df[tissue_name]
        if len(self.tissue_factor.unique()) == 1:
            # there is only one tissue
            tissue = self.tissue_factor.unique()[0]
            self.settings.analysis_type = tissue
            self.settings.directory = tissue
        else:
            # this is a PANCAN analysis
            self.settings.analysis_type = 'PANCAN'

    def _autoset_media_factor(self):

        if self.settings.analysis_type != 'PANCAN':
            self.settings.include_media_factor = False
            if self.features.found_media is True:
                # Not authorised. See 
                # http://gdsctools.readthedocs.io/en/master/anova_parttwo.html#regression-analysis 
                print("WARNING")
                print("You have only one Tissue %s " % self.features.tissues[0])
                print("When using MEDIA FACTOR, you must use MSI and a PANCAN analysis")
                print("We DO NOT include the MEDIA Factor in the analysis hereafter\n")

        elif self.features.found_media is True:
            self.settings.include_media_factor = True
            colname = self.features.colnames.media
            self.media_factor = self.features.df[colname]
        else:
            self.settings.include_media_factor = False

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

        if ctype == 'PANCAN':
            # Nothing to do, keep everything
            return

        if isinstance(ctype, str):
            ctype = [str(ctype)]

        for this in ctype:
            assert this in self.features.tissues, "%s not found" % ctype

        # keep only features that correspond to the tissue
        self.features.keep_tissue_in(ctype)

        self.ic50.df = self.ic50.df.ix[self.features.df.index]
        self.init()

    def read_settings(self, settings):
        """Read settings and update cancer type if set"""
        self.settings.from_json(settings)
        self.set_cancer_type(self.settings.analysis_type)

    def init(self):
        # Some preprocessing to speed up data access in ANOVA
        ic50_parse = self.ic50.df.copy().unstack().dropna()
        # for each drug, we store the IC50s (Y) and corresponding indices
        # of cosmic identifiers + since v0.13 the real indices
        # Create a dictionary version of the data
        # to be accessed per drug where NA have already been
        # removed. Each drug is a dictionary with 2 keys:
        # Y for the data and indices for the cosmicID where
        # there is an IC50 measured.
        self.ic50_dict = dict([
            (d, {'indices': ic50_parse.ix[d].index,
             'Y': ic50_parse.ix[d].values}) for d in self.ic50.drugIds])
        cosmicIds = list(self.ic50.df.index)
        for key in self.ic50_dict.keys():
            indices = [cosmicIds.index(this) for this in
                self.ic50_dict[key]['indices']]
            self.ic50_dict[key]['real_indices'] = indices

        # save the tissues
        self._autoset_tissue_factor()

        # and MSI (Microsatellite instability) status of the samples.
        self._autoset_msi_factor()

        # and (growth) media factor
        self._autoset_media_factor()

        # dictionaries to speed up code.
        self.msi_dict = {}
        self.tissue_dict = {}
        self.media_dict = {}
        # fill the dictionaries for each drug once for all
        for drug_name in self.ic50.drugIds:
            indices = self.ic50_dict[drug_name]['indices']

            # MSI, media and tissue are not large data files and can be stored
            # enterily
            if self.features.found_msi:
                self.msi_dict[drug_name] = self.msi_factor.ix[indices]

            if self.settings.include_media_factor:
                self.media_dict[drug_name] = self.media_factor.ix[indices]

            self.tissue_dict[drug_name] = self.tissue_factor.ix[indices]

        # some preprocessing for the OLS computation.
        # We create the dummies for the tissue factor once for all
        # Note that to agree with R convention, we have to resort the column
        # to agree with R convention that is a<B==b<c instead of
        # where A<B<C<a<b<c (in Python)
        self._tissue_dummies = pd.get_dummies(self.tissue_factor)
        columns = self._tissue_dummies.columns
        columns = sorted(columns, key=lambda s: s.lower())
        columns = ['C(tissue)[T.' + x + ']' for x in columns]
        self._tissue_dummies.columns = columns

        if self.settings.include_media_factor:
            self._media_dummies = pd.get_dummies(self.media_factor)
            columns = self._media_dummies.columns
            columns = ['C(media)[T.' + x + ']' for x in columns]
            self._media_dummies.columns = columns
            for col in columns:
                self._tissue_dummies[col] = self._media_dummies[col]

        N = len(self._tissue_dummies)
        self._tissue_dummies['C(msi)[T.1]'] = [1]*N
        self._tissue_dummies['feature'] = [1] * N
        self._tissue_dummies.insert(0, 'Intercept', [1] * N)

        # drop first feature in the tissues that seems to be used as a
        # reference in the regression
        #tissues = [x for x in self._tissue_dummies.columns if 'tissue' in x]
        #self._tissue_dummies.drop(tissues[0], axis=1, inplace=True)

        """if self.settings.include_media_factor:
            # Drop first category in the media factor ?! like for tissues.
            # What is the rationale ?
            media = [x for x in self._tissue_dummies.columns if 'media' in x]
            self._tissue_dummies.drop(media[0], axis=1, inplace=True)
        """
        # reset the buffer.
        self.individual_anova = {}

        if self.verbose and self._init_called is False:
            for this in ['tissue', 'media', 'msi', 'feature']:
                if this in self._get_analysis_mode():
                    logger.debug(this.upper() + " FACTOR : included")
                else:
                    logger.debug(this.upper() + " FACTOR : NOT included")
        self._init_called = True

    def _get_cosmics(self):
        return self.ic50.cosmicIds
    def _set_cosmics(self, cosmics):
        self.ic50.cosmicIds = cosmics
        self.features.cosmicIds = cosmics
        self.init()
        self.individual_anova = {}
    cosmicIds = property(_get_cosmics, _set_cosmics,
        doc="get/set the cosmic identifiers in the IC50 and feature matrices")

    def _get_drug_names(self):
        return self.ic50.drugIds
    def _set_drug_names(self, drugs):
        self.ic50.drugIds = drugs
        self.init()
        # not need to init this again ? self.individual_anova = {}
    drugIds = property(_get_drug_names, _set_drug_names,
            doc="Get/Set drug identifers")

    def _get_feature_names(self):
        shift = self.features.shift
        return self.features.features[shift:]
    def _set_features_names(self, features):
        self.features.features = features
        self.init()
        self.individual_anova = {}
    feature_names = property(_get_feature_names, _set_features_names,
            doc="Get/Set feature names")

    def _get_analysis_mode(self):
        modes = []
        if self.settings.analysis_type == 'PANCAN':
            modes.append('tissue')

        if self.settings.include_MSI_factor is True:
            modes.append('msi')

        if self.settings.include_media_factor is True:
            modes.append('media')

        modes.append('feature')
        return modes

    def diagnostics(self, details=False):
        """Return dataframe with information about the analysis

        """
        n_drugs = len(self.ic50.drugIds)
        n_features = len(self.features.features) - self.features.shift
        n_combos = n_drugs * n_features
        feasible = 0
        pb = Progress(n_drugs, 1)
        counter = 0

        feasibles = collections.defaultdict(int)

        for drug in self.ic50.drugIds:
            for feature in self.features.features[self.features.shift:]:
                dd = self._get_one_drug_one_feature_data(drug, feature,
                        diagnostic_only=True)
                if dd.status is True:
                    feasible += 1
                    feasibles[drug] +=1
            counter += 1
            pb.animate(counter)

        results = {
                'n_drug': n_drugs,
                'n_combos': n_combos,
                'feasible_tests': feasible,
                'percentage_feasible_tests': float(feasible)/n_combos*100}

        if details is True:
            results["feasible_tests_per_drug"] =  feasibles

        return results

    def read_drug_decode(self, filename=None):
        """Read file with the DRUG information

        .. seealso:: :class:`gdsctools.readers.DrugDecode`
        """
        # Read the DRUG decoder file into a DrugDecode/Reader instance
        self.drug_decode = readers.DrugDecode(filename)

    def __str__(self):
        txt = self.ic50.__str__()
        txt += "\n" + self.features.__str__()
        return txt

    def __repr__(self):
        txt = self.__str__()
        return txt


