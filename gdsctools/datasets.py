# coding=utf-8
# -*- python -*-
#
#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.com>
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
"""Data sets provided with GDSCTools

The datasets may be for testing purpose:

- :attr:`ic50_test`
- :attr:`drug_test`
- :attr:`cosmic_builder_test`

or informative:

- :attr:`cancer_cell_lines`

or used in analysis:

- :attr:`genomic_features`
"""
# use underscore to hide from API
import easydev
from easydev import get_share_file as _gsf


__all__ = ['Data', 'dataset', 'ic50_test', 'genomic_features', 'cosmic_info',
           "cosmic_builder_test", "cancer_cell_lines"]


class Data(object):
    """A convenience class to hold information about a dataset


    Each :class:`Data` instance contains information about :

    #. the file location (:attr:`filename`)
    #. the data description (:attr:`description`)
    #. the authors (:attr:`authors`)

    But the data is not stored and users must read the data set using
    their own tools. 
    """
    def __init__(self):
        #: where is located the data set (full path)
        self.filename = None 
        #: a short description (string)
        self.description = "No description"
        #: list of authors (string)
        self.authors = 'GDSC consortium'

    def __str__(self):
        txt = 'location: %s\n' % self.filename
        txt += 'description: %s\n' % self.description
        txt += 'authors: %s\n' % self.authors
        return txt

    def _get_location(self):
        return self.filename
    location = property(_get_location)

    def __repr__(self):
        return self.__str__()


def dataset(dataname):
    """Retrieve information about a dataset including location

    :param dataname: a data set's name (e.g., ic50_test)
    :return: a :class:`Data` holder

    Get information about a dataset and in particular its physical location
    ::

        from gdsctools.datasets import ic50_test
        print(i50_test)
        # Get its location
        ic50_test.filename

    """

    d = Data()
    if dataname == 'ic50_test':
        d.filename = _gsf('gdsctools', 'data', 'IC50_10drugs.tsv')
        d.description = 'IC50s for 10 public drugs across cell lines'
    elif dataname == 'genomic_features':
        d.filename = _gsf('gdsctools', 'data', 'genomic_features.tsv.gz')
        d.descritption = 'Set of genomic features + tissue + sample name + msi'
    elif dataname == 'cancer_cell_lines':
        d.filename = _gsf('gdsctools', 'data', 'cancer_cell_lines.csv')
        d.description = "List of cosmic identifiers with the corresponding "+\
            "name, tissue and sub tissue types"
    elif dataname == 'cosmic_builder_test':
        d.filename = _gsf('gdsctools', 'data', 'cosmic_builder_test.txt')
        d.description = "An example of flat file to be read by COSMICFetcher"
    elif dataname == 'cosmic_info':
        d.filename = _gsf('gdsctools', 'data', 'cosmic_info.csv.gz')
        d.description = "Information about 1001 cell lines including COSMIC ID"

    return d

# ALIASES

#: Dataset with IC50s for 10 drugs (for testing)
ic50_test = dataset('ic50_test')

#: Dataset with genomic features for 1001 cell lines and 680 features
genomic_features = dataset('genomic_features')

#: Dataset with cancer cell lines name / cosmic id/ tissue type and sub type
cancer_cell_lines = dataset('cancer_cell_lines')

#: Example of flat file to be read by COSMICFetcher
cosmic_builder_test = dataset("cosmic_builder_test")

#: Dataframe with COSMIC ID and their information
cosmic_info = dataset("cosmic_info")



# Build a dedicate data set for testing purposes
def _build_testing():
    testing = easydev.AttrDict()
    d = Data()
    d.filename = _gsf('gdsctools', 'data', 'test_drug_decode.tsv')
    d.description = 'drug_decode in TSV format'
    testing.drug_test_tsv = d

    d = Data()
    d.filename = _gsf('gdsctools', 'data', 'test_drug_decode.csv')
    d.description = 'drug_decode in CSV format'
    testing.drug_test_csv = d

    d = Data()
    d.filename = _gsf('gdsctools', 'data', 'test_ic50_11_50.csv')
    d.description = 'A 10drug/50 cell lines IC50 test file in CSV format'
    testing.ic50_test_csv = d

    d = Data()
    d.filename = _gsf('gdsctools', 'data', 'test_genomic_features.csv')
    d.description = 'A 50 cell lines by 20 features GenomicFeature in CSV format'
    testing.genomic_features_csv = d

    d = Data()
    d.filename = _gsf('gdsctools', 'data', 'test_IC50.csv')
    d.description = 'A 10drug/1000 cell lines IC50 test file in CSV format'
    testing.ic50_test = d
    
    d = Data()
    d.filename = _gsf('gdsctools', 'data', 'test_IC50_header2.csv')
    d.description = 'An IC50 test (header with column without Drug_ prefix)'
    testing.ic50_test_header_no_drug_prefix = d

    d = Data()
    d.filename = _gsf('gdsctools', 'data', 'test_IC50_header1.csv')
    d.description = 'An IC50 test (header with column with Drug_ prefix only)'
    testing.ic50_test_header_drug_prefix_only = d

    d = Data()
    d.filename = _gsf('gdsctools', 'data', 'test_IC50_header3.csv')
    d.description = 'An IC50 test (header with mixed prefixes i.e. Drug_ or not)'
    testing.ic50_test_header_mixed_drug_prefix = d

    d = Data()
    d.filename = _gsf('gdsctools', 'data', 'test_genomic_features_bare.csv')
    d.description = "A 50 cell lines by 17 features without MSI/tissue/sample"
    testing.genomic_features_bare_csv = d


    return testing

testing = _build_testing()
