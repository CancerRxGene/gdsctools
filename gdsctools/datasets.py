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
- :attr:`cosmic_info`

or used in analysis:

- :attr:`genomic_features_test`
- :attr:`ic50_v17`: IC50s for 1001 cell lines
- :attr:`gf_v17`: dataset with genomic features for 1001 cell lines and 680
  features (mutation, CNA)
- :attr:`ic50_v5`
- :attr:`gf_v5`
"""
# use underscore to hide from API
import easydev


__all__ = ['Data', 'ic50_test', "genomic_features_test", 
            'cosmic_info',  "cosmic_builder_test", "cancer_cell_lines"]

def _gsf(filename):
    from gdsctools import gdsctools_data
    return gdsctools_data(filename)


class Data(object):
    """A convenience class to hold information about a dataset


    Each :class:`Data` instance contains information about :

    #. the file location (:attr:`filename`)
    #. the data description (:attr:`description`)
    #. the authors (:attr:`authors`)

    But the data is not stored and users must read the data set using
    their own tools. 
    """
    def __init__(self, filename=None, description="No description",
            authors="GDSC consortium"):
        #: where is located the data set (full path)
        self.filename = filename
        #: a short description (string)
        self.description = description
        #: list of authors (string)
        self.authors = authors

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

# ============== DATA SETS DEFINITION

# Dataset with IC50s for 10 drugs (for testing)
ic50_test = Data(
        filename=_gsf("IC50_10drugs.tsv"),
        description = 'IC50s for 10 public drugs across cell lines')

# Dataset with genomic features for 1001 cell lines and 680 features
genomic_features_test = Data(
        filename = _gsf('genomic_features.tsv.gz'),
        description = 'Set of genomic features / tissue / msi')

# Dataset with cancer cell lines name / cosmic id/ tissue type and sub type
cancer_cell_lines = Data(
        filename = _gsf('cancer_cell_lines.csv'),
        description = "List of cosmic identifiers with "+\
                "the corresponding name, tissue and sub tissue types")

# Example of flat file to be read by COSMICFetcher
cosmic_builder_test = Data(
        filename = _gsf('cosmic_builder_test.txt'),
        description = "An example of flat file to be read by COSMICFetcher")

# Dataframe with COSMIC ID and their information
cosmic_info = Data(
        filename = _gsf('cosmic_info.csv.gz'),
        description = "Information about 1001 cell lines including COSMIC ID")

# IC50 from v17
ic50_v17 = Data(_gsf("IC50_v17.csv.gz"))
__all__.append("ic50_v17")

# Genomic Feature from v17
gf_v17 = Data(_gsf("genomic_features_v17.csv.gz"), 
        description="PANCAN genomic features from v17 GDSC release")
__all__.append("gf_v17")

# IC50 from v5
ic50_v5 = Data(_gsf("IC50_v5.csv.gz"))
__all__.append("ic50_v5")

# Genomic Feature from v5
gf_v5 = Data(_gsf("genomic_features_v5.csv.gz"))
__all__.append("gf_v5")



# Build a dedicate data set for testing purposes
def _build_testing():
    testing = easydev.AttrDict()
    d = Data()
    d.filename = _gsf('test_drug_decode.tsv')
    d.description = 'drug_decode in TSV format'
    testing.drug_test_tsv = d

    d = Data()
    d.filename = _gsf('test_drug_decode.csv')
    d.description = 'drug_decode in CSV format'
    testing.drug_test_csv = d

    d = Data()
    d.filename = _gsf('test_ic50_11_50.csv')
    d.description = 'A 10drug/50 cell lines IC50 test file in CSV format'
    testing.ic50_test_csv = d

    d = Data()
    d.filename = _gsf('test_genomic_features.csv')
    d.description = 'A 50 cell lines by 20 features GenomicFeature in CSV format'
    testing.genomic_features_csv = d

    d = Data()
    d.filename = _gsf('test_IC50.csv')
    d.description = 'A 10drug/1000 cell lines IC50 test file in CSV format'
    testing.ic50_test = d
    
    d = Data()
    d.filename = _gsf('test_IC50_header2.csv')
    d.description = 'An IC50 test (header with column without Drug_ prefix)'
    testing.ic50_test_header_no_drug_prefix = d

    d = Data()
    d.filename = _gsf('test_IC50_header1.csv')
    d.description = 'An IC50 test (header with column with Drug_ prefix only)'
    testing.ic50_test_header_drug_prefix_only = d

    d = Data()
    d.filename = _gsf('test_IC50_header3.csv')
    d.description = 'An IC50 test (header with mixed prefixes i.e. Drug_ or not)'
    testing.ic50_test_header_mixed_drug_prefix = d

    d = Data()
    d.filename = _gsf('test_genomic_features_bare.csv')
    d.description = "A 50 cell lines by 17 features without MSI/tissue/sample"
    testing.genomic_features_bare_csv = d


    return testing

testing = _build_testing()
