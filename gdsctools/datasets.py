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
"""Data sets provided with GDSCTools"""
import easydev


__all__ = ['dataset', 'ic50_test', 'genomic_features', 'drug_test']


class Data(object):
    """A class to hold information about a dataset


    Can be used as input to :class:`ANOVA` instance.
    """
    def __init__(self):
        #: where is located the data set
        self.filename = None
        #: a short description
        self.description = "No description"
        #: list of authors
        self.authors = 'GDSC consortium'

    def __str__(self):
        txt = 'location: %s\n' % self.filename
        txt += 'description: %s\n' % self.description
        txt += 'authors: %s\n' % self.authors
        return txt


registered = ['ic50_test', 'genomic_features', 'drug_test',
    'cancer_cell_line']


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
    easydev.check_param_in_list(dataname, registered)

    d = Data()
    if dataname == 'ic50_test':
        d.filename = easydev.get_share_file('gdsctools', 'data', 'IC50_10drugs.tsv')
        d.description = 'IC50s for 10 public drugs across cell lines'
    elif dataname == 'genomic_features':
        d.filename = easydev.get_share_file('gdsctools',
                'data', 'genomic_features.tsv')
        d.descritption = 'Set of genomic features + tissue + sample name + msi'
    elif dataname == 'drug_test':
        d.filename = easydev.get_share_file('gdsctools',
                'data', 'DRUG_DECODE.csv')
        d.description = "Mapping between drug identifiers, drug " +\
                         "name and drug target"
    elif dataname == 'cancer_cell_lines':
        d.filename = easydev.get_share_file('gdsctools',
                'data', 'cancer_cell_lines.csv')
        d.description = "List of cosmic identifiers with the corresponding "+\
            "name, tissue and sub tissue types"
    return d

#: dataset with IC50s for 10 drugs (for testing)
ic50_test = dataset('ic50_test')

#: dataset with genomic features for 1001 cell lines and 680 features
genomic_features = dataset('genomic_features')

#: dataset with drug name and targets (for testing)
drug_test = dataset('drug_test')


#: cancer cell lines table with name / cosmic id/ tissue type and sub type
cancer_cell_lines = dataset('cancer_cell_line')
