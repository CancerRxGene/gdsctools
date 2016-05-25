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
"""Sets of miscellaneous tools"""
import os
import pylab


__all__ = ['Savefig', 'get_drug_id']


class Savefig(object):
    """A simple class to save matploltib figures in the proper place

    .. note:: For developers only
    """
    def __init__(self, verbose=False):
        self.verbose = verbose
        self._directory = None
        #: directory where to save figures
        self.directory = '.'

    def _get_directory(self):
        return self._directory
    def _set_directory(self, directory):
        self._directory = directory
        try:
            if os.path.isdir(directory) is False:
                os.mkdir(self.directory)
                if self.verbose:
                    print("Created directory {}".format(directory))
        except Exception:
            if self.verbose:
                print("Could not create the directory")
    directory = property(_get_directory, _set_directory, doc="")

    def savefig(self, name, size_inches=None, **kargs):
        """Save a matplotlib figure

        :param str filename: where to save the figure.
        :param **kargs: accepts all parameters known by pylab.savefig
        """
        # if not provided, don't do anything.
        if name is None:
            return

        try:
            directory = self.settings.directory
        except:
            directory = self.directory
        filename = directory + os.sep + name

        fig = pylab.gcf()
        oldsize = fig.get_size_inches()

        if size_inches is not None:
            fig.set_size_inches(size_inches)
        else:
            fig.set_size_inches(10, 10)
        if self.verbose:
            print("saving file in %s" % filename)

        pylab.savefig(filename, **kargs)

        # reset to original size
        fig.set_size_inches(*oldsize)
        fig.canvas.draw()


def extract_drug_identifiers(drug_lists):
    """

    .. deprecated:: 0.12 Used in IC50Cluster (version 18) for now but should be
        removed
    """
    if isinstance(drug_lists, str):
        drug_lists = [drug_lists]
        string = True
    else:
        string = False

    drug_ids = []
    for drug in drug_lists:
        if drug.lower().startswith('drug_'):
            drug_ids.append(drug.split("_")[1])
        else:
            drug_ids.append(drug.split("_")[0])
    # This should not fail !   
    drug_ids = [int(x) for x in drug_ids]

    if string is True:
        return drug_ids[0]
    else:
        return drug_ids

get_drug_id = extract_drug_identifiers
