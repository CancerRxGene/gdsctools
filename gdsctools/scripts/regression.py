# -*- python -*-
# -*- coding utf-8 -*-

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
"""Main standalone application called gdsctools_regression"""
import shutil
import os
import argparse
import sys
from easydev.console import red, purple, darkgreen
import easydev
import gdsctools

import warnings
# ignore pandas warning
warnings.simplefilter(action="ignore", category=FutureWarning)


__all__ = ['main', 'RegressionOptions']


def print_color(txt, func_color=darkgreen, underline=False):
    try:
        if underline:
            print(easydev.underline(func_color(txt)))
        else:
            print(func_color(txt))
    except:
        print(txt)


def main(args=None):
    """This function is used by **gdsctools_regression**

    Type::

        gdsctools_regression --help

    to get some help.
    """
    msg = "Welcome to GDSCTools standalone (lasso, ridge, elastic net)"
    print_color(msg, purple, underline=True)

    # Keep the argument args as None by default to
    # allow testing e.g., in nosetests
    if args is None:
        args = sys.argv[:]
    elif len(args) == 1:
        args += ['--help']

    user_options = RegressionOptions(prog="gdsctools_regression")
    try:
        options = user_options.parse_args(args[1:])
    except SystemExit:
        return

    # -----------------------------------------------------------------

    if options.version is True:
        print("This is version %s of gdsctools_regression" % gdsctools.version)
        return

    if options.license is True:
        print(gdsctools.license)
        return

    # -------------------------------------------- real analysis --------
    try:
        os.mkdir(options.output_directory)
    except:
        if options.force is False:
            print_color(("directory already exists, erase it or choose a different "
               " output directory with --output-directory. Use --force to force"
                " your choice"), "red")
            return
        else:
            pass # keep going

    # Copy the regression pipeline
    filename = gdsctools.gdsctools_data("regression.rules", '../pipelines')
    #gdsctools_path = easydev.get_package_location('gdsctools')
    #filename = os.sep.join([gdsctools_path, "gdsctools", "pipelines",
    #    "regression.rules"])
    shutil.copy(filename, options.output_directory)

    # create the config
    params = {"method":options.method,
              "kfold": options.kfold,
              "ic50": os.path.realpath(options.input_ic50),
              "features": os.path.realpath(options.input_features)}

    config_template = """
# Analysis
regression:
  method: %(method)s       # lasso, elasticnet or ridge
  kfold: %(kfold)s        # Used to automatically estimate best alpha parameter
  randomness: 50
  boxplot_n: 5

# Input data sets
input:
  ic50: %(ic50)s
  genomic_features: %(features)s
    """

    with open(options.output_directory + os.sep + "config.yaml", "w") as fh:
        fh.write(config_template % params)


    print("File config.yaml and regression.rules created in ./%s" %
        options.output_directory)

    print("First go to the directory where analysis will be performed\n\n")
    print("    cd %s\n" % options.output_directory)
    msg = """You have two choices now. Either you are on a laptop, or you are on
a cluster.

1. LOCAL COMPUTER:
------------------

    snakemake -s regression.rules -p

where -p means 'print statements'

2. CLUSTERS:
------------

On a SLURM cluster, you can make use of the many cores available by typing for
instance:

    srun --qos normal snakemake -s regression.rules -j 40 --cluster "sbatch --qos normal"

For more information about snakemake commands, type

    snakemake --help

"""
    print(msg)
    with open(options.output_directory + os.sep + "README", "w") as fh:
        fh.write(msg)


class RegressionOptions(argparse.ArgumentParser):
    """Define user interface for the gdsctools_regression standalone application

    Type::

        gdsctools_regression --help

    in a shell to get detailled help about the parameters and usage.
    """
    def __init__(self, prog=None):

        usage = """

    gdsctools_regression -I <filename> -F <filename> --method lasso


    """

        epilog = """
Author(s): Thomas Cokelaer (GDSCtools) and authors from the GDSCtools repository.

How to contribute ? : Visit https://github.com/CancerRxGene/gdsctools
Issues or bug report ? Please fill an issue on
http://github.com/CancerRxGene/gdsctools/issues """

        description = """General Description:"""
        # FIXME : not robust but will work for now
        super(RegressionOptions, self).__init__(usage=usage,
                prog=prog, epilog=epilog, description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
        self.add_input_options()

    def add_input_options(self):
        """The input options to gdsctools_regression are defined here"""
        group = self.add_argument_group("General",
                                        'General options (compulsary or not)')

        group.add_argument("-I", "--input-ic50", dest='input_ic50',
                           default=None, type=str,
                           help="""A file in TSV format with IC50s.
                           First column should be the COSMIC identifiers
                           Following columns contain the IC50s for a set of
                           drugs. The header must
                           be COSMIC_ID, Drug_1_IC50, Drug_2_IC50, ... """)

        group.add_argument("-F", "--input-features", dest='input_features',
                           default=None, type=str,
                           help="""A matrix of genomic features. One column
                           with COSMIC identifiers should match
                           those from the IC50s matrix. Then columns named
                           TISSUE_FACTOR, MSI_FACTOR, MEDIA_FACTOR should
                           be provided. Finally, other columns will be
                           considered as genomic features (e.g., mutation)
                           """)

        group.add_argument("-D", "--input-drug-decode", dest='input_drug_decode',
                            default=None, type=str,
                            help="a decoder file")

        group.add_argument("-K", "--kfold", dest='kfold',
                            default=10, type=int,
                            help="kfold for regression cross validation")

        group.add_argument("-O", "--output-directory", default='analysis',
                           dest='output_directory',
                           help="""directory where to save images and HTML
                           files.""")

        group.add_argument("-v", "--verbose", dest='verbose',
                           action="store_true",
                           help="verbose option.")

        group.add_argument("-f", "--force", dest='force',
                           action="store_true",
                           help="""force creation of the directory and overwrite
                           files.""")

        group.add_argument("-M", "--method", dest="method",
                           default="lasso",
                           help="lasso, elasticnet or ridge")

        group.add_argument('--license', dest='license',
                           action="store_true",
                           help="Print the current license"
                           )

        group.add_argument('--version', dest='version',
                           action="store_true",
                           help="print current version of this application")



