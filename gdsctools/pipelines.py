# -*- python -*-
#
#  This file is part of DREAMTools software
#
#  Copyright (c) 2014-2015 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.org/dreamtools
#
##############################################################################
"""Main standalone application dreamtools"""
import os
import argparse
import sys
from easydev.console import red, purple, darkgreen
from easydev import DevTools
from gdsctools import anova, reader



# ------------------------------------------------ The User Interface
def print_color(txt, func_color, underline=False):
    import easydev
    try:
        if underline:
            print(easydev.underline(func_color(txt)))
        else:
            print(func_color(txt))
    except:
        print(txt)


def anova_pipeline(args=None):
    """This function is used by the standalone application called dreamscoring

    ::

        dreamscoring --help

    """
    d = DevTools()

    if args is None:
        args = sys.argv[:]
    user_options = ANOVAOptions(prog="gdsctools_anova")

    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
        options = user_options.parse_args(args[1:])

    r = reader.IC50(options.data)
    an = anova.GDSC_ANOVA(r.ic50, r.features)
    df = an.anova_one_drug_one_feature(options.drug, 
            feature_name=options.feature, 
            show_boxplot=options.show_boxplots)
    print(df.T)


class ANOVAOptions(argparse.ArgumentParser):
    description = "tests"
    def __init__(self, version="1.0", prog=None):

        usage = """usage: python %s --challenge d8c1 --sub-challenge sc1a --submission <filename>\n""" % prog
        usage += """      python %s --challenge d5c2 --submission <filename>""" % prog
        epilog="""Author(s): Thomas Cokelaer (GDSCtools) and authors from
the GDSCtools repository. .

Source code on: https://github.com/CancerRxGene/gdsctools
Issues or bug report ? Please fill an issue on
http://github.com/CancerRxGene/gdsctools/issues """
        description = """General Description:

            TODO

"""
        # FIXME : not robust but will work for now
        super(ANOVAOptions, self).__init__(usage=usage, version=version, prog=prog,
                epilog=epilog, description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
        self.add_input_options()

    def add_input_options(self):
        """The input oiptions.

        Default is None. Keep it that way because otherwise, the contents of
        the ini file is overwritten in :class:`apps.Apps`.
        """
        group = self.add_argument_group("General", 'General options (compulsary or not)')

        #group.add_argument("--ic50", dest='ic50',
        #                 default=None, type=str,
        #                 help="todo")
        #group.add_argument("--features", dest='features',
        #                 default=None, type=str,
        #                 help="todo")
        group.add_argument("--data", dest='data',
                         default=None, type=str,
                         help="A TSV file with cosmic Ids as rows and Drug and"
                         + " features as columns. ")
        group.add_argument("--verbose", dest='verbose',
                         action="store_true",
                         help="verbose option.")
        group.add_argument("--show-boxplots", dest='show_boxplots',
                         action="store_true",
                         help="show images.")
        group.add_argument("--drug", dest="drug",
                         help="")
        group.add_argument("--feature", dest="feature",
                         help="")


