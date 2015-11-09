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
"""Main standalone application dreamtools"""
import argparse
import sys
from easydev.console import red, purple, darkgreen
from easydev import underline


__all__ = ['anova_pipeline', 'ANOVAOptions']


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
    """This function is used by the standalone application called
    **gdsctools_anova**

    ::

        gdsctools_anova --help

    """
    msg = "Welcome to GDSCTools standalone"
    print(purple(underline(msg)))

    if args is None:
        args = sys.argv[:]

    user_options = ANOVAOptions(prog="gdsctools_anova")


    if len(args) == 1 or '--help' in args or '--test' in args:
        if '--test' in args:
            options = user_options.parse_args(args[1:])
        else:
            options = user_options.parse_args(args)
    else:
        options = user_options.parse_args(args[1:])

    if options.testing is True:
        print('Testing mode:')
        from gdsctools import ANOVA, ic50_test
        an = ANOVA(ic50_test)
        df = an.anova_one_drug_one_feature('Drug_1047_IC50', 'TP53_mut')

        assert df.loc[1,'N_FEATURE_pos'] == 554, \
            "N_feature_pos must be equal to 554"

        print(df.T)
        print(darkgreen("\nGDSCTools seems to be installed properly"))
        return 


    if options.summary is True:
        from gdsctools import anova
        an = anova.ANOVA(options.input_ic50, options.input_features)
        print(an)
        return

    # users may select a subset of drugs:
    if options.drugs is not None and options.drug is not None:
        print('Note that --include-drugs has not effect since you provided --drug')
    if options.drugs: # clean it up
        options.drugs = [x.strip(",") for x in options.drugs]
        options.drugs = [x for x in options.drugs if len(x)]



    # or a subset of features

    try:
        if options.drug is not None and options.feature is not None:
            anova_one_drug_one_feature(options)
        elif options.drug is not None:
            anova_one_drug(options)
        else: # analyse everything
            anova_all(options)
        if options.onweb is False:
            msg = "\nNote that a directory {} was created and files saved into it"
            print(purple(msg.format(options.directory)))

    except Exception as err:
        msg = """An error was caught while using gdsctools_anova.
This may be due to your input or a mis-spelled parameter (e.g. unknown feature),
or most probably a bug or lack of documentation on our side.

If you believe this is a bug, please send the error message that follows
together with your command line and input file

"""

        print(err)
        raise err


def anova_one_drug(options):
    from gdsctools import anova
    an = anova.ANOVA(options.input_ic50, options.input_features)
    an.settings.directory = options.directory
    an.settings.includeMSI_factor = options.include_msi
    an.set_cancer_type(options.tissue)

    an.settings.check()

    df = an.anova_one_drug(options.drug)

    if len(df)==0:
        print(red("\nNo valid associations tested. Please try another drug"))
        return

    df = an.add_corrected_pvalues(df)

    N = len(df)
    df.insert(0, 'assoc_ID', range(1, N+1))



    r = anova.ANOVAReport(an, results=df)
    print(darkgreen("\nCreating all figure and html documents in %s" %
            r.settings.directory))
    r.create_html_pages(onweb=options.onweb)


def anova_all(options):
    from gdsctools import anova 
    an = anova.ANOVA(options.input_ic50, options.input_features)

    if options.drugs is not None:
        an.drugIds = options.drugs
        # need to reinit, which is done when set_cancer_type is called
        # here below

    if options.features:
        an.feature_names = options.features

    an.settings.directory = options.directory
    an.settings.includeMSI_factor = options.include_msi
    an.set_cancer_type(options.tissue)
    an.settings.check()
    df = an.anova_all()
    r = anova.ANOVAReport(an, results=df)
    print("Creating all figure and html documents in %s" %
            r.settings.directory)
    r.create_html_pages(onweb=options.onweb)


def anova_one_drug_one_feature(options):
    from gdsctools import anova
    gdsc = anova.ANOVA(options.input_ic50, options.input_features)

    odof = anova.OneDrugOneFeature(gdsc,
            drug=options.drug,
            feature=options.feature,
            directory=options.directory)
    #an.factory.settings.directory = options.directory
    odof.factory.settings.includeMSI_factor = options.include_msi
    odof.factory.set_cancer_type(options.tissue)
    odof.factory.settings.check()

    # for the HTML
    odof.add_dependencies = True
    odof.add_settings = True
    df = odof.run()

    if df.ix[1]['FEATURE_IC50_effect_size'] is None:
        msg = "association %s vs %s not valid for testing (not enough" +\
              " MSI or positives for that features ? Try with "+\
              " --exclude-msi (you must then set a tissue with "+\
              " --tissue"

        print(red(msg % (options.drug, options.feature)))
    else:
        an.report(onweb=options.onweb)
    print(df.T)


class ANOVAOptions(argparse.ArgumentParser):
    """Define user interface for the gdsctools_anova standalone application


    """
    description = "tests"
    def __init__(self, version="1.0", prog=None):

        usage = """

1. analyse all IC50s data contained in <filename> and open and HTML page
   in a browser. This can be very long (5 minutes to several hours) depending
   on the size of the files:

    gdsctools_anova --input-ic50 <filename>  --onweb

2. on the same data as above, analyse only one association for a given
   drug <drug> and a given genomic features <feature>. The drug name should
   match one to be found in the header of <filename>. The feature name
   should match one of the header of the genomic feature file or if you use
   the default file, one of the feature in the default file (e.g., BRAF_mut)

   gdsctools_anova --input-ic50 <filename> --drug <drug> --feature <feature>
       --onweb

3. on the same data as above, analyse one drug across all features.

    gdsctools_anova --input-ic50 <filename> --onweb --drug <drug>

to obtain more help about the parameters, please type

    gdsctools_anova --help
    """

        epilog = """
Author(s): Thomas Cokelaer (GDSCtools) and authors from the GDSCtools repository.

How to contribute ? : Visit https://github.com/CancerRxGene/gdsctools
Issues or bug report ? Please fill an issue on
http://github.com/CancerRxGene/gdsctools/issues """

        description = """General Description:"""
        # FIXME : not robust but will work for now
        super(ANOVAOptions, self).__init__(usage=usage, version=version,
                prog=prog,
                epilog=epilog, description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
        self.add_input_options()

    def add_input_options(self):
        """The input options to gdsctools_anova are defined here"""
        group = self.add_argument_group("General", 
                                        'General options (compulsary or not)')

        group.add_argument("-I", "--input-ic50", dest='input_ic50',
                           default=None, type=str,
                           help="""A file in TSV format with IC50s.
                           First column should be the COSMIC identifiers
                           Following columns contain the IC50s for a set of
                           drugs. The header must
                           be COSMIC ID, Drug_1_IC50, Drug_2_IC50, ... """)
        group.add_argument("-F", "--input-features", dest='input_features',
                           default=None, type=str,
                           help="""A matrix of genomic features. First column
                           is made of COSMIC identifiers that should match
                           those from the IC50s matrix. Then 3 compulsary
                           columns are required named 'Sample', 'Tissue',
                           'MSI'. Other
                           columns should be the genomic features. There are
                           recognised if the (1) ends in _mut for mutation, or
                           starts with loss or gain for the CNA cases.
                           """)

        group.add_argument("--output-directory", default='html_gdsc_anova',
                           dest='directory',
                           help="""directory where to save images and HTML
                           files.""")
        group.add_argument( "--verbose", dest='verbose',
                           action="store_true",
                           help="verbose option.")
        group.add_argument("--on-web", dest='onweb',
                           action="store_true",
                           help="TODO")
        group.add_argument("--onweb", dest='onweb',
                           action="store_true",
                           help="same as -on-web")

        # if one drug one feature only
        group.add_argument("-d", "--drug", dest="drug",
                           help="""The name of a valid drug identifier to be
                           found in the header of the IC50 matrix""")
        group.add_argument("-f", "--feature", dest="feature",
                           help="""The name of a valid feature to be found in
                          the Genomic Feature matrix""")

        # selector tissue
        group.add_argument("-t", "--tissue", dest="tissue", type=str,
                           help="""The name of a specific cancer type
                          i.e., tissue to restrict the analysis
                          to """)
        group.add_argument("--exclude-msi", dest="include_msi",
                           action="store_false",
                           help="Include msi factor in the analysis")
        group.add_argument("--summary", dest="summary",
                            action="store_true",
                           help="Print summary about the data (e.g., tissue)")

        group.add_argument('--include-drugs-in', dest='drugs',
                           nargs="+", default=[],
                           help="todo"
                           )

        group.add_argument('--include-features-in', dest='features',
                           nargs="+", default=[],
                           help="todo"
                           )

        group.add_argument('--test', dest='testing',
                           action="store_true",
                           help="todo"
                           )




