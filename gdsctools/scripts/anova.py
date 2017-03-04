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
"""Main standalone application called gdsctools_anova"""
import argparse
import sys
from easydev.console import red, purple, darkgreen
import easydev
import gdsctools

import warnings
# ignore pandas warning
warnings.simplefilter(action="ignore", category=FutureWarning)
# ignore mpld3 warning
#warnings.simplefilter(action="ignore", category=UserWarning)


__all__ = ['main', 'ANOVAOptions']


def print_color(txt, func_color=darkgreen, underline=False):
    try:
        if underline:
            print(easydev.underline(func_color(txt)))
        else:
            print(func_color(txt))
    except:
        print(txt)


def main(args=None):
    """This function is used by the standalone application called
    **gdsctools_anova**

    Type::

        gdsctools_anova --help

    to get some help.
    """
    msg = "Welcome to GDSCTools standalone"
    print_color(msg, purple, underline=True)

    # Keep the argument args as None by default to
    # allow testing e.g., in nosetests
    if args is None:
        args = sys.argv[:]
    elif len(args) == 1:
        args += ['--help']

    user_options = ANOVAOptions(prog="gdsctools_anova")
    try:
        options = user_options.parse_args(args[1:])
    except SystemExit:
        return

    # -----------------------------------------------------------------
    # ---------------------------------------- options without analysis
    # -----------------------------------------------------------------

    if options.version is True:
        print("This is version %s of gdsctools_anova" % gdsctools.version)
        return

    if options.testing is True:
        print('Testing mode:')
        from gdsctools import ANOVA, ic50_test
        an = ANOVA(ic50_test)
        df = an.anova_one_drug_one_feature(1047, 'TP53_mut')

        assert df.loc[1,'N_FEATURE_pos'] == 554, \
            "N_feature_pos must be equal to 554"
        print(df.T)
        print(darkgreen("\nGDSCTools seems to be installed properly"))
        return

    if options.save_settings:
        from gdsctools import ANOVA, ic50_test
        an = ANOVA(ic50_test)
        an.settings.to_json(options.save_settings)
        print('Save a default parameter set in %s' % options.save_settings)
        return

    if options.summary is True:
        from gdsctools import anova
        an = anova.ANOVA(options.input_ic50, options.input_features)
        print(an)
        return

    if options.print_tissues is True:
        from gdsctools import anova
        an = anova.ANOVA(options.input_ic50, options.input_features)

        tissues = an.tissue_factor
        try:
            tissues = tissues.sort_values().unique()
        except Exception as err:
            print(err)
            tissues = tissues.sort().unique()
        for name in tissues:
            print(name)
        return

    if options.print_drugs is True:
        print_drugs(options)
        return

    if options.print_features is True:
        from gdsctools import anova
        gdsc = anova.ANOVA(options.input_ic50, options.input_features)
        import textwrap
        print("\n".join(textwrap.wrap(" , ".join(gdsc.feature_names))))
        return

    # -----------------------------------------------------------------
    # --------------------------------------------------- real analysis
    # -----------------------------------------------------------------
    # dispatcher to the functions according to the user parameters
    from gdsctools import ANOVA, ANOVAReport
    anova = ANOVA(options.input_ic50, options.input_features,
            options.input_drug_decode)
    anova = _set_settings(anova, options)

    if options.drug and options.drug not in anova.ic50.df.columns:
        print(red("Invalid Drug. Try one of those"))
        print_drugs(options)
        sys.exit(1)

    if options.drug is not None and options.feature is not None:
        print_color("ODOF mode", purple)
        anova_one_drug_one_feature(anova, options)
    elif options.drug is not None:
        print_color("ODAF mode", purple)
        anova_one_drug(anova, options)
    else: # analyse everything
        if options.feature is None:
            print_color("ADAF mode", purple)
            anova_all(anova, options)
        else:
            print("You provided --feature but can be used only with --drug")
            sys.exit(0)

    if options.onweb is False and options.no_html is False:
        msg = "\nNote that a directory {} was created and files saved into it"
        print(purple(msg.format(options.directory)))

    return

def _set_settings(gdsc, options):
    if options.settings is not None:
        gdsc.settings.from_json(options.settings)

    # by defaul MSI is included. It will be set to False automatically if
    # there is not enough data. One can set it to False if provided but cannot
    # be forced to True (may not be correct if there is not enough data) hence
    # the condition here below
    if options.include_msi is False:
        gdsc.settings.include_MSI_factor = False
    gdsc.settings.directory = options.directory
    gdsc.settings.FDR_threshold = options.FDR_threshold
    gdsc.settings.check()
    return gdsc


def anova_one_drug(anova, options):
    """Analyse one specific drug"""
    from gdsctools import ANOVAReport
    anova.set_cancer_type(options.tissue)

    if options.feature:
        anova.feature_names = options.features

    results = anova.anova_one_drug(int(options.drug))

    print("\nFound %s possible associations" % len(results))
    if len(results)==0:
        print(red("\nPlease try with another drug or no --drug option"))
        return

    # ?? is this required ? It looks like (May 2016)
    N = len(results)
    results.df.insert(0, 'ASSOC_ID', range(1, N+1))

    if options.no_html is True:
        return

    r = ANOVAReport(anova, results=results)
    print(darkgreen("\nCreating all figure and html documents in %s" %
            r.settings.directory))
    r.create_html_pages(onweb=options.onweb)


def anova_all(anova, options):
    """Analyse the entire data set. May be restricted to one feature"""
    from gdsctools import ANOVAReport
    if options.feature:
        anova.feature_names = [options.feature]

    # The analysis
    print(darkgreen("Starting the analysis"))
    df = anova.anova_all()
    if len(df) == 0:
        print("Found no valid association ? Check your input files")
        return
    # HTML report
    if options.no_html is True:
        return

    r = ANOVAReport(anova, results=df)
    print(darkgreen("Creating all figure and html documents in %s" %
            r.settings.directory))
    r.create_html_pages(onweb=options.onweb)


def anova_one_drug_one_feature(anova, options):
    """Analyse the entire data set"""

    from gdsctools import anova_report
    from gdsctools.report import ReportMain

    if options.tissue is not None:
        anova.set_cancer_type(options.tissue)


    # just to create the directory
    ReportMain(directory=options.directory)

    odof = anova_report.Association(anova,
            drug=int(options.drug),
            feature=options.feature)
    df = odof.run()

    if df.ix[1]['FEATURE_IC50_effect_size'] is None:
        msg = "association %s vs %s not valid for testing (not enough" +\
              " MSI or positives for that features ? Try with "+\
              " --exclude-msi (you must then set a tissue with "+\
              " --tissue"
        print(red(msg % (options.drug, options.feature)))
    else:
        print(df.T)
        # HTML report
        if options.no_html is True:
            return
        odof.create_report(onweb=options.onweb)


class ANOVAOptions(argparse.ArgumentParser):
    """Define user interface for the gdsctools_anova standalone application

    Type::

        gdsctools_anova --help

    in a shell to get detailled help about the parameters and usage.
    """
    def __init__(self, prog=None):

        usage = """

1. analyse all IC50s data contained in <filename> and open and HTML page
   in a browser. This can be very long (5 minutes to several hours) depending
   on the size of the files:

    gdsctools_anova --input-ic50 <filename>

2. on the same data as above, analyse only one association for a given
   drug <drug> and a given genomic features <feature>. The drug name should
   match one to be found in the header of <filename>. The feature name
   should match one of the header of the genomic feature file or if you use
   the default file, one of the feature in the default file (e.g., BRAF_mut)

   gdsctools_anova --input-ic50 <filename> --drug <drug> --feature <feature>

3. on the same data as above, analyse one drug across all features.

    gdsctools_anova --input-ic50 <filename> --drug <drug>

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
        super(ANOVAOptions, self).__init__(usage=usage,
                prog=prog, epilog=epilog, description=description,
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

        group.add_argument("--output-directory", default='html_gdsc_anova',
                           dest='directory',
                           help="""directory where to save images and HTML
                           files.""")

        group.add_argument( "--verbose", dest='verbose',
                           action="store_true",
                           help="verbose option.")

        group.add_argument( "--do-not-open-report", dest='onweb',
                           action="store_false",
                           help="""By default, opens the index.html page.
                           Set this option if you do not want to open the
                           html page automatically.""")
        # if one drug one feature only
        group.add_argument("-d", "--drug", dest="drug", type=int,
                           help="""The name of a valid drug identifier to be
                           found in the header of the IC50 matrix""")

        group.add_argument("-f", "--feature", dest="feature",
                           help="""The name of a valid feature to be found in
                          the Genomic Feature matrix""")

        # information
        group.add_argument("--print-drug-names", dest="print_drugs",
                           action="store_true",
                           help="Print the drug names")
        group.add_argument("--print-feature-names", dest="print_features",
                           action="store_true",
                           help="Print the features names")
        group.add_argument("--print-tissue-names", dest="print_tissues",
                           action="store_true",
                           help="Print the unique tissue names")

        # Various filters
        group.add_argument("-t", "--tissue", dest="tissue", type=str,
                           help="""The name of a specific cancer type
                          i.e., tissue to restrict the analysis
                          to. """)

        group.add_argument('--FDR-threshold', dest="FDR_threshold",
                            default=25, type=float,
                            help="""FDR (False discovery Rate) used in the
                            multitesting analysis to correct the pvalues""")

        group.add_argument("--exclude-msi", dest="include_msi",
                           action="store_false",
                           help="Include msi factor in the analysis")

        group.add_argument("--save-settings", dest="save_settings",
                            type=str,
                            help="Save settings into a json file")

        group.add_argument("--read-settings", dest="settings",
                            type=str,
                            help="""Read settings from a json file. Type
                            --save-settings <filename.json> to create
                            an example. Note that the FDR-threshold
                            and include_MSI_factor will be replaced if
                            --exclude-msi or fdr-threshold are used.""")

        # others
        group.add_argument("--summary", dest="summary",
                            action="store_true",
                           help="Print summary about the data (e.g., tissue)")

        group.add_argument('--test', dest='testing',
                           action="store_true",
                           help="""Use a small IC50 data set and run the
                           one-drug-one-feature analyse with a couple of unit
                           tests."""
                           )

        group.add_argument('--no-html', dest='no_html',
                           action="store_true",
                           help="""If set, no images or HTML are created. For
                           testing only""")

        group.add_argument('--version', dest='version',
                           action="store_true",
                           help="print current version of this application")




def print_drugs(options):
    from gdsctools import anova
    import textwrap
    gdsc = anova.ANOVA(options.input_ic50, options.input_features)
    print("\n".join(textwrap.wrap(" , ".join([str(x) for x in gdsc.drugIds]))))
