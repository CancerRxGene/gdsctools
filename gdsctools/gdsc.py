import glob
import os

from gdsctools.anova import ANOVA
from gdsctools.readers import IC50
from gdsctools.readers import DrugDecode
from gdsctools.anova_results import ANOVAResults
from gdsctools.anova_report import ANOVAReport
from gdsctools.settings import ANOVASettings
from gdsctools.anova_report import ReportMAIN

import pandas as pd

from easydev.console import purple, brown, red
from easydev import Progress

from reports import HTMLTable

__all__ = ["IC50Cluster", "GDSC"]


class IC50Cluster(IC50):
    """

    .. todo:: doc

    Used in v18 onwards to cluster the DRUG columns with same drug Identifiers.

    The convention in v18 is not great. Drug ids should be unique anyway.

    This class merges drug names (not identifiers) but duplicated identifiers
    may appear.


    ::

        ic50 = IC50Cluster(gdsctools_data("test_v18_clustering.tsv"))

        # This may not work if there are duplicated drug identifiers
        IC50(ic50)

    """
    def __init__(self, ic50, ratio_threshold=10, verbose=True):
        """

        """
        super(IC50Cluster, self).__init__(ic50, v18=True)
        self.verbose = verbose
        self.ratio_threshold = ratio_threshold

        if self.verbose:
            print(self)
        self.cluster()
        self.cleanup()
        if self.verbose:
            print(self)

    def _get_to_cluster(self):
        info = self._info()
        if len(info) > 0:
            to_cluster = info[info.ratio < 10].DRUG_ID.values
            return list(to_cluster)
        else:
            return []
    to_cluster = property(_get_to_cluster)

    def _get_mapping(self):
        from collections import defaultdict
        mapping = defaultdict(list)

        drug_ids = [self.drug_name_to_int(x) for x in self.df.columns]
        for drug_id, colname in zip(drug_ids, self.df.columns):
            mapping[drug_id].append(colname)
        return mapping

    def _get_duplicated(self):
        mapping = self._get_mapping()
        duplicated = [key for key in mapping.keys() if len(mapping[key]) > 1]
        return duplicated
    duplicated = property(_get_duplicated)

    def _info(self):
        mapping = self._get_mapping()
        duplicated = self.duplicated

        max_ids = max([len(x) for x in mapping.values()])

        results = []
        for drug_id in duplicated:
            df = self.df[mapping[drug_id]]

            total = df.mean(axis=1).count()
            individuals = list(df.count().values)

            # the number of individuals may be different. In v18 one drug had 3
            # entries (all other duplicated had only 2). so, we need to add
            # dummies (NA) when required:
            individuals += [None] * (max_ids - len(individuals))

            common = sum(df.count(axis=1) >= 2)
            result = [drug_id] + individuals + [total, common,
                    100 * common/float(total)]
            results.append(result)

        df = pd.DataFrame(results)
        if len(df):
            df.columns = ['DRUG_ID'] + [str(x) for x in range(1, max_ids+1)] +\
                ['total', 'common', 'ratio']
            try:
                df.sort_values(by='ratio', ascending=False, inplace=True)
            except:
                df.sort('ratio', ascending=False, inplace=True)
        return df

    def cluster(self):
        # get list of drug identifiers to cluster

        to_cluster = self.to_cluster
        self.clustered = to_cluster[:]
        mapping = self._get_mapping()
        self.mapped = {}

        if self.verbose:
            print('Found  %s non unique drug identifiers ' %
                    len(self.duplicated))
            print('Clustering %s of them.\n' % len(self.to_cluster))
        if len(self.to_cluster) == 0:
            return

        for identifier in to_cluster:
            drug_names = mapping[identifier]

            # Let us keep only the first concentration for now
            new_drug_name = drug_names[0]
            if len(drug_names) > 1:
                todrop = drug_names[1:]
            # add new column with new name and mean of the columns with same
            # drug id
            self.df[new_drug_name] = self.df[drug_names].mean(axis=1)
            # Remove the individual columns
            self.df.drop(todrop, axis=1, inplace=True)

    def cleanup(self):
        # Need to transform column names in proper identifiers (integer)
        # and makes sure identifiers are unique. If not, we add +10000
        # Also, for later we keep track of the original mame in a dictionary
        self.extra_mapping = {}
        new_columns = []
        for col in self.df.columns:
            identifier = self.drug_name_to_int(col)
            while identifier in new_columns:
                identifier += 10000 # not robust but would do for now
                # We use a while since ids may occur 3 times
            self.extra_mapping[identifier] = col
            new_columns.append(identifier)
        self.df.columns = new_columns


class GDSCBase(object):
    def __init__(self, genomic_feature_pattern="GF_*csv", verbose=True):
        self.verbose = True
        self.gf_filenames = glob.glob(genomic_feature_pattern)
        if len(self.gf_filenames) == 0:
            msg = "NO Genomic feature input files found. We expect files " +\
                    "with this pattern: GF_<TCGA>.csv e.g., (GF_COREAD.csv)"
            raise ValueError(msg)
        pass

    def mkdir(self, name):
        try:
            os.mkdir(name)
        except:
            if self.verbose:
                print("Note that directory %s already exists" % name)


class GDSC(GDSCBase):
    """Alias to ANOVA class with default settings to loop over all TCGA Tissues
    and comnpanies.

    Reads

    1. Nf Genomic feature files for different TCGA types.
    2. a unique DRUG DECODER file
    3. a unique IC50 file

    and perform the Nf analysis saving results in appropriate files.

    Then split the data for each companies.

    It also converts tissue names into TCGA names.

    .. todo:: How to get the GF files

    First, create all main analysis that include all drugs::

        gg = GDSC('IC50_v18.csv', 'DRUG_DECODE.txt',
            genomic_feature_pattern='GF*csv')
        # identifies all genomic features GF* that contains specific TCGA GF
        gg.analysis() # This may take a while (1 hour) depending on the number of drugs.

        # On v18, On an i7 core using 1 CPU this takes about 1 hour.30 minutes
        # The PANCAN data set is the largest and takes about 1 hour itself

    You should now have a directory called **tissue_packages** with about 20 directories for
    each TCGA GF file. Keep that in a safe place or you will have to restart
    the analysis

    Second, split those data just created for each specific proprietary
    compounds. For instance::

        gg.create_data_packages_for_companies(['AZ'])

    or for all in one go::

        gg.create_data_packages_for_companies()


    Third, create some summary pages::

        gg.create_summary_pages()

    The last step is fast (a few seconds) and create index.html in the
    tissue_package directory and each proprietary directory.


    """
    def __init__(self, ic50, drug_decode,
            genomic_feature_pattern="GF_*csv",
            main_directory="tissue_packages", verbose=True):
        """

        ic50 must be a filename (not IC50 instance) because it will be used for
        each genomic features file
        """
        super(GDSC, self).__init__(genomic_feature_pattern, verbose=verbose)
        assert isinstance(ic50, str)
        self.ic50_filename = ic50
        self.dd_filename = drug_decode
        self.main_directory = main_directory

        self.settings = ANOVASettings()
        self.settings.animate = False
        self.drug_decode = DrugDecode(drug_decode)

        print("Those settings will be used (check FDR_threshold)")
        print(self.settings)

        # figure out the cancer types:
        self.results = {}

        self.company_directory = "company_packages"

        # quick test on 15 features
        self.test = False

    def analyse(self, onweb=False, multicore=None):
        """Launch ANOVA analysis and creating data package for each tissue.

        :param bool onweb: By default, reports are created
            but HTML pages not shown. Set to True if you wish to open
            the HTML pages.
        :param multicore: number of cpu to use (1 by default)

        """
        self.mkdir(self.main_directory)
        # First analyse all TCGA cases + PANCAN once for all and
        # store all the results in a dictionary.
        self._analyse_all(onweb=onweb, multicore=multicore)

    def _analyse_all(self, onweb, multicore=None):
        for gf_filename in sorted(self.gf_filenames):
            tcga = gf_filename.split("_")[1].split('.')[0]
            print(purple('======================== Analysing %s data' % tcga))

            self.mkdir(self.main_directory + os.sep + tcga)
            # Computes the ANOVA
            try:
                self.ic50 = IC50(self.ic50_filename)
            except:
                print("Clustering IC50 (v18 released data ?)")
                self.ic50 = IC50Cluster(self.ic50_filename, verbose=False)
            an = ANOVA(self.ic50, gf_filename, self.drug_decode,
                verbose=False)

            if self.test is True:
                an.features.df = an.features.df[an.features.df.columns[0:15]]

            self.an = an
            an.settings = ANOVASettings(**self.settings)
            an.settings.analysis_type = tcga
            an.init() # This reset the directory

            results = an.anova_all(multicore=multicore)
            an.settings.directory = self.main_directory + os.sep + tcga
            # Store the results
            self.results[tcga] = results

            print('Analysing %s data and creating images' % tcga)
            self.report = ANOVAReport(an)
            self.report.settings.savefig = True

            self.report.create_html_pages(onweb=onweb)

    def create_data_packages_for_companies(self, companies=None):
        ##########################################################
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
        #                                                        #
        # DRUG_DECODE and IC50 inputs must be filtered to keep   #
        # only WEBRELEASE=Y and owner                            #
        #                                                        #
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
        ##########################################################

        # companies must be just one name (one string) or a list of strings
        # By default, takes all companies found in DrugDecode
        if isinstance(companies, str):
            companies = [companies]

        if companies is None:
            companies = self.companies

        if len(companies) == 0:
            raise ValueError("Could not find any companies in the DrugDecode file")

        # The main directory
        self.mkdir(self.company_directory)

        # Loop over all companies, retrieving information built
        # in analyse() method, selecting for each TCGA all information
        # for that company only (and public drugs)
        Ncomp = len(companies)
        for ii, company in enumerate(companies):
            print(purple("\n=========== Analysing company %s out of %s (%s)" %
                    (ii+1, Ncomp, company)))
            self.mkdir(self.company_directory + os.sep + company)

            # Handle each TCGA case separately
            for gf_filename in sorted(self.gf_filenames):
                tcga = gf_filename.split("_")[1].split('.')[0]
                print(brown("  ------- building TCGA %s sub directory" % tcga))

                # Read the results previously computed either
                try:
                    results_df = self.results[tcga].df.copy()
                except:
                    results_path = "%s/%s/OUTPUT/results.csv" % (self.main_directory, tcga)
                    results_df = ANOVAResults(results_path)

                # MAke sure the results are formatted correctly
                results = ANOVAResults(results_df)

                # Get the DrugDecode information for that company only
                drug_decode_company = self.drug_decode.df.query(
                        "WEBRELEASE=='Y' or OWNED_BY=='%s'" % company)

                # Transform into a proper DrugDecode class for safety
                drug_decode_company = DrugDecode(drug_decode_company)

                # Filter the results to keep only public drugs and that
                # company. Make sure this is integers
                results.df["DRUG_ID"] = results.df["DRUG_ID"].astype(int)

                mask = [True if x in drug_decode_company.df.index else False
                        for x in results.df.DRUG_ID]

                results.df = results.df.ix[mask]

                # We read the IC50 again
                try:
                    self.ic50 = IC50(self.ic50_filename)
                except:
                    self.ic50 = IC50Cluster(self.ic50_filename, verbose=False)

                # And create an ANOVA instance. This is not to do the analyse
                # again but to hold various information
                an = ANOVA(self.ic50, gf_filename, drug_decode_company,
                    verbose=False)

                def drug_to_keep(drug):
                    to_keep = drug in drug_decode_company.df.index
                    return to_keep
                an.ic50.df = an.ic50.df.select(drug_to_keep, axis=1)

                an.settings = ANOVASettings(**self.settings)
                an.init()
                an.settings.directory = self.company_directory + os.sep + company + os.sep + tcga
                an.settings.analysis_type = tcga

                # Now we create the report
                self.report = ANOVAReport(an, results, verbose=self.verbose)
                self.report.settings.analysis_type = tcga
                self.report.create_html_main(False)
                self.report.create_html_manova(False)
                self.report.create_html_features()
                self.report.create_html_associations()

                # Volcano plot for a given drug (ODAF)
                # do not need to be recomputed. If already done,
                # we just need to copy it (faster)
                # Although, names of the brand may not appear correctly
                from easydev import shellcmd, Progress
                drug_ids = results.df.DRUG_ID.unique()
                print("Copying drug HTML files")
                pb = Progress(len(drug_ids))
                for i, drug_id in enumerate(drug_ids):
                    # copy the HTML
                    filename = "drug_%s.html" % drug_id
                    source = "%s%s%s%sassociations%s" % (self.main_directory,
                        os.sep, tcga, os.sep , os.sep)
                    dest = os.sep.join([self.company_directory, company, tcga, "associations"])
                    dest += os.sep
                    cmd = "cp %s%s %s" % (source, filename, dest )
                    shellcmd(cmd, verbose=False)
                    #copy the images

                    if self.settings.animate:
                        pb.animate(i+1)
                if self.settings.animate:print("\n")

    def _get_tcga(self):
        return [x.split("_")[1].split(".")[0] for x in self.gf_filenames]
    tcga = property(_get_tcga)

    def _get_companies(self):
        return [x for x in self.drug_decode.companies if x != 'Commercial']
    companies = property(_get_companies)

    def create_summary_pages(self):
        """

        Once the main analyis is done (:meth:`analyse`), and the company
        packages have been created (:meth:`create_data_packages_for_companies`),
        you can run this method that will creade a summary HTML page
        (index.html) for the tissue, and a similar summary HTML page for the
        tissues of each company. Finally, an HTML summary page for the companies
        is also created. 

        The final tree direcorty looks like::


            |-- index.html
            |-- company_packages
            |   |-- Company1
            |   |   |-- Tissue1
            |   |   |-- Tissue2
            |   |   |-- index.html
            |   |-- Company2
            |   |   |-- Tissue1
            |   |   |-- Tissue2
            |   |   |-- index.html
            |-- tissue_packages
            |   |-- Tissue1
            |   |-- Tissue2
            |   |-- index.html




        """

        # First for the main directory (tissue_packages):
        print(purple("Creating summary index.html for the tissues"))
        self._create_summary_pages(self.main_directory, verbose=False)

        # Then for each companies:
        print(purple("Creating summary index.html for each company"))
        pb = Progress(len(self.companies))
        for i, company in enumerate(self.companies):
            try:
                self._create_summary_pages(self.company_directory + os.sep +
                    company, verbose=False, company=company)
            except Exception as err:
                print(red("Issue with %s. Continue with other companies" % company))
                print(err)
            pb.animate(i+1)

        # Finally, an index towards each company
        self._create_main_index()

    def _create_main_index(self):
        # We could also add a column with number of association ?
        companies = self.companies[:]
        df = pd.DataFrame({"Company": companies})
        html_page = ReportMAIN(directory=".",
                filename='index.html',
                template_filename='main_summary.html' )
        html_table = HTMLTable(df)
        html_table.add_href('Company', newtab=True, url="company_packages/",
                suffix="/index.html")
        html_page.jinja['data_table'] =  html_table.to_html(collapse_table=False)
        html_page.jinja['analysis_domain'] =  "All companies / All "
        html_page.jinja['tissue_directory'] = self.main_directory
        html_page.write()

    def _create_summary_pages(self, main_directory, verbose=True, 
            company=None):
        # Read all directories in tissue_packages

        directories = glob.glob(main_directory + os.sep + '*')

        summary = []
        for directory in sorted(directories):
            tcga = directory.split(os.sep)[-1]
            if tcga not in self.tcga:
                continue
            if verbose:
                print(directory, tcga)
            # number of hits
            path = directory + os.sep + 'OUTPUT' + os.sep
            try:
                hits = pd.read_csv(path + 'drugs_summary.csv', sep=',')
            except:
                summary.append([tcga] + [None] * 5)
                continue
            total_hits = hits.total.sum()

            drug_involved = hits['Unnamed: 0'].unique()

            results = ANOVAResults(path + 'results.csv')
            if len(results) > 0:
                drug_ids = results.df.DRUG_ID.unique()
            else:
                drug_ids = []

            # where to find the DRUG DECODE file. Should
            # have been copied
            path = directory + os.sep + 'INPUT' + os.sep
            drug_decode = DrugDecode(path + 'DRUG_DECODE.csv')
            info = drug_decode.get_info()

            webrelease = drug_decode.df.ix[drug_involved].WEBRELEASE
            drug_inv_public = sum(webrelease == 'Y')
            drug_inv_prop = sum(webrelease != 'Y')

            summary.append([tcga, total_hits,
                drug_inv_prop, info['N_prop'],
                drug_inv_public, info['N_public']])
        df = pd.DataFrame(summary)

        df.columns = ['Analysis name', 'Number of hits',
            'Number of involved proprietary compounds', 'out of',
            'Number of involved public', 'out of']

        try:
            df.sort_values(by="Number of hits", ascending=False, inplace=True)
        except:
            df.sort("Number of hits", ascending=False, inplace=True)

        output_dir = main_directory + os.sep + '..' + os.sep
        output_file = output_dir + os.sep + 'index.html'
        self.html_page = ReportMAIN(directory=main_directory,
                filename='index.html',
                template_filename='datapack_summary.html' )

        # Let us use our HTMLTable to add the HTML references
        self.html_table = HTMLTable(df)
        self.html_table.add_href('Analysis name', newtab=True, url=None,
                suffix='/index.html')
        self.html_table.add_bgcolor('Number of hits')

        self.html_page.jinja['data_table'] =  self.html_table.to_html(
                collapse_table=False)
        if company:
            self.html_page.jinja["collaborator"] = company

        self.html_page.write()

        return df





