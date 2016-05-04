import glob
import os
import easydev
import shutil

from gdsctools.anova import ANOVA
from gdsctools.readers import IC50
from gdsctools.readers import DrugDecode
from gdsctools.tools import get_drug_id
from gdsctools.anova_results import ANOVAResults
from gdsctools.anova_report import ANOVAReport
from gdsctools.settings import ANOVASettings
from gdsctools.anova_report import ReportMAIN

import pandas as pd


class IC50Cluster(IC50):
    """

    Used in v18 only

    """
    def __init__(self, ic50, ratio_threshold=10, verbose=True):
        """

        """
        super(IC50Cluster, self).__init__(ic50)
        self.verbose = verbose
        self.ratio_threshold = ratio_threshold

        if self.verbose:
            print(self)
        self.cluster()
        if self.verbose:
            print(self)

    def _get_to_cluster(self):
        info = self.info()
        if len(info)> 0:
            to_cluster = info[info.ratio < 10].DRUG_ID.values
            return list(to_cluster)
        else:
            return []
    to_cluster = property(_get_to_cluster)

    def _get_mapping(self):
        from collections import defaultdict
        mapping = defaultdict(list)

        drug_ids = get_drug_id(self.df.columns)
        for drug_id, colname in zip(drug_ids, self.df.columns):
            mapping[drug_id].append(colname)
        return mapping

    def _get_duplicated(self):
        mapping = self._get_mapping()
        duplicated = [key for key in mapping.keys() if len(mapping[key])>1]
        return duplicated
    duplicated = property(_get_duplicated)

    def info(self):
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

            common = sum(df.count(axis=1)>=2)
            result = [drug_id] + individuals + [total, common,
                    100*common/float(total)]
            results.append(result)

        df = pd.DataFrame(results)
        if len(df):
            df.columns = ['DRUG_ID'] + [str(x) for x in range(1, max_ids+1)] +\
                ['total', 'common', 'ratio']
            df.sort_values(by='ratio', ascending=False, inplace=True)
        return df

    def cluster(self):
        # get list of drug identifiers to cluster

        to_cluster = self.to_cluster
        mapping = self._get_mapping()

        if self.verbose:
            print('Found  %s non unique drug identifiers ' %
                    len(self.duplicated))
            print('Clustering %s of them.\n' % len(self.to_cluster))
        if len(self.to_cluster) == 0:
            return

        for identifier in to_cluster:
            drug_names = mapping[identifier]
            new_drug_name = str(identifier) + "_" + \
                    "_".join([x.split("_",1)[1] for x in  drug_names])
            # add new column with new name and mean of the columns with same
            # drug id
            self.df[new_drug_name] = self.df[drug_names].mean(axis=1)
            # Remove the individual columns
            self.df.drop(drug_names, axis=1, inplace=True)


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
                print("directory %s already exists" % name)


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
        gg.run() # This will take hours depending on the number of drugs.


        # On v18, On an i7 core using 1 CPU this tqkes about 1 hour.30 minutes
        # THe PANCAN data set is the largest and takes about 1 hour itself

    You should now have a directory called **tissue_packages** with about 20 directories for
    each TCGA GF file. Keep that in a safe place or you will have to restart
    the analysis

    Second, split those data just created for each specific proprietary 
    compounds. For instance::

        gg.create_data_packages_for_companies(['AZ'])
    
    or for all in one go::

        gg.create_data_packages_for_companies()


    Third, create some summary pages::

        from gdsctools.gdsc import GDSCDirectorySummary()
        gs = GDSCDirectorySummary()
        gs.create_summary_pages('tissue_packages')
        for company in gg.companies:
            gs.create_summary_pages(company)

    The last step is fast but the whole process of analyse and image 
    creation is very long. 



    """
    def __init__(self, ic50, drug_decode,
            genomic_feature_pattern="GF_*csv",
            mode='standard', main_directory="tissue_packages"):
        super(GDSC, self).__init__(genomic_feature_pattern, verbose=True)
        self.debug = False
        self.ic50_filename = ic50
        self.dd_filename = drug_decode
        self.main_directory = main_directory

        if mode == 'v18':
            self.ic50 = IC50Cluster(ic50)
        else:
            self.ic50 = IC50(ic50)

        self.drug_decode = DrugDecode(drug_decode)
        self.settings = ANOVASettings()
        self.settings.low_memory = True

        if mode == 'v18':
            self.settings.FDR_threshold = 35

        print("Those settings will be used (note that low_memory is set "
              "to True and check the value of FDR_threshold (set to 35 in v18)")
        print(self.settings)

        # figure out the cancer types:
        self.results = {}

    def run(self):
        self.mkdir(self.main_directory)
        # First analyse all case of TCGA + PANCAN once for all and
        # store all results in a dictionary.
        self._analyse_all()

    def _analyse_all(self):
        for gf_filename in sorted(self.gf_filenames):
            tcga = gf_filename.split("_")[1].split('.')[0]
            print('================================ Analysing %s data' % tcga)

            self.mkdir(self.main_directory + os.sep + tcga)

            # Computes the ANOVA
            an = ANOVA(self.ic50_filename, gf_filename, self.drug_decode)
            self.an = an
            an.settings = ANOVASettings(**self.settings)
            an.init() # reset the analysis_type automatically
            results = an.anova_all()

            # Store the results
            self.results[tcga] = results

            print('Analysing %s data and creating images' % tcga)
            self.report = ANOVAReport(an, self.results[tcga])
            self.report.settings.savefig = True
            self.report.settings.directory = self.main_directory + os.sep + tcga
            self.report.settings.analysis_type = tcga
            self.report.create_html_pages()

    def create_data_packages_for_companies(self, companies=None):
        ##########################################################
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
        #                                                        #
        # DRUG_DECODE and IC50 inputs must be filtered to keep   #
        # only WEBRELEASE=Y and owner                            #
        #                                                        #
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
        ##########################################################
        if isinstance(companies, str):
            companies = [companies]

        if companies is None:
            companies = self.companies

        Ncomp = len(companies)
        for ii, company in enumerate(companies):
            print("\n\n========= Analysing company %s out of %s (%s)" %
                    (ii+1, Ncomp, company))
            self.mkdir(company)
            for gf_filename in sorted(self.gf_filenames):
                tcga = gf_filename.split("_")[1].split('.')[0]
                print("---------------- for TCGA %s" % tcga)

                # Read the results previously computed
                try:
                    results_df = self.results[tcga].df.copy()
                except:
                    results_path = "%s/%s/OUTPUT/results.csv" % (self.main_directory, tcga)
                    print("Downloading results from %s" % results_path)
                    results_df = ANOVAResults(results_path)

                results = ANOVAResults(results_df)

                # Get a DrugDecode for that company
                drug_decode_company = self.drug_decode.df.query(
                        "WEBRELEASE=='Y' or OWNED_BY=='%s'" % company)
                # Transform into a proper DrugDecode class for safety
                drug_decode_company = DrugDecode(drug_decode_company)

                # filter results using the new drug decode
                drug_ids_in_results = get_drug_id(results.df.DRUG_ID)

                mask = [True if x in drug_decode_company.df.index else False
                        for x in drug_ids_in_results]

                results.df = results.df.ix[mask]

                # Just to create an instance with the subset of drug_decode
                # and correct settings. This is also used to store
                # the entire input data set. So, we must remove all drugs
                # not relevant for the analysis of this company
                an = ANOVA(self.ic50_filename, gf_filename, drug_decode_company)

                def drug_to_keep(drug):
                    to_keep = get_drug_id(drug) in drug_decode_company.df.index
                    return to_keep
                an.ic50.df = an.ic50.df.select(drug_to_keep, axis=1) 

                an.settings = ANOVASettings(**self.settings)
                an.init()
                an.settings.directory = company + os.sep + tcga
                an.settings.analysis_type = tcga
                self.report = ANOVAReport(an, results)
                self.report.settings.analysis_type = tcga
                self.report.create_html_main(False)
                self.report.create_html_manova(False)

                if self.debug is False:
                    self.report.create_html_features()
                    self.report.create_html_associations()

                    # For now, we just copy all DRUG images from 
                    # the analysis made in tissue_packages 
                    from easydev import shellcmd, Progress
                    print("\nCopying drug files")
                    drug_ids = results.df.DRUG_ID.unique()
                    pb = Progress(len(drug_ids))
                    for i, drug_id in enumerate(drug_ids):
                        # copy the HTML
                        filename = "%s.html" % drug_id
                        source = "%s%s%s%s" % (self.main_directory, os.sep, tcga, os.sep)
                        dest = "%s%s%s%s" % (company, os.sep, tcga, os.sep )
                        cmd = "cp %s%s %s" % (source, filename, dest )
                        shellcmd(cmd, verbose=False)
                        #copy the images
                        filename = "volcano_%s.*" % drug_id
                        source = "%s%s%s%simages%s" % (self.main_directory, os.sep, tcga,
                                os.sep, os.sep)
                        dest = "%s%s%s%simages%s" % (company, os.sep,
                                tcga, os.sep , os.sep)
                        cmd = "cp %s%s %s" % (source, filename, dest )
                        shellcmd(cmd, verbose=False)
                        pb.animate(i+1)

    def _get_tcga(self):
        return [x.split("_")[1].split(".")[0] for x in self.gf_filenames]
    tcga = property(_get_tcga)

    def _get_companies(self):
        return [x for x in self.drug_decode.companies if x != 'Commercial']
    companies = property(_get_companies)

    def create_summary_pages(self):
        # Read in ALL all directories
        main_directory = self.main_directory

        # create directories and copy relevant files
        self.mkdir(main_directory + os.sep + 'images')
        self.mkdir(main_directory + os.sep + 'css')
        self.mkdir(main_directory + os.sep + 'js')
        from gdsctools import gdsctools_data
        for filename in ['gdsc.css', 'github-gist.css']:
            target = os.sep.join([main_directory, 'css', filename ])
            if os.path.isfile(target) is False:
                filename = gdsctools_data(filename)
                shutil.copy(filename, target)
        for filename in ['highlight.pack.js']:
            target = os.sep.join([main_directory, 'js', filename ])
            if os.path.isfile(target) is False:
                filename = gdsctools_data(filename)
                shutil.copy(filename, target)
        for filename in ['EBI_logo.png', 'sanger-logo.png']:
            target = os.sep.join([main_directory, 'images', filename ])
            if os.path.isfile(target) is False:
                dire = 'data' + os.sep + 'images'
                filename = gdsctools_data("images" + os.sep +filename)
                shutil.copy(filename, target)
        directories = glob.glob(self.main_directory + os.sep + '*')
        directories = [x for x in directories if os.path.isdir(x)]

        summary = []
        for directory in sorted(directories):

            tcga = directory.split(os.sep)[1]
            if tcga in ['css', 'images']:
                continue

            # number of hits
            path = directory + os.sep + 'OUTPUT' + os.sep
            try:
                hits = pd.read_csv(path + 'drugs_summary.csv', sep=',')
            except:
                summary.append([tcga] + [None] * 5)
                continue
            total_hits = hits.total.sum()

            drug_involved = get_drug_id(hits['Unnamed: 0'].unique())

            results = ANOVAResults(path + 'results.csv')
            if len(results)>0:
                drug_ids = get_drug_id(results.df.DRUG_ID.unique())
            else:
                drug_ids = []

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


        # FIXME include css and images of logo
        # FIXME save in the proper directory
        output_dir = main_directory + os.sep + '..' + os.sep
        output_file = output_dir + os.sep + 'index.html'
        self.html_page = ReportMAIN(directory=self.main_directory, filename='index.html',
                template_filename='datapack_summary.html' )

        # Let us use our HTMLTable to add the HTML references
        from gdsctools.report import HTMLTable
        self.html_table = HTMLTable(df)
        self.html_table.add_href('Analysis name', newtab=True, url=None,
                suffix='/index.html')
        #html_table.add_bgcolor('Number of hits')

        self.html_page.jinja['data_table'] = self.html_table.to_html()
        self.html_page.write()

        return df

    def load_results(self):
        """Find the files results.csv in all TCGA directories"""
        for tcga in self.tcga:
            print(tcga)
            self.results[tcga] = ANOVAResults(self.main_directory + os.sep + tcga + os.sep +
                    'OUTPUT' + os.sep + 'results.csv')




