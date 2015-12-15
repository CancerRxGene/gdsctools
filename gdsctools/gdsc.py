import glob
import os

from gdsctools.anova import ANOVA
from gdsctools.readers import IC50
from gdsctools.readers import DrugDecode
from gdsctools.tools import extract_drug_identifiers
from gdsctools.anova_results import ANOVAResults
from gdsctools.anova_report import ANOVAReport
from gdsctools.settings import ANOVASettings

import pandas as pd


class IC50Cluster(IC50):
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
        to_cluster = info[info.ratio < 10].DRUG_ID.values
        return list(to_cluster)
    to_cluster = property(_get_to_cluster)

    def _get_mapping(self):
        from collections import defaultdict
        mapping = defaultdict(list)

        drug_ids = extract_drug_identifiers(self.df.columns)
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





class GDSC(object):
    """Alias to ANOVA class with default settings

    Reads

    1. Nf  Genomic feature files for different TCGA types.
    2. a unique DRUG DECODER file
    3. a unique IC50 file

    and perform the Nf analysis saving results in appropriate files.

    Then split the data for each comapnies.

    It also converts tissue names into TCGA names.
    """

    def __init__(self, ic50, drug_decode, genomic_feature_pattern="GF_*csv",
            mode='standard'):
        self.verbose = True

        self.ic50_filename = ic50
        self.dd_filename = drug_decode

        if mode == 'v18':
            self.ic50 = IC50Cluster(ic50)
        else:
            self.ic50 = IC50(ic50)

        self.drug_decode = DrugDecode(drug_decode)
        self.settings = ANOVASettings()
        if mode == 'v18':
            self.settings.FDR_threshold = 35

        print("Those settings will be used")
        print(self.settings)

        self.gf_filenames = glob.glob(genomic_feature_pattern)
        if len(self.gf_filenames) == 0:
            msg = "NO Genomic feature input files found. We expect files " +\
                    "with this pattern: GF_<TCGA>.csv e.g., (GF_COREAD.csv)"
            raise ValueError(msg)
        # figure out the cancer types:
        self.results = {}

    def run(self):
        self.mkdir('ALL')
        # First analyse all case of TCGA + PANCAN once for all and
        # store all results in a dictionary.
        self._analyse_all()
        self._create_reports( )
        # Now that all data has been analysed, we can split
        # for each company. Note that some pictures will be created
        # but not the drug-related ones, which will emcompass all
        # drugs across all companies. 
        #self._create_data_packages_for_companies


    def _analyse_all(self):
        for gf_filename in sorted(self.gf_filenames):
            tcga = gf_filename.split("_")[1].split('.')[0]
            print('================================ Analysing %s data' % tcga)

            self.mkdir('ALL' + os.sep + tcga)

            low_memory = True
            if tcga != 'PANCAN':
                low_memory = False

            an = ANOVA(self.ic50_filename, gf_filename, self.drug_decode,
                    low_memory = low_memory)
            self.an = an
            an.settings = ANOVASettings(**self.settings)
            an.init() # reset the analysis_type automatically

            # this line may not be required.
            an.settings.low_memory = low_memory
            print(an)
            results = an.anova_all()
            self.results[tcga] = results

    def _create_reports(self):
        from gdsctools import ANOVAReport

        for gf_filename in sorted(self.gf_filenames):
            tcga = gf_filename.split("_")[1].split('.')[0]
            print('Analysing %s data' % tcga)

            an = ANOVA(self.ic50_filename, gf_filename, self.drug_decode,
                       low_memory=True)
            # FIXME use os.sep
            an.settings = ANOVASettings(**self.settings)
            an.init()
            #an.settings.analysis_type = tcga
            an.settings.directory = 'ALL/' + tcga
            self.report = ANOVAReport(an, self.results[tcga])
            self.report.settings.analysis_type = tcga
            self.report.create_html_pages()

    def _create_data_packages_for_companies(self):
        Ncomp = len(self.companies)
        for ii, company in enumerate(self.companies):
            print("\n\n========= Analysing company %s out of %s (%s)" %
                    (ii+1, Ncomp, company))
            self.mkdir(company)
            for gf_filename in self.gf_filenames:
                tcga = gf_filename.split("_")[1].split('.')[0]
                print("---------------- for TCGA %s" % tcga)

                # Read the results previously computed
                results_df = self.results[tcga].df.copy()
                results = ANOVAResults(results_df)

                # Get a DrugDecode for that company
                drug_decode_company = self.drug_decode.df.query(
                        "WEBRELEASE=='Y' or OWNED_BY=='%s'" % company)
                # Transform into a proper DrugDecode class for safety
                drug_decode_company = DrugDecode(drug_decode_company)

                # filter results using the new drug decode
                drug_ids_in_results = extract_drug_identifiers(
                        results.df.DRUG_ID)

                mask = [True if x in drug_decode_company.df.index else False
                        for x in drug_ids_in_results]

                results.df = results.df.ix[mask]

                # just to create an instance with the subset of drug_decode
                # and correct settings
                an = ANOVA(self.ic50_filename, gf_filename,
                        drug_decode_company, low_memory=True)
                an.settings = ANOVASettings(**self.settings)
                an.settings.directory = company + os.sep + tcga
                an.settings.analysis_type = tcga

                r = ANOVAReport(an, results)
                r.settings.analysis_type = tcga
                r.create_html_main(False)
                r.create_html_manova(False)
                r.create_html_features()
                r.create_html_associations()

                # FIXME. Most of the HTML and images may have already been
                # created in the _create_reports when using all the data
                # For now, we just recreate everything but could be speed up
                # significantly.
                from easydev import shellcmd
                print("Copying drug files")
                from easydev import Progress
                drug_ids = results.df.DRUG_ID.unique()
                pb = Progress(len(drug_ids))
                for i, drug_id in enumerate(drug_ids):
                    # copy the HTML
                    filename = "%s.html" % drug_id
                    source = "ALL%s%s%s" % (os.sep, tcga, os.sep)
                    dest = "%s%s%s%s" % (company, os.sep, tcga, os.sep )
                    cmd = "cp %s%s %s" % (source, filename, dest )
                    shellcmd(cmd, verbose=False)
                    #copy the images
                    filename = "volcano_%s.*" % drug_id
                    source = "ALL%s%s%simages%s" % (os.sep, tcga,
                            os.sep, os.sep)
                    dest = "%s%s%s%simages%s" % (company, os.sep,
                            tcga, os.sep , os.sep)
                    cmd = "cp %s%s %s" % (source, filename, dest )
                    shellcmd(cmd, verbose=False)
                    pb.animate(i+1)

    def _get_tcga(self):
        return [x.split("_")[1].split(".")[0] for x in self.gf_filenames]
    tcga = property(_get_tcga)

    def mkdir(self, name):
        try:
            os.mkdir(name)
        except:
            if self.verbose:
                print("directory %s already exists" % name)

    def _get_companies(self):
        return [x for x in self.drug_decode.companies if x != 'Commercial']
    companies = property(_get_companies)



    def create_summary_pages(self):
        # Read in ALL all directories
        import glob
        directories = glob.glob('ALL' + os.sep + '*')

        summary = []
        for directory in sorted(directories):
            tcga = directory.split(os.sep)[1]
            print(tcga)

            # number of hits
            path = directory + os.sep + 'OUTPUT' + os.sep
            try:
                hits = pd.read_csv(path + 'drugs_summary.csv', sep=',')
            except:
                summary.append([tcga] + ['?'] * 5)
                continue
            total_hits = hits.total.sum()

            drug_involved = extract_drug_identifiers(
                    hits['Unnamed: 0'].unique())


            results = ANOVAResults(path + 'results.csv')
            if len(results)>0:
                drug_ids = extract_drug_identifiers(results.df.DRUG_ID.unique())
            else:
                drug_ids = []

            path = directory + os.sep + 'INPUT' + os.sep
            drug_decode = DrugDecode(path + 'DRUG_DECODE.csv')
            info = drug_decode.get_info()

            drug_inv_public = sum(drug_decode.df.ix[drug_involved].WEBRELEASE == 'Y')
            drug_inv_prop = sum(drug_decode.df.ix[drug_involved].WEBRELEASE != 'Y')

            summary.append([tcga, total_hits,
                drug_inv_prop, info['N_prop'], 
                drug_inv_public, info['N_public']])
        print summary
        df = pd.DataFrame(summary)
        df.columns = ['Analysis name', 'Number of hits', 
        'Number of involved proprietary compounds', 'out of',
        'Number of involved public', 'out of']
        #html = df.to_html()
        
        return df
            








