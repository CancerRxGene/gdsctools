import os
import glob

from gdsctools.report import HTMLTable
from gdsctools.report import ReportMAIN
from gdsctools.tissues import TCGA
from gdsctools.readers import DrugDecode, get_drug_id
from gdsctools.anova_results import ANOVAResults

import pandas as pd
from gdsctools.gdsc import GDSCBase


class DataPackagesSummary(object):
    """Scan a set of data packages and create a summary report.


    In a given directory, identifies data packages and create a summary page

    What are the initial files ?
    gf = GenomicFeatures('ANOVA_input.txt')

    """

    def __init__(self, collaborator=None, directories=None,
            directory='.'):
        pass
        # Identify all companies
        # r.df.groupby("OWNED_BY").groups.keys()


class GDSCDirectorySummary(GDSCBase):
    """


    First, one need to create the GF file for each TCGA cancer



    """
    def __init__(self, genomic_feature_pattern="GF_*csv"):
        super(GDSCDirectorySummary, self).__init__(genomic_feature_pattern)

    def create_summary_pages(self, main_directory="tissue_packages"):
        # Read in ALL all directories

        # create directories and copy relevant files
        directories = glob.glob(main_directory + os.sep + '*')
        directories = [x for x in directories if os.path.isdir(x)]

        summary = []
        for directory in sorted(directories):

            tcga = directory.split(os.sep)[1]
            if tcga in ['css', 'images', 'INPUT', 'OUTPUT', 'code', 'js' ]:
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
            if len(results) > 0:
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
        output_dir = main_directory + os.sep
        output_file = output_dir + os.sep + 'index.html'
        self.html_page = ReportMAIN(directory=main_directory, 
                filename='index.html',
                template_filename='datapack_summary.html')

        # Let us use our HTMLTable to add the HTML references
        self.html_table = HTMLTable(df)
        self.html_table.add_href('Analysis name', newtab=True, url=None,
                suffix='/index.html')
        #html_table.add_bgcolor('Number of hits')

        self.html_page.jinja['data_table'] = self.html_table.to_html()
        self.html_page.jinja['collaborator'] = main_directory
        #self.html_page.jinja['analysis_domain'] = 'v18'
        self.html_page.write()

        return df


