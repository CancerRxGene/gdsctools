import os

from gdsctools.report import ReportMAIN, HTMLTable
from gdsctools.tissues import TCGA
from gdsctools.readers import DrugDecoder


import pandas as pd


class DataPackagesSummary(object):
    """Scan a set of data packages and create a summary report.


    In a given directory, identifies data packages and create a summary page

    """

    def __init__(self, collaborator=None, directories=None,
            directory='.'):
        self.directory = directory
        self.filename = 'index.html'
        self.collaborator =  collaborator

        if directories is None:
            self.directories = sorted(TCGA().keys())
        else:
            self.directories = directories

        data_table = []
        for directory in self.directories:
            if os.path.exists(directory) is False:
                print(directory +' not found. skipping...')
                continue
            else:
                print("Scanning " + directory)

            
            datapack = DataPackage(directory)

            data_table.append([directory,
                datapack.get_number_of_hits(),
                datapack.get_number_prop_involved(), 
                datapack.get_proprietary_compounds(),
                datapack.get_number_public_involved(), 
                datapack.get_public_compounds()])

        df = pd.DataFrame(data_table,
                columns=['Analysis name', 'Number of hits',
                    'N. proprietary compounds', 'out of',
                    'N. Public compounds', 'out of'])
        html_table = HTMLTable(df)
        html_table.add_href('Analysis name', suffix='/index.html', 
                newtab=True)
        self.html = html_table.to_html()

    def create_report(self, onweb=True):
        report = ReportMAIN(directory=self.directory,
                    template_filename="datapack_summary.html",
                    mode='summary')

        report.jinja['data_table'] = self.html
        report.jinja['collaborator'] = self.collaborator
        report.create_report(onweb=onweb)


class DataPackage(object):

    def __init__(self, directory):
        self.directory = directory
        self.dd = DrugDecoder(os.sep.join([directory,'INPUT',
            'DRUG_DECODE.csv']))
        filename = os.sep.join([self.directory, 'OUTPUT',
            'drugs_summary.csv'])
        self.associations = pd.read_csv(filename)
        if len(self.associations)>0:
            cols = list(self.associations.columns)
            self.associations.columns = ['DRUG_ID'] + cols[1:]

    def get_number_prop_involved(self):
        if len(self.associations) == 0:
            return 0
        else:
            data = [self.dd.get_public(x) for x in self.associations.DRUG_ID]
            return data.count('N')

    def get_number_public_involved(self):
        if len(self.associations) == 0:
            return 0
        else:
            data = [self.dd.get_public(x) for x in self.associations.DRUG_ID]
            return data.count('Y')

    def get_number_of_hits(self):
        if len(self.associations) == 0:
            return 0
        else:
            return self.associations.sum().total

    def get_proprietary_compounds(self):
        return self.dd.get_info()['N_prop']

    def get_public_compounds(self):
        return self.dd.get_info()['N_public']

    def get_total_compounds(self):
        return self.dd.get_info()['N']

