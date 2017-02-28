"""Download GDSC drug data sets (IC50s 

Module for downloading IC50 data from the gdsc1000 paper.
Using one line, the script also formats the data so that it can be directly imported into the ANOVA script and also creates a drug decoder.
E.D.Chen 2016-06-06



"""
import urllib.request
import pandas as pd
import os


class DownloadDrugIC50s(object):
    """Download IC50 data from GDSC website

    ::

        from gdsctools.gdsc1000 import DownloadDrug
        fdd = fetch_drug_data.DownloadDrugIC50s()
        fdd.download_data()

    data is formatted automatically and saved as './data/gdsc_1000_data/ic50.csv'
    if wanted, drug data frame can be filtered by calling:

        fdd.filter_drugs( [ 221, 211, 152, ... ] )

    Where drug IDs are given as a list of integers. The resulting data frame is 
    then stored as './data/gdsc_1000_data/ic50_filtered.csv'

    """

    def __init__(self):
        self.url_base = "http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/suppData/"
        self.data_folder_name = "./data/gdsc_1000_data/"

    def fetch_ic50s(self):

        if os.path.exists(self.data_folder_name):
            pass
        else:
            from easydev import mkdirs
            mkdirs(self.data_folder_name)

        urllib.request.urlretrieve(self.url_base + "TableS4A.xlsx", 
            self.data_folder_name + "ic50s.xlsx" )

        self._format_data()

    def _format_data(self):
        xls = pd.ExcelFile( self.data_folder_name + "ic50s.xlsx" )
        df = xls.parse( header = 4 )
        df.columns.values[ 0 ] = 'COSMIC_ID'
        df = df.drop(0)
        df = df.drop( 'Sample Names', axis = 1 )
        self.drug_matrix = df
        df.to_csv( self.data_folder_name + 'ic50.csv', index = False )

        drug_dict = xls.parse(header=4)[0:1]
        drug_dict.drop( drug_dict.columns[[0, 1]], axis = 1, inplace = True )
        drug_dict = drug_dict.T.reset_index()
        drug_dict.columns = [ 'DRUG_ID', 'DRUG_NAME' ]
        drug_dict[ 'DRUG_TARGET' ] = 'Unknown'
        drug_dict.to_csv( self.data_folder_name + 'drug_decoder.csv', index = False )

    def filter_drugs(self, drug_list):
        drug_list = [ x for x in drug_list if x in self.drug_matrix.columns.values ]

        self.filter_df = self.drug_matrix[ [ 'COSMIC_ID' ] + drug_list ]
        self.filter_df.to_csv( self.data_folder_name + 'ic50_filtered.csv', index = False )
        print( "Following Drug IDs now saved as 'ic50_filtered.csv': \n" + str(drug_list) )
