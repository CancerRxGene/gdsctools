""" This file is a dump for code that is designed to fetch the GDSC1000 data sets and to convert them into the appropriate format.
E.D.Chen 2016-06-01
"""

import urllib.request
import pandas as pd
import os

class FetchGDSCData(object):
    """

    These commands will download the CNA, Methylation, Variant and CellLines
    data from GDSC website::

        fd = FetchGDSCData()
        fd.download()

    """
    def __init__(self):
        self.url_base ="http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/suppData/"
        self.data_folder_name = "./data/gdsc_1000_data/"

        self.methylation_path = self.data_folder_name + "methylation.xlsx"
        self.cna_path = self.data_folder_name + "cna.xlsx"
        self.variant_path = self.data_folder_name + "variant.xlsx"
        self.cell_dict_path = self.data_folder_name + "cell_line_dict.xlsx"

    def download(self):
        self._download_data()
        self._read_cell_line_data()
        self._read_cna_data()
        self._read_methylation_data()
        self._read_variant_data()

    # This block downloads the data locally.
    def _download_data(self):
        print("Downloading data...")

        methylation_url = self.url_base +  "TableS2J.xlsx"
        cna_url = self.url_base + "TableS2G.xlsx"
        variant_url = self.url_base + "TableS2C.xlsx"
        cell_line_url = self.url_base + "TableS1E.xlsx"

        if os.path.exists(self.data_folder_name):
            pass
        else:
            os.mkdir(self.data_folder_name)

        urllib.request.urlretrieve(methylation_url, self.methylation_path)
        urllib.request.urlretrieve(cna_url, self.cna_path)
        urllib.request.urlretrieve(variant_url, self.variant_path)
        urllib.request.urlretrieve(cell_line_url, self.cell_dict_path)

    def _read_cna_data(self):
        print("Processing CNA data.")
        # This block reads in and formats CNA data
        xls = pd.ExcelFile( self.cna_path )
        df = xls.parse(header=2)
        df.columns = ["BLANK", "CELL_LINE", "TISSUE_ID1", "TISSUE_ID2", "TISSUE_FACTOR", "ALTERATION_TYPE", "IDENTIFIER" ]
        df = df.drop(['BLANK', 'TISSUE_ID1', 'TISSUE_ID2'], axis = 1)
        df.IDENTIFIER = df.IDENTIFIER.str.split(' ').str[0]
        df.to_csv(self.data_folder_name + 'cna.csv', index=False)

    def _read_methylation_data(self):
        print("Processing methylation data.")
        # This block deals with the methylation data
        xls = pd.ExcelFile( self.methylation_path )
        df = xls.parse(header=2)
        df.columns = ["BLANK", "CELL_LINE", "TISSUE_ID1", "TISSUE_ID2", "TISSUE_FACTOR", "IDENTIFIER"]
        df = df.drop(["BLANK", "TISSUE_ID1", "TISSUE_ID2"], axis = 1)
        df.to_csv(self.data_folder_name + 'methylation.csv', index=False)

    def _read_variant_data(self):
        print("Processing variant data.")
        # This block deals with variant data
        xls = pd.ExcelFile( self.variant_path )
        df = xls.parse( skiprows=20 )
        df = df[ [ 'SAMPLE', 'COSMIC_ID', 'Cancer Type', 'Gene', 'Recurrence Filter' ] ]
        df = df.loc[ df[ 'Recurrence Filter'] == 'Yes' ]
        df.columns = [ 'CELL_LINE', 'COSMIC_ID', 'TISSUE_FACTOR', 'IDENTIFIER', 'RECURRENCE' ]
        df.to_csv(self.data_folder_name + 'variant.csv', index=False )

    def _read_cell_line_data(self):
        print("Processing cell line dictionary.")
        # This block imports the cell line dictionary
        xls = pd.ExcelFile( self.cell_dict_path )
        df = xls.parse( header = 2 )
        df.columns = [ "CELL_LINE", "COSMIC_ID", "WES", "CNA", "GEX", "METH", "DRUGDATA", "T_ID1", "T_ID2", "TISSUE_FACTOR", "MSI_FACTOR", "MEDIA_FACTOR", "GROWTH_FACTOR" ]
        df.index = range( len( df.index ) )
        df = df.drop(0)
        df = df[ [ 'CELL_LINE', 'COSMIC_ID', 'TISSUE_FACTOR', 'MSI_FACTOR', 'MEDIA_FACTOR', 'GROWTH_FACTOR' ] ]
        df.to_csv( self.data_folder_name + 'cell_line_dict.csv', index = False )

