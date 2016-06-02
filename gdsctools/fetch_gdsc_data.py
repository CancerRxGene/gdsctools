""" This file is a dump for code that is designed to fetch the GDSC1000 data sets and to convert them into the appropriate format. 
E.D.Chen 2016-06-01
"""


# This block downloads the data locally. 
# NOTE: in future, try to figure out how to define separate folder for download.
import urllib.request
import pandas as pd

class FetchGDSCData(object):
    """
    
    These commands will download the CNA, Methylation, Variant and CellLines 
    data from GDSC website::
    
        fd = FetchGDSCData()
        fd.download()
    
    """
    def __init__(self):
        self.url_base = "http://www.cancerrxgene.org/gdsc1000//Data/suppData/"

    def download(self):
        self._download_data()
        self._read_cell_line_data()
        self._read_cna_data()
        self._read_methylation_data()
        self._read_variant_data()
    
    def _download_data(self):
        methylation_url = self.url_base +  "TableS2J.xlsx"
        cna_url = self.url_base + "TableS2G.xlsx"
        variant_url = self.url_base + "TableS2C.xlsx"
        cell_line_url = self.url_base + "TableS1E.xlsx"

        urllib.request.urlretrieve(methylation_url, "methylation.xlsx")
        urllib.request.urlretrieve(cna_url, "cna.xlsx")
        urllib.request.urlretrieve(variant_url, "variant.xlsx")
        urllib.request.urlretrieve(cell_line_url, "cell_line_dict.xlsx")

    def _read_cna_data(self):
        # This block reads in and formats CNA data
        xls = pd.ExcelFile('cna.xlsx')
        df = xls.parse(header=2)
        df.columns = ["BLANK", "CELL_LINE", "TISSUE_ID1", "TISSUE_ID2", "TCGA", "ALTERATION_TYPE", "IDENTIFIER" ]
        df = df.drop(['BLANK', 'TISSUE_ID1', 'TISSUE_ID2'], axis = 1)
        df.to_csv('cna.csv', index=False)

    def _read_methylation_data(self):
        # This block deals with the methylation data
        xls = pd.ExcelFile('methylation.xlsx')
        df = xls.parse(header=2)
        df.columns = ["BLANK", "CELL_LINE", "TISSUE_ID1", "TISSUE_ID2", "TCGA", "IDENTIFIER"]
        df = df.drop(["BLANK", "TISSUE_ID1", "TISSUE_ID2"], axis = 1)
        df.to_csv('methylation.csv', index=False)

    def _read_variant_data(self):
        # This block deals with variant data
        xls = pd.ExcelFile('variant.xlsx')
        df = xls.parse( skiprows=20 )
        df = df[ [ 'SAMPLE', 'COSMIC_ID', 'Cancer Type', 'Gene', 'Recurrence Filter' ] ]
        df = df.loc[ df[ 'Recurrence Filter'] == 'Yes' ]
        df.columns = [ 'CELL_LINE', 'COSMIC_ID', 'TCGA', 'IDENTIFIER', 'RECURRENCE' ]
        df.to_csv( 'variant.csv', index=False )

    def _read_cell_line_data(self):
        # This block imports the cell line dictionary
        xls = pd.ExcelFile( "cell_line_dict.xlsx" )
        df = xls.parse( header = 2 )
        df.columns = [ "CELL_LINE", "COSMIC_ID", "WES", "CNA", "GEX", "METH", "DRUGDATA", "T_ID1", "T_ID2", "TISSUE_FACTOR", "MSI_FACTOR", "MEDIA_FACTOR", "GROWTH_FACTOR" ]
        df.index = range( len( df.index ) )
        df = df.drop(0)
        df = df[ [ 'CELL_LINE', 'COSMIC_ID', 'TISSUE_FACTOR', 'MSI_FACTOR', 'MEDIA_FACTOR', 'GROWTH_FACTOR' ] ]
        
        
        