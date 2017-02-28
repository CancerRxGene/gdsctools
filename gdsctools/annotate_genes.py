"""Gene annotations (gene_annotate module)

This module concatenates previously fetched GDSC1000 files into a
single data frame that can be converted into a BEM and fed into the
ANOVA analysis step of gdsctools. 

If needed, this script has the option to convert CNA and
methylation data from their ID values to their gene annotations.

E.D.Chen 2016-06-02

"""
import urllib.request as urlrequest
import pandas as pd

__all__ = ["MergeGDSCData"]

class MergeGDSCData(object):
    """

    This class depends on having previously run the fetch_gdsc_data class in the same folder::

        an = annotate_genes.MergeGDSCData()
        an.load_files() # loads previously fetched files into class
        an.annotate_all() # Optional. adds gene-level annotation to CNA and methylation data
        an.merge_all() # Concatantes variant_df, cna_df and methyl_df into a single data frame (complete_df)

    """

    def __init__(self):
        self.url_base = "http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/suppData/"
        self.data_folder_name = "./data/gdsc_1000_data/"
        self.annotate = False

    def load_files( self, annotate = None ):
        if annotate == None:
            annotate = self.annotate
        self.variant_df = pd.read_csv( self.data_folder_name + "variant.csv" )
        self.variant_df['GENE'] = self.variant_df.IDENTIFIER
        self.cell_line_df = pd.read_csv( self.data_folder_name + "cell_line_dict.csv" )
        if annotate:
            self.annotate_all()
        else:
            self.cna_df = pd.read_csv( self.data_folder_name + "cna.csv" )
            self.methyl_df = pd.read_csv( self.data_folder_name + "methylation.csv" )

    def annotate_all(self):
        self.annotate = True
        self._download_annotation()
        self.methyl_df = self._annotate_methylation()
        self.cna_df = self._annotate_cna()

    def merge_all( self ):
        """
        Concatanates genomic annotations and adds cell line information.
        If CNA and methylation are annotated, an extra 'GENE' column will be included in self.complete_df
        """
        self.variant_df['TYPE'] = 'GENE_VARIANT'
        self.cna_df['TYPE'] = self.cna_df.ALTERATION_TYPE.str.upper()
        self.methyl_df['TYPE'] = 'METHYLATION'
        self.complete_df = pd.concat( [self.variant_df, self.cna_df, self.methyl_df] , axis = 0, join = 'inner' )
        self.complete_df = self.complete_df.merge( self.cell_line_df, how = 'left', on = ['CELL_LINE', 'TISSUE_FACTOR'])

    def _download_annotation(self):
        """
        Downloads the annotation files from the cancerrxgene.org/gdsc1000 website.
        """
        methylation_annotation_url = self.url_base + "TableS2H.xlsx"
        cna_annotation_url = self.url_base + "TableS2D.xlsx"
        urlrequest.urlretrieve( methylation_annotation_url, self.data_folder_name + "methylation_annotation.xlsx")
        urlrequest.urlretrieve( cna_annotation_url, self.data_folder_name + "cna_annotation.xlsx")

    def _annotate_methylation(self):
        methyl_df = pd.read_csv(self.data_folder_name + "methylation.csv")
        xls = pd.ExcelFile(self.data_folder_name + "methylation_annotation.xlsx")
        methyl_annotation = xls.parse(header = 16)
        methyl_annotation = methyl_annotation[ ['Cancer Types', 'Genomic Coordinates', 'GN' ] ]
        methyl_annotation.index = range( len( methyl_annotation.index ) )
        methyl_annotation.columns = [ 'TISSUE_FACTOR', 'IDENTIFIER', 'GENE' ]
        methyl_df = methyl_df.merge( methyl_annotation, how = 'left', on = ['IDENTIFIER', 'TISSUE_FACTOR'] )
        methyl_df = self._string_split( methyl_df, 'GENE', '; ')
        return methyl_df

    def _annotate_cna(self):
        cna_df = pd.read_csv(self.data_folder_name + "cna.csv")
        xls = pd.ExcelFile(self.data_folder_name + "cna_annotation.xlsx")
        cna_annotation = xls.parse(header = 19)
        cna_annotation = cna_annotation[ [ 'Identifier', 'Cancer Type', 'Recurrent\nAmplification/\nDeletion', 'nGenes', 'Contained genes' ] ]
        cna_annotation.columns = ['IDENTIFIER', 'TISSUE_FACTOR', 'CNA_TYPE', 'NGENE', 'GENE']
        cna_df = cna_df.merge( cna_annotation, how = 'left', on = ['IDENTIFIER', 'TISSUE_FACTOR'] )
        cna_df = self._string_split( cna_df, 'GENE', ',' )
        return cna_df


    def _string_split(self, to_be_split, col_name, separator):
        """
        This function splits multiple genes contained within the same genomic region across multiple rows, duplicating all other values.
        """
        s = to_be_split[col_name].str.split( separator ).apply(pd.Series, 1).stack()
        s.index = s.index.droplevel(-1) # Flattens levels to retrieve indices from original frame
        s.name = col_name # Join requires Series name
        del to_be_split[col_name]
        split_df = to_be_split.join(s)
        return split_df
