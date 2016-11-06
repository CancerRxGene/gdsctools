"""

This module of gdsctools provides an easy interface for anybody to download and play with the GDSC1000 data, as published by Iorio et al (2016).

Author: Elisabeth D. Chen
Date: 2016-06-17

"""
import pandas as pd
import os
import urllib.request
import csv


class GDSC1000(object ):
    """Main GDSC class for data retrival

    1. Load data:
        GDSC1000.download_data()
            Downloads data from cancerrxgene.org for methylation, cna, variant and cell lines.
            Extracts relevant columns, re-names column names and saves as csv files.
        GDSC1000.convert_raw_gdsc1000_data()
            Loads excel frames locally for methylation, cna, variant and cell lines.
            Extracts relevant columns, re-names column names and saves as csv files.
        GDSC1000.load_data()
            Loads methylation, cna and variant data frames from csv files.
    2. Annotate / Filter data:
        GDSC1000.annotate_all()
            Annotates methylation and copy number alteration on a gene level
            Also sets self.annotate = True (only really used when calling GDSC1000.make_matrix() )
        GDSC1000.filter(...):
            GDSC1000.filter_by_cell_line([ "AsPC-1", "U-2-OS", "MDA-MB-231"... ] )
            GDSC1000.filter_by_cosmic_id([ 948121, 2839818, ...] )
            GDSC1000.filter_by_gene(["ATM", "TP53", ...])
            GDSC1000.filter_by_recurrence(min_recurrence = 3 )
            GDSC1000.filter_by_tissue_type([ "BRCA", "COAD/READ", ... ] )
            GDSC1000.filter_by_alteration_type(["METHYLATION", "AMPLIFICATION", 
                                                "DELETION", "GENE_VARIANT" ] )
    3. Make genomic matrix:
        GDSC1000.make_matrix()
            Makes matrix and stores it under GDSC1000.genomic_matrix_path
            Filters GDSC1000.ic50_df to only contain cell lines also 
            featured in GDSC1000.genomic_matrix

    """
    def __init__(self, annotate=False, data_folder_name="./data/gdsc_1000_data/" ):
        """.. rubric:: constructor

        :param bool annotate:
        :param str data_folder_name:


        """
        self.annotate = annotate
        self.recurrence = 3
        self.data_folder_name = data_folder_name
        self.url_base = "http://www.cancerrxgene.org/gdsc1000/Data/suppData/"
        print("Welcome to GDSCTools 1000.\n",
               "To download the GDSC 1000 data, import" )

        self._methylation_path = self.data_folder_name + "methylation.xlsx"
        self._cna_path = self.data_folder_name + "cna.xlsx"
        self._variant_path = self.data_folder_name + "variant.xlsx"
        self._cell_dict_path = self.data_folder_name + "cell_line_dict.xlsx"
        self._ic50_path = self.data_folder_name + "ic50.xlsx"

        self.genomic_matrix_path = self.data_folder_name + "genomic_features_matrix.csv"
        self.ic50_path = self.data_folder_name + "ic50_full_matrix.csv"
        self.ic50_matrix_path = self.data_folder_name + "ic50_matrix.csv"
        self.drug_decoder_original_path = self.data_folder_name + "drug_decoder.csv"

    def download_data(self, save_csv=False ):
        self._download_raw()
        self.convert_raw_gdsc1000_data(save_csv=save_csv )

    def convert_raw_gdsc1000_data(self, save_csv=True ):
        self._read_cell_line_data(save_csv = save_csv)
        self._read_cna_data(save_csv = save_csv)
        self._read_methylation_data(save_csv = save_csv)
        self._read_variant_data(save_csv = save_csv)
        self._read_ic50_data(save_csv = save_csv )
        self._combine_genomic_annotation(save_csv = save_csv )
        self._read_annotation(save_csv = save_csv )

    def _read_csv(self, filename):
        return pd.read_csv(self.data_folder_name + filename)

    def load_data(self):
        self.variant_df =       self._read_csv("variant.csv")
        self.methylation_df =   self._read_csv("methylation.csv")
        self.cna_df =           self._read_csv("cna.csv")
        self.cell_line_dict =   self._read_csv("cell_line_dict.csv")
        self.ic50_df =          pd.read_csv(self.ic50_path)
        self.drug_decoder_df =  pd.read_csv(self.drug_decoder_original_path)
        self.genomic_df =       self._read_csv("genomic.csv")
        self.annotation_df =    self._read_csv("annotation.csv")

        self.backup_genomic_df = self.genomic_df
        self.backup_ic50_df = self.ic50_df
        self.backup_drug_decoder_df = self.drug_decoder_df

    def _download_raw(self):
        print("Downloading data...")

        methylation_url = self.url_base +  "TableS2J.xlsx"
        cna_url = self.url_base + "TableS2G.xlsx"
        variant_url = self.url_base + "TableS2C.xlsx"
        cell_line_url = self.url_base + "TableS1E.xlsx"
        ic50_url = self.url_base + "TableS4A.xlsx"
        methylation_annotation_url = self.url_base + "TableS2H.xlsx"
        cna_annotation_url = self.url_base + "TableS2D.xlsx"

        if os.path.exists(self.data_folder_name):
            pass
        else:
            os.makedirs(self.data_folder_name)

        urllib.request.urlretrieve(methylation_url, self._methylation_path)
        urllib.request.urlretrieve(cna_url, self._cna_path)
        urllib.request.urlretrieve(variant_url, self._variant_path)
        urllib.request.urlretrieve(cell_line_url, self._cell_dict_path)
        urllib.request.urlretrieve(ic50_url, self.data_folder_name + self._ic50_path )
        urlrequest.urlretrieve(methylation_annotation_url,
                               self.data_folder_name + "methylation_annotation.xlsx")
        urlrequest.urlretrieve(cna_annotation_url,
                               self.data_folder_name + "cna_annotation.xlsx")

    def _read_cna_data(self, save_csv):
        print("Processing CNA data.")
        # This block reads in and formats CNA data
        xls = pd.ExcelFile(self._cna_path )
        df = xls.parse(header=2)
        df.columns = ["BLANK", "CELL_LINE", "TISSUE_ID1", "TISSUE_ID2",
                      "TISSUE_FACTOR", "ALTERATION_TYPE", "IDENTIFIER" ]
        df = df.drop(['BLANK', 'TISSUE_ID1', 'TISSUE_ID2'], axis = 1)
        df.IDENTIFIER = df.IDENTIFIER.str.split(' ').str[0]
        df.ALTERATION_TYPE = df.ALTERATION_TYPE.str.upper()
        if save_csv:
            df.to_csv(self.data_folder_name + 'cna.csv', index=False)
        self.cna_df = df

    def _read_methylation_data(self, save_csv):
        print("Processing methylation data.")
        # This block deals with the methylation data
        xls = pd.ExcelFile(self._methylation_path )
        df = xls.parse(header=2)
        df.columns = ["BLANK", "CELL_LINE", "TISSUE_ID1", "TISSUE_ID2",
                      "TISSUE_FACTOR", "IDENTIFIER"]
        df = df.drop(["BLANK", "TISSUE_ID1", "TISSUE_ID2"], axis=1)
        df['ALTERATION_TYPE'] = "METHYLATION"
        if save_csv:
            df.to_csv(self.data_folder_name + 'methylation.csv', index=False)
        self.methylation_df = df

    def _read_variant_data(self, save_csv):
        print("Processing variant data.")
        # This block deals with variant data
        print("Reading variant data.")
        xls = pd.ExcelFile(self._variant_path )
        df = xls.parse(skiprows=20 )
        df = df[ ['SAMPLE', 'COSMIC_ID', 'Cancer Type', 'Gene', 'Recurrence Filter' ] ]
        df = df.loc[df['Recurrence Filter'] == 'Yes' ]
        df.columns = ['CELL_LINE', 'COSMIC_ID', 'TISSUE_FACTOR', 'IDENTIFIER',
                      'RECURRENCE' ]
        df = df.drop('RECURRENCE', axis=1)
        df['ALTERATION_TYPE'] = "GENETIC_VARIATION"
        if save_csv:
            print("Writing variant csv." )
            df.to_csv(self.data_folder_name + 'variant.csv', index=False )
        self.variant_df = df

    def _read_cell_line_data(self, save_csv):
        print("Processing cell line dictionary.")
        # This block imports the cell line dictionary
        xls = pd.ExcelFile(self._cell_dict_path )
        df = xls.parse(header = 2 )
        df.columns = ["CELL_LINE", "COSMIC_ID", "WES", "CNA", "GEX",
                      "METH", "DRUGDATA", "T_ID1", "T_ID2", "TISSUE_FACTOR",
                      "MSI_FACTOR", "MEDIA_FACTOR", "GROWTH_FACTOR" ]
        df.index = range(len(df.index ) )
        df = df.drop(0)
        df = df[['CELL_LINE', 'COSMIC_ID', 'TISSUE_FACTOR', 'MSI_FACTOR',
                 'MEDIA_FACTOR', 'GROWTH_FACTOR' ] ]
        if save_csv:
            df.to_csv(self.data_folder_name + 'cell_line_dict.csv', index = False )
        self.cell_line_dict = df

    def _read_ic50_data(self, save_csv ):
        print("Processing IC50 data." )
        xls = pd.ExcelFile(self._ic50_path )
        df = xls.parse(header = 4 )
        df.columns.values[ 0 ] = 'COSMIC_ID'
        df = df.drop(0)
        df = df.drop('Sample Names', axis = 1 )
        self.ic50_df = df

        drug_dict = xls.parse(header = 4 )[ 0:1 ]
        drug_dict.drop(drug_dict.columns[[ 0, 1]], axis = 1, inplace = True )
        drug_dict = drug_dict.T.reset_index()
        drug_dict.columns = [ 'DRUG_ID', 'DRUG_NAME' ]
        drug_dict[ 'DRUG_TARGET' ] = 'Unknown'
        self.drug_decoder_df = drug_dict
        if save_csv:
            df.to_csv(self.ic50_path, index = False )
            drug_dict.to_csv(self.drug_decoder_original_path, index = False )

        self.backup_ic50_df = self.ic50_df
        self.backup_drug_decoder_df = self.drug_decoder_df

    def _combine_genomic_annotation(self, save_csv ):
        print("Combining data frames." )
        self.genomic_df = pd.concat([self.variant_df, self.cna_df,
                                     self.methylation_df] ,
                                     axis=0, join='inner')
        self.genomic_df = self.genomic_df.merge(self.cell_line_dict,
                                    how='left', on=['CELL_LINE', 'TISSUE_FACTOR'])
        print("Saving genomic data frame." )
        if save_csv:
            self.genomic_df.to_csv(self.data_folder_name + 'genomic.csv', index = False )
        self.backup_genomic_df = self.genomic_df

    def reset_genomic_data(self ):
        self.genomic_df = self.backup_genomic_df
        self.annotate = False
        self.ic50_df = self.backup_ic50_df
        self.drug_decoder_df = self.backup_drug_decoder_df


    def annotate_all(self ):
        # Annotate after read_annotation
        print("Merging" )
        self.genomic_df = self.genomic_df.merge(self.annotation_df,
            how='left', on=['IDENTIFIER'])
        print("Splitting strings" )
        self.genomic_df = self._string_split(self.genomic_df, 'GENE', ',' )
        self.annotate = True


    def _string_split(self, to_be_split, col_name, separator):
        """This function splits multiple genes contained within the 
        same genomic region across multiple rows, duplicating all other values.
        """
        print("Making s" )
        s = to_be_split[col_name].str.split(separator ).apply(pd.Series, 1).stack()
        print("Making index" )
        s.index = s.index.droplevel(-1) # Flattens levels to retrieve indices from original frame
        print("Defining name" )
        s.name = col_name # Join requires Series name
        print("Delete column in original" )
        del to_be_split[col_name]
        print("Join data frames." )
        split_df = to_be_split.join(s)
        print("Return split data frame" )
        return split_df

    def _read_annotation(self, save_csv = True ):
        # Read and merge annotations
        xls = pd.ExcelFile(self.data_folder_name + "methylation_annotation.xlsx")
        methylation_annotation = xls.parse(header = 16)
        methylation_annotation = methylation_annotation[ ['Genomic Coordinates', 'GN' ] ]
        methylation_annotation.columns = [ 'IDENTIFIER', 'GENE' ]
        methylation_annotation.GENE = methylation_annotation.GENE.replace(to_replace = "; ", value = ",", regex = True )
        methylation_annotation = methylation_annotation[ methylation_annotation.IDENTIFIER.duplicated() == False ]

        xls = pd.ExcelFile(self.data_folder_name + "cna_annotation.xlsx")
        cna_annotation = xls.parse(header = 19)
        cna_annotation = cna_annotation[ [ 'Identifier', 'Contained genes' ] ]
        cna_annotation.columns = ['IDENTIFIER', 'GENE']

        variant_annotation = pd.DataFrame({ 'IDENTIFIER': self.variant_df.IDENTIFIER.unique(), 'GENE': self.variant_df.IDENTIFIER.unique() } )
        # CHECK FOR UNIQUENESS

        self.annotation_df = pd.concat([ methylation_annotation, cna_annotation, variant_annotation ], axis = 0, join = "inner" )

        if save_csv:
            self.annotation_df.to_csv(self.data_folder_name + "annotation.csv", index = False )

    def filter_by_alteration_type(self, type_list = None ):
        if type_list == None:
            print("Please enter list of types to keep. Acceptable types include: 'GENE_VARIANT', 'AMPLIFICATION', 'DELETION', 'METHYLATION'." )
        else:
            type_list = [ x.upper() for x in type_list ]
            self.genomic_df = self.genomic_df.query("TYPE in @type_list" )

    def filter_by_gene(self, gene_list = None ):
        if gene_list == None:
            print("Please enter gene list. To use the GDSC1000 core genes, add argument gene_list = 'Core Genes'." )
            if self.dicts in globals():
                print("Possible dictionary keys include: ", self.dicts.gene_dict.keys() )
        elif isinstance(gene_list, str ):
            gene_list = self.dicts.gene_dict[ gene_list ]
        else:
            pass
        gene_list = [ x.upper() for x in gene_list ]
        self.genomic_df = self.genomic_df.query("GENE in @gene_list" )

    def filter_by_tissue_type(self, tissue_list = None ):
        if tissue_list == None:
            print("Please enter list of tissues. To use all tissues, add argument tissue_list = 'Complete'" )
            if self.dicts in globals():
                print("Possible dictionary keys include: ", self.dicts.tissue_dict.keys() )
        elif isinstance(tissue_list, str ):
            tissue_list = self.dicts.tissue_dict[ tissue_list ]
        else:
            pass
        tissue_list = [ x.upper() for x in tissue_list ]
        self.genomic_df = self.genomic_df.query("TISSUE_FACTOR in @tissue_list" )

    def filter_by_cosmic_id(self, cosmic_list = None ):
        if cosmic_list == None:
            print("Please enter list of COSMIC IDs. To use an example set of COSMIC IDs, add argument cosmic_list = 'Dummy'" )
            if self.dicts in globals():
                print("Possible dictionary keys include: ", self.dicts.cosmic_dict.keys() )
        elif isinstance(cosmic_list, str ):
            cosmic_list = self.dicts.cosmic_dict[ cosmic_list ]
        else:
            pass
        self.genomic_df = self.genomic_df.query("COSMIC_ID in @cosmic_list" )

    def filter_by_cell_line(self, cell_line_list = None ):
        if cell_line_list == None:
            print("Please enter list of cell line names. To use an example set of names, add argument cell_line_list = 'Dummy'" )
            if self.dicts in globals():
                print("Possible dictionary keys include: ", self.dicts.cell_line_dict.keys() )
        elif isinstance(cell_line_list, str ):
            cell_line_list = self.dicts.cell_line_dict[ cell_line_list ]
        else:
            pass
        self.genomic_df = self.genomic_df.query("CELL_LINE in @cell_line_list" )

    def filter_by_recurrence(self, min_recurrence = None ):
        if min_recurrence == None:
            min_recurrence = self.recurrence

        if self.annotate == True:
            feature = "GENE"
        else:
            feature = "IDENTIFIER"

        flatten_cell_lines = self.genomic_df.groupby([ feature, "COSMIC_ID", "TISSUE_FACTOR" ], as_index = False )[[ "ALTERATION_TYPE" ]].count()
        feature_count = flatten_cell_lines.groupby([feature], as_index = False )[[ "COSMIC_ID" ]].count()
        self.feature_count = feature_count[ feature_count["COSMIC_ID"] >= min_recurrence ]

        self.genomic_df = self.genomic_df[ self.genomic_df[ feature ].isin(self.feature_count[ feature ] ) ]


    # Make into matrix
    def make_matrix(self, min_recurrence = None ):
        """
        Return a dataframe compatible with ANOVA analysis
        """

        msi_dictionary = { "MSI-H": 1, "MSS/MSI-L": 0 }

        if min_recurrence == None:
            min_recurrence = self.recurrence
        self.filter_by_recurrence(min_recurrence = 3 )

        if self.annotate == True:
            feature = "GENE"
        else:
            feature = "IDENTIFIER"

        genomic_matrix = pd.crosstab(self.genomic_df[ feature ], columns=[ self.genomic_df["COSMIC_ID"], self.genomic_df['TISSUE_FACTOR'], self.genomic_df["CELL_LINE"], self.genomic_df[ "MSI_FACTOR" ], self.genomic_df[ "MEDIA_FACTOR"] ] )
        genomic_matrix[ genomic_matrix > 1 ] = 1
        self.genomic_matrix = genomic_matrix.T.reset_index()
        self.genomic_matrix.MSI_FACTOR = self.genomic_matrix.MSI_FACTOR.map(msi_dictionary )

        self.genomic_matrix.to_csv(self.genomic_matrix_path, index = False, quoting = csv.QUOTE_NONNUMERIC )

        ic50_filtered_matrix = self.ic50_df.query("COSMIC_ID in @self.genomic_matrix.COSMIC_ID.tolist()" )
        ic50_filtered_matrix.to_csv(self.ic50_matrix_path, index = False )
