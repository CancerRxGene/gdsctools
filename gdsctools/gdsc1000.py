"""

This module of gdsctools provides an easy interface for anybody to 
download and play with the GDSC1000 data, as published by Iorio et al (2016).

Author: Elisabeth D. Chen
Date: 2016-06-17

Please see gdsctools.gdsc1000 module for other data sets

"""
import os

try:     #python 3
    from urllib.request import urlretrieve
except:
    from urllib import urlretrieve 

import csv

import pandas as pd

from gdsctools.default_sets import DefaultDictionaries

import colorlog as logger

class GDSC1000(object ):
    """Help data retrieval from GDSC website (version v17)

    ::

        from gdsctools import GDSC1000
        data = GDSC1000()
        data.download_data()

    .. autosummary::

        download_data
        load_data
        make_matrix


    With this class, you may (1) download the data, (2) annotate or filter the
    data and (3) create a genomic matrix.

    To load the data, either download it from the website using
    :meth:`download_data` (downloads data from cancerrxgene.org loading methylation, cna, 
    variant and cell lines datasets). It extract relevant columns, re-names column names 
    and saves as csv files in the :attr:`data_folder_name` directory.

    If you have already downloaded the data, you may just load it using
    :meth:`load_data`.

    Then, you can annotate or filter the genomic data using :meth:`annotate_all` 
    or one of the filter methods. ::

        filter_by_cell_line([ "AsPC-1", "U-2-OS", "MDA-MB-231"... ] )
        filter_by_cosmic_id([ 948121, 2839818, ...] )
        filter_by_gene(["ATM", "TP53", ...])
        filter_by_recurrence(min_recurrence = 3 )
        filter_by_tissue_type([ "BRCA", "COAD/READ", ... ] )
        filter_by_alteration_type(["METHYLATION", "AMPLIFICATION",
                                        "DELETION", "GENE_VARIANT" ] )

    Finally, you can create the genomic matrix using :meth:`make_matrix`. 

    """
    def __init__(self, annotate=False, data_folder_name="./data/gdsc_1000_data/" ):
        """.. rubric:: constructor

        :param bool annotate:
        :param str data_folder_name:

        """
        self._save_csv = True
        self.annotate = annotate
        self.recurrence = 3
        self._data_folder_name = data_folder_name
        self._url_base = "http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/suppData/"

        self._mapper = {"methylation": ["TableS2J.xlsx", "methylation.xlsx"],
                        "cna": ["TableS2G.xlsx", "cna.xlsx"],
                        "variant": ["TableS2C.xlsx", "variant.xlsx"],
                        "cell_line": ["TableS1E.xlsx", "cell_line_dict.xlsx"],
                        "methylation_annotation": ["TableS2H.xlsx", "methylation_annotation.xlsx"],
                        "cna_annotation": ["TableS2D.xlsx", "cna_annotation.xlsx"],
                        "ic50": ["TableS4A.xlsx", "ic50.xlsx"],
                        }

        self.defaults = DefaultDictionaries()

    def _get_path(self, key):
        return self._data_folder_name + self._mapper[key][1]

    def download_data(self):
        """Download and load data in memory"""
        self._download_raw()
        self._convert_raw_gdsc1000_data()

    def _convert_raw_gdsc1000_data(self):
        """Converts raw Excel document into CSV"""
        self._read_cell_line_data()
        self._read_cna_data()
        self._read_methylation_data()
        self._read_variant_data()
        self._read_ic50_data()
        self._read_annotation()
        self._combine_genomic_annotation()

    def _to_csv(self, df, filename, quoting=None):
        if self._save_csv:
            logger.info("Saving %s" % filename)
            df.to_csv(self._data_folder_name + filename, index=False)

    def _read_csv(self, filename):
        return pd.read_csv(self._data_folder_name + filename)

    def load_data(self):
        """If CSV files are already downloaded, just load them"""
        self.variant_df = self._read_csv("variant.csv")
        self.methylation_df = self._read_csv("methylation.csv")
        self.cna_df = self._read_csv("cna.csv")
        self.cell_line_dict = self._read_csv("cell_line_dict.csv")
        self.ic50_df = self._read_csv("ic50_full_matrix.csv")

        self.drug_decoder_df = self._read_csv("drug_decoder.csv")
        self.genomic_df = self._read_csv("genomic.csv")
        self.annotation_df = self._read_csv("annotation.csv")

        self.backup_genomic_df = self.genomic_df
        self.backup_ic50_df = self.ic50_df
        self.backup_drug_decoder_df = self.drug_decoder_df

    def _urlretrieve(self, filename, target):
        urlretrieve(
            self._url_base + filename,
            self._data_folder_name + target)

    def _download_raw(self):

        if os.path.exists(self._data_folder_name):
            pass
        else:
            os.makedirs(self._data_folder_name)

        # download all data sets
        for key, value in self._mapper.items():
            logger.info('Downloading %s ' % key)
            self._urlretrieve(*value)

    def _read_cna_data(self):
        logger.info("Processing CNA data.")
        # This block reads in and formats CNA data
        xls = pd.ExcelFile(self._get_path('cna'))
        df = xls.parse(header=2)

        # in version from feb 2017, there are 4 extra blank columns to be
        # drop
        dummies = ["dum1", "dum2", "dum3", "dum4"]
        df.columns = ["BLANK", "CELL_LINE", "TISSUE_ID1", "TISSUE_ID2",
                      "TISSUE_FACTOR", "ALTERATION_TYPE", "IDENTIFIER" ] + dummies

        df = df.drop(['BLANK', 'TISSUE_ID1', 'TISSUE_ID2'] + dummies, axis = 1)
        df.IDENTIFIER = df.IDENTIFIER.str.split(' ').str[0]
        df.ALTERATION_TYPE = df.ALTERATION_TYPE.str.upper()

        self._to_csv(df, 'cna.csv')
        self.cna_df = df

    def _read_methylation_data(self):
        logger.info("Processing methylation data.")
        # This block deals with the methylation data
        xls = pd.ExcelFile(self._get_path('methylation'))
        df = xls.parse(header=2)
        df.columns = ["BLANK", "CELL_LINE", "TISSUE_ID1", "TISSUE_ID2",
                      "TISSUE_FACTOR", "IDENTIFIER"]
        df = df.drop(["BLANK", "TISSUE_ID1", "TISSUE_ID2"], axis=1)
        df['ALTERATION_TYPE'] = "METHYLATION"
        self._to_csv(df, 'methylation.csv')
        self.methylation_df = df

    def _read_variant_data(self):
        logger.info("Processing variant data.")
        # This block deals with variant data
        xls = pd.ExcelFile(self._get_path('variant'))
        df = xls.parse(skiprows=20 )
        df = df[ ['SAMPLE', 'COSMIC_ID', 'Cancer Type', 'Gene', 'Recurrence Filter' ] ]
        df = df.loc[df['Recurrence Filter'] == 'Yes' ]
        df.columns = ['CELL_LINE', 'COSMIC_ID', 'TISSUE_FACTOR', 'IDENTIFIER',
                      'RECURRENCE' ]
        df = df.drop('RECURRENCE', axis=1)
        df['ALTERATION_TYPE'] = "GENETIC_VARIATION"
        self._to_csv(df, 'variant.csv')
        self.variant_df = df

    def _read_cell_line_data(self):
        logger.info("Processing cell line dictionary.")
        # This block imports the cell line dictionary
        xls = pd.ExcelFile(self._get_path("cell_line"))
        df = xls.parse(header = 2 )
        df.columns = ["CELL_LINE", "COSMIC_ID", "WES", "CNA", "GEX",
                      "METH", "DRUGDATA", "T_ID1", "T_ID2", "TISSUE_FACTOR",
                      "MSI_FACTOR", "MEDIA_FACTOR", "GROWTH_FACTOR" ]
        df.index = range(len(df.index ) )
        df = df.drop(0)
        df = df[['CELL_LINE', 'COSMIC_ID', 'TISSUE_FACTOR', 'MSI_FACTOR',
                 'MEDIA_FACTOR', 'GROWTH_FACTOR' ] ]
        self._to_csv(df, 'cell_line_dict.csv')
        self.cell_line_dict = df

    def _read_ic50_data(self):
        logger.info("Processing IC50 data." )
        xls = pd.ExcelFile(self._get_path('ic50') )
        df = xls.parse(header = 4 )
        df.columns.values[ 0 ] = 'COSMIC_ID'
        df = df.drop(0)
        df = df.drop('Sample Names', axis = 1 )
        self.ic50_df = df

        drug_dict = xls.parse(header=4)[0:1]
        drug_dict.drop(drug_dict.columns[[0, 1]], axis=1, inplace=True)
        drug_dict = drug_dict.T.reset_index()
        drug_dict.columns = ['DRUG_ID', 'DRUG_NAME']
        drug_dict['DRUG_TARGET'] = 'Unknown'
        self.drug_decoder_df = drug_dict

        self._to_csv(df, "ic50_full_matrix.csv")
        self._to_csv(drug_dict, "drug_decoder.csv")

        self.backup_ic50_df = self.ic50_df
        self.backup_drug_decoder_df = self.drug_decoder_df

    def _combine_genomic_annotation(self):
        logger.info("Combining data frames." )
        self.genomic_df = pd.concat([self.variant_df, self.cna_df,
                                     self.methylation_df] ,
                                     axis=0, join='inner')
        self.genomic_df = self.genomic_df.merge(self.cell_line_dict,
                                    how='left', on=['CELL_LINE', 'TISSUE_FACTOR'])
        self._to_csv(self.genomic_df, 'genomic.csv')
        self.backup_genomic_df = self.genomic_df

    def reset_genomic_data(self ):
        self.genomic_df = self.backup_genomic_df
        self.annotate = False
        self.ic50_df = self.backup_ic50_df
        self.drug_decoder_df = self.backup_drug_decoder_df

    def annotate_all(self ):
        """Annotates dataframes after read_annotation"""
        logger.info("Annotating data" )
        self.genomic_df = self.genomic_df.merge(self.annotation_df,
            how='left', on=['IDENTIFIER'])
        self.genomic_df = self._string_split(self.genomic_df, 'GENE', ',' )
        self.annotate = True

    def _string_split(self, to_be_split, col_name, separator):
        """This function splits multiple genes contained within the
        same genomic region across multiple rows, duplicating all other values.
        """
        s = to_be_split[col_name].str.split(separator ).apply(pd.Series, 1).stack()
        # Making index
        s.index = s.index.droplevel(-1) # Flattens levels to retrieve indices from original frame
        # Defining name
        s.name = col_name # Join requires Series name
        #  Delete column in original dataframe
        del to_be_split[col_name]
        # Join data frames.
        split_df = to_be_split.join(s)
        # Return split data frame
        return split_df

    def _read_annotation(self):
        # Read and merge annotations
        xls = pd.ExcelFile(self._get_path("methylation_annotation"))

        methylation_annotation = xls.parse(header = 16)
        methylation_annotation = methylation_annotation[ ['Genomic Coordinates', 'GN' ] ]
        methylation_annotation.columns = [ 'IDENTIFIER', 'GENE' ]
        methylation_annotation.GENE = methylation_annotation.GENE.replace(
            to_replace = "; ", value=",", regex=True)
        methylation_annotation = methylation_annotation[
            methylation_annotation.IDENTIFIER.duplicated() == False ]

        xls = pd.ExcelFile(self._get_path("cna_annotation"))
        cna_annotation = xls.parse(header = 19)
        cna_annotation = cna_annotation[ [ 'Identifier', 'Contained genes' ] ]
        cna_annotation.columns = ['IDENTIFIER', 'GENE']

        variant_annotation = pd.DataFrame({
            'IDENTIFIER': self.variant_df.IDENTIFIER.unique(),
            'GENE': self.variant_df.IDENTIFIER.unique()})
        # CHECK FOR UNIQUENESS

        self.annotation_df = pd.concat([ methylation_annotation, 
            cna_annotation, variant_annotation ], axis = 0, join = "inner" )

        self._to_csv(self.annotation_df, "annotation.csv")

    def _check_filter(self, values, valid_values):
        if values is None:
            raise TypeError("Please enter list of types to keep. Valid types: " + 
                ", ".join(valid_values))
        for this in values:
            if this not in valid_values:
                raise ValueError("%s not valid. Use one of " % this + 
                    " ".join(valid_values))

    def filter_by_alteration_type(self, type_list=None):
        valid_types = self.genomic_df.ALTERATION_TYPE.unique()
        self._check_filter(type_list, valid_types)
        type_list = [ x.upper() for x in type_list ]
        self.genomic_df = self.genomic_df.query("ALTERATION_TYPE in @type_list" )

    def filter_by_gene(self, gene_list=None):
        if gene_list is None:
            print("Please enter a valid gene list or set to 'Core Genes' to use default GDSC1000 core genes" )
            return
        elif isinstance(gene_list, str ):
            gene_list = self.defaults.gene_dict[gene_list]
        gene_list = [ x.upper() for x in gene_list ]
        self.genomic_df = self.genomic_df.query("IDENTIFIER in @gene_list" )

    def filter_by_tissue_type(self, tissue_list=None):
        if tissue_list == None:
            print("Please enter list of tissues. To use all tissues, add argument tissue_list = 'Complete'" )
            print("Possible dictionary keys include: ", self.defaults.tissue_dict.keys() )
            return
        elif isinstance(tissue_list, str ):
            tissue_list = self.defaults.tissue_dict[ tissue_list ]
        tissue_list = [ x.upper() for x in tissue_list ]
        self.genomic_df = self.genomic_df.query("TISSUE_FACTOR in @tissue_list" )

    def filter_by_cosmic_id(self, cosmic_list=None):
        if cosmic_list == None:
            print("Please enter list of COSMIC IDs. To use an example set of COSMIC IDs, add argument cosmic_list = 'Dummy'" )
            print("Possible dictionary keys include: ", self.defaults.cosmic_dict.keys() )
            return
        elif isinstance(cosmic_list, str ):
            cosmic_list = self.defaults.cosmic_dict[ cosmic_list ]
        self.genomic_df = self.genomic_df.query("COSMIC_ID in @cosmic_list" )

    def filter_by_cell_line(self, cell_line_list = None ):
        if cell_line_list == None:
            print("Please enter list of cell line names. To use an example set of names, add argument cell_line_list = 'Dummy'" )
            print("Possible dictionary keys include: ", self.defaults.cell_line_dict.keys() )
            return
        elif isinstance(cell_line_list, str ):
            cell_line_list = self.defaults.cell_line_dict[cell_line_list]

        self.genomic_df = self.genomic_df.query("CELL_LINE in @cell_line_list")

    def filter_by_recurrence(self, min_recurrence=None):
        if min_recurrence == None:
            min_recurrence = self.recurrence

        if self.annotate == True:
            feature = "GENE"
        else:
            feature = "IDENTIFIER"

        flatten_cell_lines = self.genomic_df.groupby([ feature, "COSMIC_ID", "TISSUE_FACTOR" ], 
            as_index = False )[[ "ALTERATION_TYPE" ]].count()
        feature_count = flatten_cell_lines.groupby([feature], as_index = False )[[ "COSMIC_ID" ]].count()
        self.feature_count = feature_count[ feature_count["COSMIC_ID"] >= min_recurrence ]

        self.genomic_df = self.genomic_df[self.genomic_df[feature].isin(self.feature_count[feature])]

    # Make into matrix
    def make_matrix(self, min_recurrence=None):
        """Return a dataframe compatible with ANOVA analysis"""
        msi_dictionary = {"MSI-H": 1, "MSS/MSI-L": 0}

        if min_recurrence == None:
            min_recurrence = self.recurrence
        self.filter_by_recurrence(min_recurrence=min_recurrence)

        if self.annotate == True:
            feature = "GENE"
        else:
            feature = "IDENTIFIER"

        genomic_matrix = pd.crosstab(self.genomic_df[ feature ], 
            columns=[ self.genomic_df["COSMIC_ID"], self.genomic_df['TISSUE_FACTOR'], 
                      self.genomic_df["CELL_LINE"], self.genomic_df[ "MSI_FACTOR" ],
                      self.genomic_df[ "MEDIA_FACTOR"] ] )
        genomic_matrix[ genomic_matrix > 1 ] = 1
        self.genomic_matrix = genomic_matrix.T.reset_index()
        self.genomic_matrix.MSI_FACTOR = self.genomic_matrix.MSI_FACTOR.map(msi_dictionary)

        self._to_csv(self.genomic_matrix, "genomic_features_make_matrix.csv", quoting=csv.QUOTE_NONNUMERIC)


        ic50_filtered_matrix = self.ic50_df.query("COSMIC_ID in @self.genomic_matrix.COSMIC_ID.tolist()" )
        self._to_csv(ic50_filtered_matrix, "ic50_make_matrix.csv")
