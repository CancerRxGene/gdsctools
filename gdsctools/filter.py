"""

To filter complete_df that is outputted from annotate_genes.py
E.D.Chen 2016-06-05


See gdsctools.gdsc1000 instead
"""

import pandas as pd

class FilterFeatures( object ):
    """

    fl = FilterFeatures( annotated_class ) # This function takes in an object of class annotate_genes that contains at least values for self.annotate and self.complete_df

    Filters that can then be applied include:

    filter_by_type( [ 'GENE_VARIANT', 'AMPLIFICATION', 'DELETION', 'METHYLATION' ] )
    filter_by_tissue_list( ['BRCA', 'SCLSC', ...] )
    filter_by_gene_list( [ 'TP53', 'APC', ... ] ) # Takes list of gene names in HUGO format
    filter_by_cell_line_list( [ 'U-2-OS', 'AsPC-1', '...' ] ) # Takes list of cell line names. Case sensitive!!
    filter_by_cosmic_list( [ ... ] ) # Takes list of COSMIC IDs (must be integer or floats)

    filter_by_recurrence( min_recurrence = 3 ) # Only keeps genomic features that are altered in at least x cell lines. By default, x = min_recurrence = 3


    At the end, call the function:

    fl.make_matrix()
    # Automatically calls the function filter_by_recurrence with min_recurrence = 3.
    # Stores final matrix that can be fed into gdsctools ANOVA part under fl.final_matrix

    """
    def __init__( self, annotated_class ):
        self.annotate = annotated_class.annotate
        self.starting_df = annotated_class.complete_df
        self.starting_df.TISSUE_FACTOR.fillna("UNKNOWN", inplace=True)
        self.working_df = self.starting_df
        self.load_filter_dicts() # loads in default dictionaries (decide if we should keep or load optionally)

    def reset_all_filters( self ):
        self.working_df = self.starting_df

    def load_filter_dicts( self ):
        import default_sets
        self.dicts = default_sets.DefaultDictionaries()

    def filter_by_type( self, type_list = None ):
        if type_list == None:
            print( "Please enter list of types to keep. Acceptable types include: 'GENE_VARIANT', 'AMPLIFICATION', 'DELETION', 'METHYLATION'." )
        else:
            type_list = [ x.upper() for x in type_list ]
            self.working_df = self.working_df.query( "TYPE in @type_list" )

    def filter_by_gene_list( self, gene_list = None ):
        if gene_list == None:
            print( "Please enter gene list. To use the GDSC1000 core genes, add argument gene_list = 'Core Genes'." )
            if self.dicts in globals():
                print( "Possible dictionary keys include: ", self.dicts.gene_dict.keys() )
        elif isinstance( gene_list, str ):
            gene_list = self.dicts.gene_dict[ gene_list ]
        else:
            pass
        gene_list = [ x.upper() for x in gene_list ]
        self.working_df = self.working_df.query( "GENE in @gene_list" )


    def filter_by_tissue_list( self, tissue_list = None ):
        if tissue_list == None:
            print( "Please enter list of tissues. To use all tissues, add argument tissue_list = 'Complete'" )
            if self.dicts in globals():
                print( "Possible dictionary keys include: ", self.dicts.tissue_dict.keys() )
        elif isinstance( tissue_list, str ):
            tissue_list = self.dicts.tissue_dict[ tissue_list ]
        else:
            pass
        tissue_list = [ x.upper() for x in tissue_list ]
        self.working_df = self.working_df.query( "TISSUE_FACTOR in @tissue_list" )

    def filter_by_cosmic_list( self, cosmic_list = None ):
        if cosmic_list == None:
            print( "Please enter list of COSMIC IDs. To use an example set of COSMIC IDs, add argument cosmic_list = 'Dummy'" )
            if self.dicts in globals():
                print( "Possible dictionary keys include: ", self.dicts.cosmic_dict.keys() )
        elif isinstance( cosmic_list, str ):
            cosmic_list = self.dicts.cosmic_dict[ cosmic_list ]
        else:
            pass
        self.working_df = self.working_df.query( "COSMIC_ID in @cosmic_list" )


    def filter_by_cell_line_list( self, cell_line_list = None ):
        if cell_line_list == None:
            print( "Please enter list of cell line names. To use an example set of names, add argument cell_line_list = 'Dummy'" )
            if self.dicts in globals():
                print( "Possible dictionary keys include: ", self.dicts.cell_line_dict.keys() )
        elif isinstance( cell_line_list, str ):
            cell_line_list = self.dicts.cell_line_dict[ cell_line_list ]
        else:
            pass
        self.working_df = self.working_df.query( "CELL_LINE in @cell_line_list" )



    def filter_by_recurrence( self, min_recurrence=3):
        if self.annotate == True:
            feature = "GENE"
        else:
            feature = "IDENTIFIER"

        flatten_cell_lines = self.working_df.groupby( [ feature, "COSMIC_ID", "TISSUE_FACTOR" ], as_index = False )[[ "TYPE" ]].count()
        feature_count = flatten_cell_lines.groupby( [feature], as_index = False )[[ "COSMIC_ID" ]].count()
        self.feature_count = feature_count[ feature_count["COSMIC_ID"] >= min_recurrence ]

        self.working_df = self.working_df[ self.working_df[ feature ].isin( self.feature_count[ feature ] ) ]



    # Make into matrix
    def make_matrix( self, min_recurrence=3):
        """Return a dataframe compatible with ANOVA analysis"""
        self.filter_by_recurrence( min_recurrence=3)

        if self.annotate == True:
            feature = "GENE"
        else:
            feature = "IDENTIFIER"

        final_matrix = pd.crosstab( self.working_df[ feature ], columns=[ self.working_df["COSMIC_ID"], self.working_df['TISSUE_FACTOR'], self.working_df["CELL_LINE"], self.working_df[ "MSI_FACTOR" ], self.working_df[ "MEDIA_FACTOR"] ] )
        final_matrix[ final_matrix > 1 ] = 1
        self.final_matrix = final_matrix.T.reset_index()
        self.final_matrix.to_csv( "./data/gdsc_1000_data/features_matrix.csv", index = False )
