from gdsctools.anova import ANOVA


class GDSC(ANOVA):
    """Alias to ANOVA class with default settings 
    
    It also converts tissue names into TCGA names.
    """



    def __init__(self, ic50, genomic_features,
        drug_decode, verbose=True, low_memory=False):
        

        super(GDSC, self).__init__(ic50, genomic_features=genomic_features, 
            drug_decoder=drug_decode, verbose=True, low_memory=low_memory)
        from gdsctools.tissues import map_to_tcga
        def tcga_names(x):
            if x in map_to_tcga.keys():
                return map_to_tcga[x]
            else:
                return x
        newnames =  self.features.df.TISSUE_FACTOR.apply(
                            lambda x: tcga_names(x))

        self.features.df['TISSUE_FACTOR'] = newnames
