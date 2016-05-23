import pandas as pd


class OmniBEM(object):
    def __init__(self, genomic_alteration):
        self.df = pd.read_csv(genomic_alteration, sep=",")


    def filter(self, gene_list, minimum_gene=3):
        
        pass 
        """self.df
        # filtering happening here
        results = ... 

        return results
        """

