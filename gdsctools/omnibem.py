import pandas as pd
from gdsctools import Reader
import pylab


__all__ = ["OmniBEMBuilder"]



class OmniBEMBuilder(object):
    """Create a BEM (Genomic Feature) file

    ::

        bem = OmniBEM("data.csv")
        # optional
        bem.filter_by_gene_list(gene_list)
        mobem = bem.get_mobem()

        # Sort genes by number of times they have an alteration
        bem.unified.sort_values(by="IDENTIFIER")


        # top genes represented multiple times i
        bem.unified.sort_values(by="IDENTIFIER", ascending=False).head(20)


    You can filter by, gene, cosmic identifiers, sample name, type of
    alterations.


    - merge features or not.
    """
    def __init__(self, genomic_alteration):
        self.df = Reader(genomic_alteration).df
        # Note that if we keep NA, the groupby will not work as expected.
        self.df.fillna("", inplace=True)
        self.df.rename(columns={"COSMIC.ID": "COSMIC_ID",
            "TISSUE.TYPE": "TISSUE_TYPE"}, inplace=True)
        self._update_unified()

    def _update_unified(self):
        # In R this is an aggregate function. In Pandas a groupby + aggregate
        # http://pandas.pydata.org/pandas-docs/stable/comparison_with_r.html
        groups = self.df.groupby(by=["COSMIC_ID", "GENE", "SAMPLE",
            "TISSUE_TYPE"], as_index=False)
        self.unified = groups[['IDENTIFIER']].aggregate(len)

        # Building a summary on cell lines
        unique = self.unified[self.unified['COSMIC_ID'].duplicated() == False]
        self.frequency = unique.groupby("TISSUE_TYPE")["TISSUE_TYPE"].count()
        self.total = self.unified.groupby("TISSUE_TYPE")["TISSUE_TYPE"].count()

        df = pd.concat([self.total, self.frequency, self.total/self.frequency], axis=1)
        df.columns = ['total', 'grouped', 'fraction']
        self.summary = df

    def get_mobem(self, minimum_gene=3):
        # Select gene that appear at least a minimum number of times
        agg = self.unified.groupby("GENE")["GENE"].count()
        self.selection = agg[agg>=minimum_gene]

        # keep only gene in the selection
        df = self.unified.query("GENE in @self.selection.index")

        this = pd.crosstab(df['GENE'], columns=[
            df["COSMIC_ID"], df['TISSUE_TYPE'], df["SAMPLE"]])
        this = this.T
        this = this.reset_index()
        return this

    def get_significant_genes(self, N=20):
        df = self.unified.groupby("GENE")["GENE"]
        return df.count().sort_values(ascending=False).head(N)

    def filter_by_gene_list(self, genes):
        """

        :param genes: a list of genes or a filename


        """
        if isinstance(genes, str):
            df = pd.read_csv(genes)
            genes = df.GENE.values.tolist()
        self.df = self.df.query("GENE in @genes")
        # Update unified matrix
        self._update_unified()

    def filter_by_tissue_list(self, tissue_list):
        self.df = self.df.query("TISSUE_TYPE in @tissue_list")
        self._update_unified()

    def filter_by_type_list(self, type_list):
        self.df = self.df.query("TYPE in @type_list")
        self._update_unified()

    def filter_by_sample_list(self, sample_list):
        self.df = self.df.query("SAMPLE in @sample_list")
        self._update_unified()

    def filter_by_cosmic_list(self, cosmic_list):
        self.df = self.df.query("COSMIC_ID in @cosmic_list")
        self._update_unified()

    def plot_number_alteration_by_tissue(self):
        count = self.unified.groupby(['TISSUE_TYPE'])['GENE'].count()
        count.sort_values(inplace=True, ascending=False)
        count.plot(kind="bar")
        pylab.grid()
        pylab.xlabel("Tissue Type")
        pylab.ylabel("Total number of alterations in cell lines")
        try:pylab.tight_layout()
        except:pass

    def plot_alterations_per_cellline(self):
        df = self.summary.sort_values("fraction", ascending=False)
        df.plot(y="fraction", legend=False, kind='bar', fontsize=10, width=0.8)
        pylab.ylabel("Alterations per cell line")
        try:pylab.tight_layout()
        except:pass




