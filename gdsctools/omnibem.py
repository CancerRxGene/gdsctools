import pandas as pd
from gdsctools import Reader
import pylab


__all__ = ["OmniBEMBuilder"]



class OmniBEMBuilder(object):
    """Utility to create :class:`~gdsctools.readers.GenomicFeatures` instance from BEM data

    Starting from an example provided in GDSCTools
    (test_omnibem_genomic_alteration.csv.gz), here is the code to get a data
    structure compatible with the :class:`GenomicFeature`, which can then be
    used as input to the :class:`~gdsctools.anova.ANOVA` class.

    ::

        from gdsctools import gdsctools_data, OmniBEM, GenomicFeatures

        input_data = gdsctools_data("test_omnibem_genomic_alteration.csv.gz")
        bem = OmniBEM(input_data)

        # You may filter the data for instance to keep only a set of genes
        bem.filter_by_gene_list(gene_list)

        # Finally, create a MoBEM dataframe
        mobem = bem.get_mobem()

        # You may filter with other functions(e.g., to keep only Methylation
        # features
        bem.filter_by_type_list(["Methylation"])
        mobem = bem.get_mobem()

        # Then, let us create a dataframe that is compatible with
        # GenomicFeature. We just need to make sure the columns are correct
        mobem.columns = [x.replace("TISSUE_TYPE", "TISSUE_FACTOR") for x in mobem.columns]
        mobem[[x for x in mobem.columns if x!="SAMPLE"]]
        gf = GenomicFeatures(mobem[[x for x in mobem.columns if x!="SAMPLE"]])


        # Now, we can perform an ANOVA analysis:
        an = ANOVA(ic50_test, gf)
        results = an.anova_all()
        results.volcano()

    .. note:: The underlying data is stored in the attribute :attr:`df`.

    """
    def __init__(self, genomic_alteration):
        """.. rubric:: Constructor


        The input must be a CSV file (gzipped or not) with the following
        columns::

            COSMIC_ID:
            TISSUE_TYPE: e.g. Methylation,
            SAMPLE: this should correspond to the COSMID_ID
            GENE:
            IDENTIFIER: required but may be removed (rows can be used as
                identifier indeed


        """
        self.df = Reader(genomic_alteration).df
        # Note that if we keep NA, the groupby will not work as expected.
        self.df.fillna("", inplace=True)
        self.df.rename(columns={"COSMIC.ID": "COSMIC_ID",
            "TISSUE.TYPE": "TISSUE_TYPE"}, inplace=True)
        self._update_unified()

    def _update_unified(self):
        """# Sort genes by number of times they have an alteration
        bem.unified.sort_values(by="IDENTIFIER")

        # top genes represented multiple times i
        bem.unified.sort_values(by="IDENTIFIER", ascending=False).head(20)
        """
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

    def get_mobem(self):
        """Return a dataframe compatible with ANOVA analysis



        """
        # Select gene that appear at least a minimum number of times
        #agg = self.unified.groupby("GENE")["GENE"].count()
        #self.selection = agg[agg>=minimum_gene]

        # keep only gene in the selection
        #df = self.unified.query("GENE in @self.selection.index")
        df = self.unified
        this = pd.crosstab(df['GENE'], columns=[
            df["COSMIC_ID"], df['TISSUE_TYPE'], df["SAMPLE"]])
        this = this.T
        this = this.reset_index()
        return this

    def get_significant_genes(self, N=20):
        df = self.unified.groupby("GENE")["GENE"]
        return df.count().sort_values(ascending=False).head(N)

    def filter_by_gene_list(self, genes, minimum_gene=3):
        """Keeps only required genes

        :param genes: a list of genes or a filename (a CSV file with a column
            named GENE).
        :param minimum_gene: genes with not enough entries are remove
            (defaults to 3)

        """
        if isinstance(genes, str):
            df = pd.read_csv(genes)
            genes = df.GENE.values.tolist()
        self.df = self.df.query("GENE in @genes")
        # Update unified matrix
        agg = self.unified.groupby("GENE")["GENE"].count()
        self.selection = agg[agg>=minimum_gene]

        # keep only gene in the selection
        self.unified = self.unified.query("GENE in @self.selection.index")
        #self._update_unified()

    def filter_by_tissue_list(self, tissue_list):
        self.df = self.df.query("TISSUE_TYPE in @tissue_list")
        self._update_unified()

    def filter_by_type_list(self, type_list):
        """Keeps only required type"""
        self.df = self.df.query("TYPE in @type_list")
        self._update_unified()

    def filter_by_sample_list(self, sample_list):
        """Filter data to keep sample in sample_list"""
        self.df = self.df.query("SAMPLE in @sample_list")
        self._update_unified()

    def filter_by_cosmic_list(self, cosmic_list):
        """Filter data to keep cosmic ids in cosmic_list"""
        self.df = self.df.query("COSMIC_ID in @cosmic_list")
        self._update_unified()

    def plot_number_alteration_by_tissue(self):
        """Plot number of alterations

        .. plot::
            :width: 50%
            :include-source:

            from gdsctools import *
            data = gdsctools_data("test_omnibem_genomic_alterations.csv.gz")
            bem = OmniBEMBuilder(data)
            bem.filter_by_gene_list(gdsctools_data("test_omnibem_genes.txt"))
            bem.plot_number_alteration_by_tissue()

        """
        count = self.unified.groupby(['TISSUE_TYPE'])['GENE'].count()
        try:
            count.sort_values(inplace=True, ascending=False)
        except:
            count.sort(inplace=True, ascending=False)
        count.plot(kind="bar")
        pylab.grid()
        pylab.xlabel("Tissue Type")
        pylab.ylabel("Total number of alterations in cell lines")
        try:pylab.tight_layout()
        except:pass

    def plot_alterations_per_cellline(self):
        """Plot number of alterations

        .. plot::
            :width: 50%
            :include-source:

            from gdsctools import *
            data = gdsctools_data("test_omnibem_genomic_alterations.csv.gz")
            bem = OmniBEMBuilder(data)
            bem.filter_by_gene_list(gdsctools_data("test_omnibem_genes.txt"))
            bem.plot_alterations_per_cellline()
        """
        try:
            df = self.summary.sort_values("fraction", ascending=False)
        except:
            df = self.summary.sort("fraction", ascending=False)
        df.plot(y="fraction", legend=False, kind='bar', fontsize=10, width=0.8)
        pylab.ylabel("Alterations per cell line")
        pylab.grid()
        try:pylab.tight_layout()
        except:pass




