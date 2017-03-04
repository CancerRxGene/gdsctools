# -*- python -*-
# -*- coding utf-8 -*-
#
#  This file is part of GDSCTools software
#
#  Copyright (c) 2016 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): GDSCTools authors
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
"""OmniBEM functions"""
from gdsctools import Reader, GenomicFeatures

import pandas as pd
import pylab


__all__ = ["OmniBEMBuilder"]


class OmniBEMBuilder(object):
    """Utility to create :class:`~gdsctools.readers.GenomicFeatures` instance from BEM data

    Starting from an example provided in GDSCTools
    (test_omnibem_genomic_alteration.csv.gz), here is the code to get a data
    structure compatible with the :class:`GenomicFeature`, which can then be
    used as input to the :class:`~gdsctools.anova.ANOVA` class.

    See the constructor for the header format.

    .. plot::
        :include-source:

        from gdsctools import gdsctools_data, OmniBEMBuilder, GenomicFeatures

        input_data = gdsctools_data("test_omnibem_genomic_alterations.csv.gz")
        bem = OmniBEMBuilder(input_data)

        # You may filter the data for instance to keep only a set of genes.
        # Here, we keep everything
        gene_list = bem.df.GENE.unique()
        bem.filter_by_gene_list(gene_list)

        # Finally, create a MoBEM dataframe
        mobem = bem.get_mobem()

        # You may filter with other functions(e.g., to keep only Methylation
        # features
        bem.filter_by_type_list(["Methylation"])
        mobem = bem.get_mobem()

        # Then, let us create a dataframe that is compatible with
        # GenomicFeature. We just need to make sure the columns are correct
        mobem[[x for x in mobem.columns if x!="SAMPLE"]]
        gf = GenomicFeatures(mobem[[x for x in mobem.columns if x!="SAMPLE"]])


        # Now, we can perform an ANOVA analysis:
        from gdsctools import ANOVA, ic50_test
        an = ANOVA(ic50_test, gf)
        results = an.anova_all()
        results.volcano()

    .. note:: The underlying data is stored in the attribute :attr:`df`.

    """
    def __init__(self, genomic_alteration):
        """.. rubric:: Constructor

        :param str genomic_alteration: a filename in CSV format (gzipped or
            not). The format is explained here below.

        The input must be a 5-columns CSV file with the following columns::

            COSMIC_ID: an integer
            TISSUE_TYPE: e.g. Methylation,
            SAMPLE: this should correspond to the COSMID_ID
            TYPE: Methylation,
            GENE: gene name
            IDENTIFIER: required for now but may be removed (rows can be
                used as identifier indeed

        .. warning:: If GENE is set to NA, we drop it. In the resources shown in
            the example here above, this corresponds to 380 rows (out of 60703).
            Similarly, if TISSUE is NA, we also drop these rows that is about
            3000 rows.
        """
        self.df = Reader(genomic_alteration).df

        # We drop genes and tissues that are set to NA
        self.df.dropna(inplace=True)

        # Some column names are changed for convenience
        self.df.rename(columns={"COSMIC.ID": "COSMIC_ID",
            "TISSUE.TYPE": "TISSUE_TYPE"}, inplace=True)
        self._update_unified()

    def __len__(self):
        return len(self.df)

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

        The returned object can be read
        by :class:`~gdsctools.readers.GenomicFeatures`.

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

        if "TISSUE_TYPE" in this.columns:
            this.rename(columns={"TISSUE_TYPE":"TISSUE_FACTOR"}, inplace=True)
        else:
            print("Expected TISSUE_TYPE column. Not found.")

        return this

    def get_significant_genes(self, N=20):
        """Return most present genes

        :param int N: the maximum number of genes to return
        :return: list of genes with the number of occurences

        The genes returned by be used to filter the data::

            genes = bem.get_significant_genes(N=20)
            bem.filter_by_gene_list(genes.index)

        """
        df = self.unified.groupby("GENE")["GENE"]
        return df.count().sort_values(ascending=False).head(N)

    def filter_by_gene_list(self, genes, minimum_gene=3):
        """Keeps only required genes

        :param genes: a list of genes or a filename (a CSV file with a column
            named GENE).
        :param minimum_gene: genes with not enough entries are removed
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
        """Filter the data by tissue

        :param list tissue: the tissues to be kept. The data is update
            in place.

        """
        self.df = self.df.query("TISSUE_TYPE in @tissue_list")
        self._update_unified()

    def filter_by_type_list(self, type_list):
        """Filter the data by type

        :param list tissue: the types to be kept. The data is update
            in place. Here are some examples of types: Point.mutation,
            Amplification, Deletion, Methylation.

        """
        self.df = self.df.query("TYPE in @type_list")
        self._update_unified()

    def filter_by_sample_list(self, sample_list):
        """Filter the data by sample name

        :param list tissue: the samples to be kept. The data is updated
            in place.

        """
        self.df = self.df.query("SAMPLE in @sample_list")
        self._update_unified()

    def filter_by_cosmic_list(self, cosmic_list):
        """Filter the data by cosmic identifers

        :param list cosmic: the cosmic identifiers to be kept. The data is
            updated in place.

        """
        self.df = self.df.query("COSMIC_ID in @cosmic_list")
        self._update_unified()

    def plot_number_alteration_by_tissue(self, fontsize=10, width=0.9):
        """Plot number of alterations

        .. plot::
            :width: 100%
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
        count.plot(kind="bar", width=width)
        pylab.grid()
        pylab.xlabel("Tissue Type", fontsize=fontsize)
        pylab.ylabel("Total number of alterations in cell lines",
                     fontsize=fontsize)
        try:pylab.tight_layout()
        except:pass

    def plot_alterations_per_cellline(self, fontsize=10, width=0.9):
        """Plot number of alterations

        .. plot::
            :width: 100%
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
        df.plot(y="fraction", legend=False, kind='bar', fontsize=fontsize, 
            width=width)
        pylab.ylabel("Alterations per cell line",
                     fontsize=fontsize)
        pylab.grid()
        try:pylab.tight_layout()
        except:pass

    def get_genomic_features(self):
        mobem = self.get_mobem()
        gf = GenomicFeatures(mobem[[x for x in mobem.columns if x!="SAMPLE"]])
        return gf
