
.. _omnibem:

OmniBEM Builder
================

**OmniBEM Builder** is an optional tool within **GDSCTools** designed to give the user more flexibility on the levels of genomic annotation probed by **GDSCTools**.

By default, **GDSCTools** is based on the genomic annotation of 1001 cell lines represented in COSMIC and published by *Iorio et al* (Cell Resource). The annotation includes genetic variants as identified by exome sequencing, copy number alterations and differentially methylated CPG islands. 

OmniBEM Builder allows the user to merge the different levels of annotations into a single gene-level view that queries whether a given gene has been altered at any level of annotations.

The user can additionally specify which sets of genomic annotations to integrate as well as upload and integrate their own sets of genomic annotations.

::

    from gdsctools import OmniBEM


Input Data Structure
----------------------

At its most basic, OmniBEM Builder requires a cell line ID (e.g. COSMIC ID), a gene name and an alteration type (e.g. Point mutation, Amplification or Copy Number Alteration). A simple example table is given below:

==========  ==========  ============ ======== ============= =============
COSMIC ID   GENE        TYPE          SAMPLE   TISSUE_TYPE    IDENTIFIER
==========  ==========  ============ ======== ============= =============
111         ZAP         Methylation    NM        breast           1
==========  ==========  ============ ======== ============= =============

Example
----------

We provide a data set available in **GDSCTools** that can be loaded as follows :: 

    from gdsctools import gdsctools_data, OmniBEMBuilder

    input_data = gdsctools_data("test_omnibem_genomic_alterations.csv.gz")
    bem = OmniBEMBuilder(input_data)

The data is stored in the :attr:`df` attribute::

    bem.df

From this data frame, one can filter, group and perform various data mangling
operations. For instance these commands group the data by tissue type, count the
number of row per tissue and return the 5 most representative tissues::

    >>> count = bem.df.groupby("TISSUE_TYPE").count()
    >>> list(count.sort_values("COSMIC_ID").index[-5:])
    ['SKCM', 'BRCA', 'HNSC', 'LUAD', 'COAD/READ']

In :class:`~gdsctools.omnibem.OmniBEMBuilder`, we provide convenient methods to
filter the data by genes, cosmic IDs, sample names, tissues and types.

For instance::

    >>> # You may filter the data for instance to keep only a set of genes.
    >>> # Here, we keep the 100 most present genes:
    >>> bem = OmniBEMBuilder(input_data)
    >>> len(bem)
    56943 
    >>> gene_list = bem.get_significant_genes(100)
    >>> bem.filter_by_gene_list(gene_list.index)  # gene_list is a dataframe
    >>> len(bem)
    6141

Once you BEM data is filtered as expected, save it::

    gf = bem.get_genomic_features()
    gf.to_csv("Your_genomic_features.csv")

This file can now be used with tha ANOVA or Regression analysis.


.. seealso:: :class:`gdsctools.omnibem.OmniBEMBuilder`



