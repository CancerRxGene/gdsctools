OmniBEM Builder
================

**OmniBEM Builder** is an optional tool within **GDSCtools** designed to give the user more flexibility on the levels of genomic annotation probed by **GDSCtools**.

By default, **GDSCtools** is based on the genomic annotation of 1001 cell lines represented in COSMIC and published by *Iorio et al* (Cell Resource). The annotation includes genetic variants as identified by exome sequencing, copy number alterations and differentially methylated CPG islands. 

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


.. seealso:: :class:`gdsctools.omnibem.OmniBEMBuilder`
