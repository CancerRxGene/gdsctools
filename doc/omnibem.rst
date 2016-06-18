OmniBEM Builder
================

**OmniBEM Builder** is an optional tool within gdsctools designed to give the user more flexibility on the levels of genomic annotation probed by gdsctools.

By default, gdsctools is based on the genomic annotation of 1001 cell lines represented in COSMIC and published by *Iorio et al* (Cell Resource, in press). The annotation includes genetic variants as identified by exome sequencing, copy number alterations and differentially methylated cpg islands. 

OmniBEM Builder allows the user merge the different levels of annotation into a single gene-level view that queries whether a given gene has been altereted at any level of annotation.

The user can additionally specify which sets of genomic annotations to integrate as well as upload and integrate their own sets of genomic annotation.

::

    from gdsctools import OmniBEM


Input Data Structure
----------------------

At its most basic, OmniBEM Builder requires a cell line ID (e.g. COSMIC ID), a gene name and an alteration type (e.g. Point mutation, Amplification or Copy Number Alteration). A simple example table is given below:

==========  ==========  ==========
COSMIC ID   GENE        TYPE
==========  ==========  ==========
