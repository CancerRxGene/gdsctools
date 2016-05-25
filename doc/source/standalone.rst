
.. index:: gdsctools_anova
.. _standalone:

Standalone applications
==========================

Although we would encourage you to use the Python shell to have as much
flexibility as possible, we also provide a standalone application called **gdsctools_anova**. This standalone application should be installed with **GDSCTools** automatically. It focuses on the ANOVA analysis only, and can be used to analysis one set of IC50 and genomic feature at a time. 


You can obtain the help by typing::

    gdsctools_anova --help


The main goal is to provide an interface to the python library and consequently, one be able to redo the analysis as shown in the quickstart::


* One drug One Feature with figure(s) and HTMLs
* One Drug All Feature with figure and HTMLs
* All Drug All Feature with figures and HTMLs

We suppose the input data file is called IC50_10drugs.tsv

ODOF
-----------

::

    gdsctools_anova --input-ic50 IC50_10drugs.tsv --drug
        Drug_999_IC50 --feature TP53_mut --onweb


ODAF
----------
::

    gdsctools_anova --input-ic50 IC50_10drugs.tsv --drug
        Drug_999_IC50 --onweb



ADAF
---------

::

    gdsctools_anova --input-ic50 IC50_10drugs.tsv --onweb



Some other settings
----------------------


Again, you can use the --help to get up-to-date information about the available
arguments. However, let us give a couple of interesting ones.

If you are interesting in a specific association of drug and feature, it is
convenient to get the valid drug names::

    --print-drug-names

or feature names::
    
    --print-feature-names

By default the analysis is :term:`PANCAN` (includes all tissues) but you can restrict the analysis to a set of tissues (or just one)::
    
    --tissues breast, cervix 

To know the names of the tissues, use::

    --print-tissue-names






