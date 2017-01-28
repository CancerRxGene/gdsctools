
.. index:: gdsctools_anova
.. _standalone:

Standalone application
==========================

Although we would encourage you to use the Python shell to have as much
flexibility as possible, we also provide a standalone applications.

Currently, there are two standalones. The **gdsctools_anova** and
**gdsctools_regression**. The first one is a pure Python implementation while
the second one is snakemake-based.


gdsctools_anova application
---------------------------------


 called **gdsctools_anova**. This standalone application should be installed with **GDSCTools** automatically. It focuses on the ANOVA analysis only, and can be used to analysis one set of IC50 and genomic feature at a time.


You can obtain the help by typing::

    gdsctools_anova --help


The main goal is to provide an interface to the python library and consequently, one be able to redo the analysis as shown in the quickstart::


* One drug One Feature with figure(s) and HTMLs
* One Drug All Feature with figure and HTMLs
* All Drug All Feature with figures and HTMLs

We suppose the input data file is called IC50_10drugs.tsv

:term:`ODOF`
~~~~~~~~~~~~~~~~~~

::

    gdsctools_anova --input-ic50 IC50_10drugs.tsv --drug
        Drug_999_IC50 --feature TP53_mut --onweb


:term:`ODAF`
~~~~~~~~~~~~~~~~~~
::

    gdsctools_anova --input-ic50 IC50_10drugs.tsv --drug
        Drug_999_IC50 --onweb



:term:`ADAF`
~~~~~~~~~~~~~~~~~~

::

    gdsctools_anova --input-ic50 IC50_10drugs.tsv --onweb



Some other settings
~~~~~~~~~~~~~~~~~~~~~~~~~~


Again, you can use the ``--help`` to get up-to-date information about the available
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


gdsctools_regression application
---------------------------------

Let us consider the case where you have an IC50 file and a genomic file. The
first step consists in preparing the working directory::

    gdsctools_regression -I IC50_v17.csv.gz -F genomic_features_v17.csv.gz
        --method lasso -O lasso_analysis

::

    cd lasso_analysis

On a local computer::

    snakemake -s regression.rules -j 4

On a distributed-computing system using e.g SLURM framework, use::

    srun --qos normal snakemake -s regression.rules -j 40 --cluster "sbatch --qos normal"




























