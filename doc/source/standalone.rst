
.. index:: gdsctools_anova
.. _standalone:

Standalone applications
==========================

The standalone application is called **gdsctools_anova**. You can obtain the
help typing::

    gdsctools_anova --help


The main goal is to provide an interface to the python library. One should be
able to redo the analysis shwon in the quickstart.


* ODOF with figure(s)
* ODAF with figure
* ADAF with figure
* All analysis with HTML report

We supposer the input data file is called IC50.tsv

ODOF
-----------

::

    gdsctools_anova --input-ic50 ../share/data/IC50_10drugs.tsv --drug
    Drug_999_IC50 --feature TP53_mut --onweb


ODAF
----------


ADAF
---------




other settings
------------------
