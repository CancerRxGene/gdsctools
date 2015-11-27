ChangeLog
==============



version 0.9 Nov 2015
--------------------------

.. rubric:: 0.9.10

* NEW:
  - add settings as json file in the HTML report

* CHANGES:

* BUG FIXES:


.. rubric:: 0.9.9

* CHANGES:
   - add new regression method: Ridge/Lasso/ElasticNet in
     :class:`gdsctools.anova.ANOVA`
   - Rename some of the settings to have a more uniform naming convention in
     :class:`gdsctools.settings.ANOVASettings`
   - Add new module related to fitting ot logistic function  parameters
     (:mod:`gdsctools.logistics`)

.. rubric:: 0.9.8

* BUG: javascript were not included in version 0.9.7 had to rename js directory
  into javascript to avoid known bug in distutils. Maybe solved in the future
  but for bow just renamed the directory.

.. rubric:: 0.9.7

- MSI/Sample/Tissue columns in the genomic features are not required anymore.
- FDR lines in volcano plots are now using interpolation and 
  therefore more  precisily placed. Fixes #57
- volcano plot improvments. Fixes #79, #80, #81
- Fixes issue #72 to get the drug_decoder information from the ANOVA class.
- Fixes issue #76  to drop IC50 cosmic Id not found in the genomic feature
  matrix
- Readers (e.g. IC50) can now read CSV files with commented lines (# character)
  issue #78
- Readers can now ignored columns that are not named (usually first column of
  index exported by excel document)
- IC reader figure out automatically if the prefix "Drug" has been used. It so, 
  it drops other irrelevant columns. Useful if genomic features and IC50 are
  mixed together.
- IC50 and GenomicFeatures, DrugDecode now accepts both TSV and CSV format
  (gziped or not)
- add more datasets for testing purposes
- double checked results on BLCA tissue v17 and v18
- Finalise a first version of the standalone application 
- ReadTheDocs documentation is now on line gdsctools.readthedocs.org
- GDSCTools has now all features of the original R version
- With in addition:
  - a standalone application
  - test suite
  - documentation
- benchmarking for the analysis in about 20 minutes 265 drugs and 680 features
  across 980 cell lines. HTML report takes as much time. 

Before version 0.9
------------------------

v0.3
~~~~
- Cancer specific now included and tested on BRCA and BLCA cases.


v0.2 23 Oct 2015
~~~~~~~~~~~~~~~~~~~~

First working version with HTML output.

v0.1 14 Oct 2015
~~~~~~~~~~~~~~~~~~~~~

First working version of gdsctools with test and documenation. 
Tested against version17. A standalone app is also provide as a command
line argument (named **gdsctools_anova**).
