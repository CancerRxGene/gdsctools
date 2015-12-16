ChangeLog
==============

.. contents::


Version 0.10
--------------------------

Lots of changes in this version but for users the API should be very similar.

.. rubric:: 0.10

* NEW:
    - Add a new factor called MEDIA_FACTOR. If not provided, genomic
      feature matrix can populated the MEDIA_FACTOR column automatically.
    - add a class COSMICInfo and a related data file called
      cosmic_info.csv.gz to get information about COSMIC ids. Replaces
      COSMIC class, which was removed.
    - add new class GDSC to perform the entire analysis splitting data across
      companies found in DrugDecode and across cancer types.


* CHANGES:
    - COSMIC class removed and replaced by COSMICInfo class
    - column name convention:
        - FEATURE_ANOVA_pval --> ANOVA_FEATURE_pval
        - MSI_ANOVA_pval --> ANOVA_MSI_pval
        - TISSUE_ANOVA_pval --> ANOVA_TISSUE_pval
        - FEATURE_ANOVA_FDR_% -->  ANOVA_FEATURE_FDR
        - new column named ANOVA_MEDIA_pval
        - to be constistent, names such as FEATURE_pos have now underscores
          to separate words e.g., (FEATUREpos --> FEATURE_pos, FEATUREneg 
          --> FEATURE_neg, deltaMEAN --> delta_MEAN).
    - refactor :mod:`gdsctools.volcano` module to use new naming convention.
    - SAMPLE_NAME is not required anymore in the genomic features. This is
      indeed just an annotation and is now encoded in the flat file
      cosmic_info.csv.gz (see above)
    - :mod:`~gdsctools.anova`, anova_results modules:
        - Implement new factor (MEDIA) in the regression
        - Uses new naming convention for the columns as described above
        - When initialising a ANOVA instance, prints the factor that will be
          included.
        - add new option (set_media_factor) to populate the MEDIA column
          automatically
    - :mod:`~gdsctools.readers` module:
        - 'Sample Name' or SAMPLE_NAME are deprecated.
          There are removed from the genomic_feature matrix if found.
    - Uses MEDIA_FACTOR column in addition to MSI and tissue columns
    - shift attribute is now read-only and set automatically
    - add a function to fill media column automatically
    - print function is  more verbose
    - volcano: uses new naming convention for the columns as described above.
    - split :mod:`~gdsctools.anova` module (create
      :mod:`~gdsctools.anova_report`) (issue #98).
    - :mod:`~gdsctools.readers`: improved DrugDecoder and renamed into
      DrugDecode (issue #102 and #101)
    - add new settings and code to apply pvalue correction at drug level
      rather than global level.
    - add new module to find chemblId/ChemSpider from drug name.

Version 0.9
--------------------------

.. rubric:: 0.9.10

* NEW:
   - add settings as json file in the HTML report
   - ANOVAResults has now a volcano() method
   - add read_settings method in ANOVA
   - add code in the HTML tree directory to reproduce HTML report and results

* CHANGES:
   - anova_one_drug now returns an ANOVAResults object
   - Restructure data package tree directory (#83)
   - Default header have changed:
       - COSMIC ID --> COSMID_ID
       - Sample Name --> SAMPLE_NAME
       - MS-instability Factor Value --> MSI_FACTOR
       - Tissue Factor Value --> TISSUE_FACTOR

     Previous values will still be accepted but deprecation warning added.

* BUG FIXES:
    - Fixes #89 (tight layout buggy under MAC)

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

Version 0.3
------------------------

- Cancer specific now included and tested on BRCA and BLCA cases.


Version 0.2
---------------

First working version with HTML output.

Version 0.1
---------------

First working version of gdsctools with test and documenation.
Tested against version17. A standalone app is also provide as a command
line argument (named **gdsctools_anova**).
