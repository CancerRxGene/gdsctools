ChangeLog
==============

.. contents::


Version 0.17.1
---------------

* NEWS:

    - add test for module gdsc1000, which has been refactored and cleanup

* BUGS: 

    - Fix GDSC links to fetch drug/genomic feature data sets
    - Fix bug in regression.rules pipeline (create missing directory) 

Version 0.17
------------
:summary: Fixing bunch of deprecated warnings, working regression 
    pipeline based on snakemake

* CHANGES and BUG Fixes:

     - regression (snakemake pipeline)
     - add missing init and regression modules
     - starting to use colorlog

* NEWS:

    - add regression script to create the snakefile workflow easily


Version 0.16.4
-----------------

* NEWS:

    - a snakemake pipeline for the regression analysis
    - Report for the regression analysis

Version 0.16.3
----------------

* BUG Fix: #158 infinity are removed in ANOVAReport in the constructor

Version 0.16.2
----------------

* BUG Fix:

    - Fix volcano plot and anova results (avoid infinite values)
    

Version 0.16.1
---------------------

* BUG Fixes

    - Fix regression dendogram plot
    - Fix bug leading to NA in effect size reported by Carlos P. (private communication)"

* CHANGES:

    - Omnibem: drops the NA instead of replacing by ""


Version 0.16 (Nov 2016)
----------------------------
* BUG Fixes

    - Fix pending issue related to #153 
      ic_input.csv", "test_gf_input.csv" for the case where some media factor
      are missing. This is now handle properly.
    - Fix #151 : large integer are not cast properly with consequence that
      indices are strings, not integers leading to further issues in the 
      HTML pages. 
 
- NEWS:

     - anova_one_drug_one_feature_custom allows to perform any regression using
       formula like in R. This is not for production but should be useful to
       perform custom analysis
     - Include the MEDIA factor boxplot in the library and reports


Version 0.15 (Nov 2016)
-------------------------------------

* BUG Fixes

    - Fix issue #156 (GDSC failures in some cases). This was due to special
      MOBEM input files for which no tests are performed. In such cases, some
      codes were failing in ANOVAReport and ANOVAResults, which have been fixed
      in this release.
    - Fix buggy volcano plot (no plot if no data); if FDR threshold below
      minimum value, set fdr_threshold to minimum value so that it scaled the
      plots properly. Now, users can add as many lines as desired using
      settings.additional_fdr. pvalues found to be NAN are set to 0 to prevent
      plotting issues.
    - anova module: regression in the case Y ~ C(TISSUE) + C(MSI) + feature,
        the tissue sum of squares was using N-1 tissues (one missing).
    - anova module: regression in the case Y ~ C(TISSUE) + C(MEDIA) + C(MSI) + feature,
        the media sum of squares was not normalised properly.

* CHANGES:

   - elastic_net module renamed into regression
   - ANOVASetting prints keys in alphabetical order (instead of randomly)
   - ANOVASetting: regression_formula is added; other regression_XX settings
     are not removed but not used anymore.
   - GenomicFeature reader: if a tissue is empty, it is replaced by UNDEFINED.
   - standalone: the drug option must be an integer. THis is now caught are the
     option level, not later in the code.
   - anova module: remove code related to elastic net, ridge, and lasso. This
     won't be used in production with ANOVA. EN, Ridge and Lasso are used 
     in the regression module and will be part of an independent type of
     analysis. See NEWS

* NEWS:

   - Add Ridge + Lasso + LassoLars classes in addition to ElasticNet regression
     method into the regression module.
   - Add more features in regression module (boxplot, dendogram)
   - New regression notebook in the notebooks directory
   - anova module: We can now use any combo of regression formula using
     statsmodels. This is slower but one can do use any formula accepted
     by statsmodels. The previous faster code is still used for the standard
     analysis.


Version 0.14 (20th June 2016)
---------------------------------

* NEWS:

    - ElasticNet: new method elastic_all()
    - plot_elastic_weight in the gallery

* CHANGES:

    - ElasticNet plot_weights is now split into plot_weights 
      and plot_importance.
   
* BUGS: 

    - Fixes missing files in the pypi distributino (MANIFEST changed) 


Version 0.13 (27th May 2016)
-------------------------------

.. rubric:: 0.13.1

* CHANGES:

    - **DrugDecode**: In brief, the DRUG ID in the IC50 input file and the
      DrugDecode files should be integers. Some old data sets use the
      following convention to refer to a drug Drug_<ID>_IC50 and so DrugDecode
      was using the same convention. However, we now convert this type of
      identifier into integers. This is done internally for the IC50 file,
      however, was not done inside the DrugDecoder file. This is now effective.
    - HTML reports when using the GDSC class:
      - Company names now appear systematically in the top of the company data
        packages.
      - Drug Names were missing and do now appear in top of the relevant HTML
        pages.
    - Boxplots: If a DrugDecode file is provided Boxplots show the DRUG ID 
      and the real drug name in the matplotlib and JS boxplots


.. rubric:: 0.13.0

* CHANGES:

    - Reader class simplification and improvments: files can now be compressed
      using gzip but also xz, zip and bz2 formats. The NA can be encoded as NA
      or NaN strings. Spaces are interpreted as NA.
    - Sort DrugDecode's dataframe columns
    - Updated all documentation

* BUG:

      - Fix scaling of the data with newest version of scikit-learn
      - fix typo in the setup.py file. Passed travis + all tests before main
        release.

Version 0.12 (9th May 2016)
-------------------------------


.. rubric:: 0.12.1

* BUG:

    - add missing CSS in the distribution


.. rubric:: 0.12

* CHANGES:

    - SPEEDUP:
      - tissue specific analysis computational time decreased by 50%
        by dropping the creation of dataframe and using a simple numpy array
        inside ANOVA.anova_one_drug_one_feature
      - Creation of volcano plots uses pure javascript for the data packages
        and the creation of the volcano plots was dramatically sped up by a
        factor between 10 and 100e. One can still create volcano plot manually
        in pure matplotlib.
      - Similarly, boxplots for tissue, MSI and all associations are now
        created using JS.
    - Data packages have been refactored. The major difference concerns
      the HTML layout (most HTML files are now in the sub-directory
      called associations) so that is it cleaner at the top level. The volcano
      plots are not in PNG format anymore but pure HTML/JS, which can be
      exported manually. The consequences is that the creation of data
      packages is 10 times faster.
    - The standalone application had 2 options removed: --feature (alone)
      and --fast options
    - Drug Identifier are now handled as pure integer. For back
      compatibility, old files that mix up IC50 and Genomic Features (e.g. v17
      data) are still interpreted; the DRUG ID in that case are written as
      Drug_ID_IC50 and are transformed as just <ID> everywhere.
    - associations output were named 1.html, 2.html... and are now named
      a1.html, a2.html...
    - Because DRUG_ID are now integer and all HTML stored in the same directory
      the naming of the HTML files have been altered (e.g., associations starts
    - Report now accepts only one argument (the anova isntance). Second
      argument (results) is now optional. If not provided, ANOVA are computed on
      the fly
    - Multicore module removed but ANOVA.anova_all has multicore option. This
      seems to work on Linux systems. Not tested on windows or MacOsX
    - IC50 may have duplicated drug ids (at different concentrations). Not good
      practice but that the format of e.g. v18, v19 IC50 files. A class
      IC50Cluster was created to interepret those files. ANOVA will switch to
      IC50Cluster automatically if there are duplicated files.
    - Settings: low_memory option has been removed


Version 0.11 (April-May 2016)
--------------------------------

.. rubric:: 0.11.3

* CHANGES:

    - The parameter **pvalue_threshold** in the general settings was changed
      from infinite to 10e-3. This has an effect on the numlber of significant
      hits reported in the HTML reports and volvano plots. This should not have
      a strong impact on the number of hits but guarantees a reasonably low
      pvalue before multiple testing
    - If an input file named with .csv extension but the content is tabulated,
      there was no immediate error but lead to errors later (e.g. in ANOVA), which
      is difficult to debug. Now, in such cases, an error will occur immediately
      when reading the file.
    - The warnings about MEDIA factor is removed since most of the files do not
      contain that column.

* BUG

    - The data packages were stored in the "ALL" directory, which may be a  TCGA
      tissue by itself. This has been renamed into "tissue_packages".

.. rubric:: 0.11.2

* BUG:

    - add missing file in the setup.py

.. rubric:: 0.11.1

* BUG:

    - Fixes the missing data package in the setup for pip installation

.. rubric:: 0.11.0

* NEWS:

    - Elastic notebook and module implemented
    - GenomicFeatures has now a compression method

* CHANGES:

    - anova module was split into modules + anova so that elastic_net
      module can inherit from module
    - all share/data moved to gdsctools data
    - add scikit-learn dependencies

* BUGS:

    - Fix onevent picking in the volcano plot and use 4 digit for the FDR plot




Version 0.10
--------------------------

.. rubric:: 0.10.2

* BUG:

    - Fixes issue #127 (If MSI factor missing, the anova still tries to use it)
    - Fixes issue #126 (--out-directory ignored in gdsctools-anova pipeline)
    - Fixes issue #125 and #124 (HTML report links broken)

.. rubric:: 0.10.1

* BUG:

    - Fix set_cancer_type to accept lists of tissues again

* CHANGES:

    - Fixes #119 by adding more tests.
    - reactivate get_significant hits functions.
    - rename ANOVAResults.get_significant_hits into get_html_table

.. rubric:: 0.10

Lots of changes in this version but for users the API should be very similar.

* NEWS:

    - Add a new factor called MEDIA_FACTOR. If not provided, genomic
      feature matrix can populated the MEDIA_FACTOR column automatically.
    - add a class COSMICInfo and a related data file called
      cosmic_info.csv.gz to get information about COSMIC ids. Replaces
      COSMIC class, which was removed.
    - add new class GDSC to perform the entire analysis splitting data across
      companies found in DrugDecode and across cancer types.

* CHANGES:

    - COSMIC class removed and replaced by COSMICInfo class
    - Column name convention:
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

* BUGS:

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

* BUG:

    - javascript were not included in version 0.9.7 had to rename js directory
      into javascript to avoid known bug in distutils. Maybe solved in the
      future but for bow just renamed the directory.

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
