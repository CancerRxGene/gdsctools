.. _data_packages:

Data Packages
=================

.. warning:: Only one IC50 files should be provided. It is be cut 
    according to the GF file.

Definition
--------------

A **Data Package** is a terminology used to speak about a directory that
contains the results of an analysis. For example, the data package called
**BLCA** has a tree directory that looks like::

    BLCA/
    |-- code
    |-- css
    |-- images
    |-- INPUT
    |-- js
    `-- OUTPUT

The main directory contains a file named **index.html**. Many other HTML files
may be present but the **index** is your entry point to browse the content of
data package.


THe directory **css** and **js** contains resources required by the HTML
documents.

The **INPUT** directory contains 3 data files and the settings used during the
analysis:

#. ANOVA_input.csv
#. DRUG_DECODE.csv
#. genomic_features.csv
#. settings.json

See section on :ref:`data` for details about the data format and the
:ref:`settings` section.

Finally, the **OUTPUT** directory contains::

- drugs_summary.csv
- features_summary.csv
- results.csv

The **images** directory contains a mix of images and HTML for each
significant associations.

Create your own package
-----------------------------

In fact, we have already seen how to create a pacakge. This is covered in
:ref:`html` when we used the :class:`ANOVAReport` class but let us look at
the code again::

    from gdsctools import ANOVA, ic50_test, ANOVAReport
    gdsc = ANOVA(ic50_test)
    results = gdsc.anova_all()

    report = ANOVAReport(gdsc, results)
    report.create_html_pages()

Here, we have not yet mentionned the type of cancers or tissues since we used a
simple genomic feature file but one we need to repeat this analysis across man
y different genomic features files. The :class:`gdsctools.gdsc.GDSC` class will
help us for that purpose.

Create data packages across TCGA
--------------------------------------

When we do a full GDSC analysis, the cell lines span a set of TCGA tissues
(e.g., COREAD, BLCA) and generally we want to perform the analysis not on all
cellines at the same time but each type of tissues independently.

Besides, you may then wish to have data packages not only for a given TCGA
tissue but also for a given company (if your DrugDecode file is filled properly; see later).

The recommended way is to used the :class:`gdsctools.gdsc.GDSC` class that will
help you in this task.


First, you need to prepare the input data. Create a directory and add these
files::

    - a unique IC50 file
    - The genomic features files for each type of tissues.
    - The DrugDecode file


The genomic feature must be named as follows::

  <prefix>_BLCA.csv
  <prefix>_COREAD.csv
  ...

The name of the TCGA can include **ALL**, **PANCAN** and will be used later to
create the directories for each data paakage.

The important point being that there must be an underscore only and followed by
the TCGA tag.

The GDSC class will then loop over the TCGA cases and create data packages.

::

    from gdsctools import GDSC
    gg = GDSC("IC50.csv", "DrugDecode.csv", "GF_*.csv")
    gg.anaalyse()


This may take hours to finalise: the ANOVA and creation of all images will be
done for each TCGA.

This may be parallelised since each input Genomic Feature analysis is
independent::

    gg_blca = GDSC("IC50.csv", "DrugDecode.csv", "GF_BLCA.csv")
    gg_blca.analyse()

    gg_coread = GDSC("IC50.csv", "DrugDecode.csv", "GF_COREAD.csv")
    gg_coread.analyse()


In an error occurs for one Genomic Feature file, the analysis we jump to the
next file. You may need to check re-run the specific TCGA tissue analysis your
self when an error occured (meaning you do not need to re-run everything).


Once done, you should have all data packages locally in the directory where you
ran the scripts.


The next step is to read back all those results and create data pacakges
dedicated to a company. Based on the DRUG_DECODE file::

    gg = GDSC("IC50.csv", "DrugDecode.csv", "GF_*.csv")
    gg.create_data_packages_for_companies()


For each companies, which names can be checked with::

    gg.companies

a new directory (data package) is created locally


For now, it is important to run this in the same directory where previous
pacakges were created.


Again thiis may be parallelised::

    for each company in gg.companies:
        single = GDSC("IC50.csv", "DrugDecode.csv", "GF_*.csv")
        single.create_data_packages_for_companies([company])


Create summary pages
-----------------------

Following the creating of the "all" TCGA packages and the dedicated packages for
all companies, you end up with quite a few directories. This command will create
summary HTML page to ease your life::

    gg.create_summary_pages()


This must be called after :meth:`analyse` and :meth:`create_data_packages_for_companies`.


































































