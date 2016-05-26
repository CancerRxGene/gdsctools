.. index:: IC50

.. _data:

Data Format and Readers
============================

CSV and TSV formats
--------------------
The formats used in **GDSCTools** are :term:`CSV` based but :term:`TSV` formatted files may be accepted. However, we would encourage to use CSV formats as much as possible. Results will be saved in CSV format.


There are different CSV files used in the analysis. They all have dedicated
readers and expect specific names in the headers as explained hereafter.
So far, we have  3 types of input files (defined in the
:mod:`~gdsctools.readers` module):

- :class:`~gdsctools.readers.IC50`
- :class:`~gdsctools.readers.GenomicFeatures`
- :class:`~gdsctools.readers.DrugDecode`

To be identified CSV/TSV files must have the proper extension that is **.csv** or **.tsv**; they may be compressed e.g., in gzip format but keep the **.csv** or **.tsv** extension in the name since it is interpreted to recognised the format of the file (e.g, ic50.csv.gz)


.. note:: with CSV files, if a column's name is ambiguous in the header (i.e, contains already a comma), then it should be enclosed within double quotes.


IC50
------

The most important input data is what we call the IC50 file. It is
simply a CSV (or TSV) file with a header and a set of rows. The header's column name must be labelled after each drug considered plus an additional column named **COSMIC_ID**. The order of the columns is not important. Each row contains a unique COSMIC identifier and a set of IC50s. Note that we speak of IC50s but one can populate the file with anything (e.g., AUCs).

Although the name in the header does not matter note that if one column's name
starts with **Drug** then all column that do not start with **Drug** are ignored (except the special one that contains the COSMIC IDs). This feature was implemented to account for old data files that stores all Drug identifiers and also all genomic features within a single file.

Here is an example of IC50 input::

    COSMIC_ID, Drug_1_IC50, Drug_20_IC50, Other
    111111,    0.5,         0.8,          10
    222222,    1,           2,            10

If you save that example in a file, you can read it with the
:class:`~gdsctools.readers.IC50` class as follows:

.. doctest::

    >>> from gdsctools import IC50
    >>> r = IC50('source/ic50_tiny.csv')
    >>> r.drugIds
    [1, 20]


.. note:: the columns' names should be identifiers (not drug names). There
    are two main reasons. The first one is that it allows us to keep anonymous
    all drug names and targets. The second reason is that many characteristics
    such as plate number and drug concentration may be associated with a drug
    identifier. This should be stored in a different table rather than in
    the name. It can then be handled and interpreted using the DrugDecode
    file (see below).

.. note:: column without a name are ignored.


.. seealso:: developers should look at the references for more
    functionalities of the :class:`~gdsctools.readers.IC50`
    class (e.g., filter by tissues, removing drugs, visualisation of IC50s).



Genomic Features
---------------------

The **ANOVA** analysis computes the associations between the Drug IC50s and
genomic features. The mapping between these two data sets is performed on a common column named **COSMIC_ID**, which should contain the same COSMIC identifiers. If not, only the intersection will be kept in sub-sequent analysis.

In addition to the COSMIC identifiers, the following columns may be provided but are not strictly speaking required::

    - TISSUE_FACTOR
    - MSI_FACTOR
    - MEDIA_FACTOR

If not provided, the tissue, :term:`MSI` and :term:`MEDIA` factors will not be taken into account in the regression analysis. If the :term:`TCGA` tissue is not provided, it is created and set to *unidentified*.

.. note::
    .. versionchanged:: 0.9.11
        A column called 'Sample Name' was interpreted if found. This is not
        the case anymore. It is actually removed now.


All remaining columns are assumed to be genomic features.

.. warning:: In the current version, all columns starting 
    with `Drug_` are removed without warning.


Here is a simple example::

    COSMIC_ID, TISSUE_FACTOR, MSI_FACTOR, BRAF_mut, gain_cna
    111111, lung_NSCLC,  1, 1, 0
    222222, prostate,    1, 0, 1

It can be saved and read as follows:

.. doctest::

    >>> from gdsctools import GenomicFeatures
    >>> gf = GenomicFeatures('source/gf_tiny.csv')
    >>> gf
    GenomicFeatures <Nc=2, Nf=2, Nt=2>

In **GDSCTools**, we provide a :download:`zipped Genomic Features file<../../gdsctools/data/genomic_features.tsv.gz>`. It contains about 1000 cell lines and 47 genomic features (gene mutations). A more complex file tagged v17 it also provided with about 600 features :download:`v17 genomic feature <../../gdsctools/data/genomic_features_v17.csv.gz>`.

By default, the creation of an ANOVA class we read that file automatically. Of
course, you may provide your own. The :class:`~gdsctools.readers.GenomicFeatures` if created without input contains the default file mentionned here above::


    >>> from gdsctools import GenomicFeatures
    >>> gf = GenomicFeatures()
    >>> print(gf)
    Genomic features distribution
    Number of unique tissues 27
    Number of unique features 47 with
    - Mutation: 47
    - CNA (gain): 0
    - CNA (loss): 0

Drug Decoder
----------------

The :class:`~gdsctools.readers.DrugDecode` class reads a CSV file that contains information about a drug and its target(s). It must contain 3 columns named as
follows::

    DRUG_ID,        DRUG_NAME,  DRUG_TARGET
    Drug_999_IC50,  Erlotinib,  EGFR
    Drug_1039_IC50, SL 0101-1,  "RSK, AURKB, PIM3"


This column will be used if provided::

    - WEBRELEASE
    - OWNED_BY

In addition, these column may be populated for later use::

    - CHEMSPIDER_ID
    - PUBCHEM_ID
    - CHEMBL_ID

An example can be read as follows:

.. doctest::

    >>> from gdsctools import DrugDecode, datasets
    >>> drug_filename = datasets.testing.drug_test_csv.location
    >>> dd = DrugDecode(drug_filename)
    >>> dd.get_name(1047)
    'Nutlin-3a'
    >>> dd.df.ix[999]
    CHEMBL_ID              NaN
    CHEMSPIDER_ID          NaN
    DRUG_NAME        Erlotinib
    DRUG_TARGET           EGFR
    OWNED_BY               NaN
    PUBCHEM_ID             NaN
    WEBRELEASE             NaN
    Name: 999, dtype: object



DrugDecode files are not required for the analysis but are used by
:class:`gdsctools.anova_report.ANOVAReport` to fill the HTML reports.


You can also run the analysis and set the drug names and target later on as
follows using the :class:`~gdsctools.readers.drug_annotations` method::

    from gdsctools import *
    an = ANOVA(ic50_test)
    an.anova_all()
    results = an.anova_all()
    dd = DrugDecode("v19_drug_decode.csv")
    newdf = dd.drug_annotations(results.df)











