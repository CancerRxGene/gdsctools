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
- :class:`~gdsctools.readers.DrugDecoder`
- **Concentrations** (not yet used).

To be identified CSV/TSV files must have the proper extension that is **.csv** or **.tsv**; they may be compressed e.g., in gzip format but keep the **.csv** or **.tsv** extension in the name since it is interpreted to recognised the format of the file (e.g, ic50.csv.gz)


.. note:: with CSV files, if a column's name is ambiguous in the header (i.e, contains already a comma), then it should be enclosed within double quotes.


IC50
------

The most important input is the file that contains the IC50s. That input file
contains a set of row where each row must have 1 unique COSMIC identifier and a
set of IC50s; one IC50 per drug. So, the compulsary header must indicate the name of the drugs and the column that contains the cosmic identifiers, which must be named **COSMIC ID** (note the space). All columns that starts with **Drug** are interpreted as the IC50s for each drug considered. Othe columns are ignored. The order of the columns is not important. Here is a valid example::

    COSMIC ID, Drug_1_IC50, Drug_20_IC50
    111111,    0.5,         0.8
    222222,    1,           2

If you save that example in a file, you can read it with the
:class:`~gdsctools.readers.IC50` class as follows:

.. doctest::

    >>> from gdsctools import IC50
    >>> r = IC50('source/ic50_tiny.csv')
    >>> r.drugIds
    ['Drug_1_IC50', 'Drug_20_IC50']

    
.. seealso:: developers should look at the references for more 
    functionalities of the :class:`~gdsctools.readers.IC50`  
    class (e.g., filter by tissues, removing drugs, visualisation of IC50s).



Genomic Features
---------------------

The **ANOVA** analysis computes the associations between the Drug IC50s and
genomic features. The file containing the Genomic Features must map to the IC50s file that is it must contains a column named **COSMIC ID** with the same COSMIC identifiers. Besides, 2 compulsary columns are required. One that contains the tissue names and one with information about the :term:`MSI` factor. Those 2 columns must be named ::

    - 'Tissue Factor Value'
    - 'MS-instability Factor Value'

We may provide alternative (simple) names in the futures. You may also have an additional informative column named:: 

    - 'Sample Name'

Finally, remaining columns are assumed to be related to genomic features. 
Note that columns starting with `Drug_` are removed without warning for now. 
Here is a simple example::
    
    COSMIC ID, Tissue Factor Value, Sample Name, MS-instability Factor Value, BRAF_mut, gain_cna
    111111, lung_NSCLC, 201T,  1, 1, 0
    222222, prostate,   22RV1, 1, 0, 1

It can be saved and read as follows:

.. doctest::

    >>> from gdsctools import GenomicFeatures
    >>> gf = GenomicFeatures('source/gf_tiny.csv')
    >>> gf
    GenomicFeatures <Nc=2, Nf=2, Nt=2>

In **GDSCTools**, we provide a :download:`zipped Genomic Features file<../../share/data/genomic_features.tsv.gz>`. It contains about 1000 cell lines and 700 genomic features. 

By default, the creation of an ANOVA class we read that file automatically. Of
course, you may provide your own. The :class:`~gdsctools.readers.GenomicFeatures` if created without input contains the default file mentionned here above::


    >>> from gdsctools import GenomicFeatures
    >>> gf = GenomicFeatures()
    >>> print(gf)
    Genomic features distribution
    Number of unique tissues 27
    Number of unique features 677 with
    - Mutation: 270
    - CNA (gain): 116
    - CNA (loss): 291

Drug Decoder
----------------

The :class:`~gdsctools.readers.DrugDecoder` class reads a CSV file that contains information about a drug and its target(s). It must contain 3 columns named as
follows::

    DRUG_ID,        DRUG_NAME,  DRUG_TARGET
    Drug_999_IC50,  Erlotinib,  EGFR
    Drug_1039_IC50, SL 0101-1,  "RSK, AURKB, PIM3"

An example can be read as follows:

.. doctest::

    >>> from gdsctools import DrugDecoder, datasets
    >>> drug_filename = datasets.testing.drug_test_csv.location
    >>> dd = DrugDecoder(drug_filename)
    >>> dd.get_name('Drug_999_IC50')
    'Erlotinib'

















