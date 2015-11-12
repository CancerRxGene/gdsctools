.. index:: IC50

.. _data:

Data Format and Readers
============================

You will need several type of files to perform the analysis. 

- :class:`~gdsctools.readers.IC50` 
- :class:`~gdsctools.readers.GenomicFeatures`
- :class:`~gdsctools.readers.DrugDecoder`
- **Concentrations**

The first one is compulsary, the second and third ones are optional. 
The last one is not used yet.

The format used for all input and output files is the comma-separated values (CSV). Tabular-separeted values (TSV) may be understood as an input file. The format choice may not be optimal but has the advantage of being human readable. Note that you may compress the file using gzip.

.. note:: the interpreter with look for the .csv or .tsv extension to 
    interpret the file format so you must keep the extension in the 
    compressed files (e.g., data.csv.gz or data.tsv.gz)


.. note:: if the header's names are ambiguous (contains already a comma), then it should be enclosed enclosed within double quotes.

.. autosummary::

    gdsctools.readers
    gdsctools.readers.IC50
    gdsctools.readers.DrugDecoder
    gdsctools.readers.GenomicFeatures



IC50
------

.. autosummary:: gdsctools.readers.IC50

The input matrix may be a tab-separated value file (TSV) but we would encourage
to use comma separated value. The matrix must have at least 2 columns and 2 rows. The first row is the header describing the columns’ contents. One column must be named "COSMIC ID" or "COSMIC_ID". Other columns must be named after the drug. There is no restriction about the naming of the drug but it should be an identifier (string or integer). The order of the column does not matter.

The column "COSMIC ID" contains the cosmic identifiers (cell line). The other
columns should be filled with the IC50s corresponding to a pair of COSMIC ID and
a Drug.

Extra columns (e.g., tissue, sample name, MSI, features) will be ignored.

Here is a simple example of a valid TSV file::

    COSMIC ID   Drug_1_IC50 Drug_20_IC50
    111111      0.5         0.8
    222222      1           2





Genomic Features
---------------------


Drug Decoder
----------------
