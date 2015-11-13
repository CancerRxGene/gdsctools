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


IC50
------

The IC50s file should be a CSV file where each column correpond to a
drug, each row to cell line. The header is compulsary and indicate the name of
the drugs. An additional column with the cosmic identifier must be provided.
This column must be named **COSMIC ID** (note the space).

The name of the drugs' columns is not interpreted but used for representing a
drug in the different plots' title or labels.

The order of the columns' names is not important.

Extra columns (e.g., tissue, sample name, MSI, features) will be ignored.

Here is a simple example of a valid CSV file::

    COSMIC ID, Drug_1_IC50, Drug_20_IC50
    111111,    0.5,         0.8
    222222,    1,           2

.. seealso:: developers should look at the references for more 
    functionalities of the :class:`~gdsctools.readers.IC50`  
    class (e.g., filter by tissues, removing drugs, visualisation of IC50s).



Genomic Features
---------------------


Drug Decoder
----------------
