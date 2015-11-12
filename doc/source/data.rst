.. index:: IC50

.. _data:

Data format and readers
============================

You will need several type of files to perform the analysis. 

- **IC50**: 
- **Genomic Features**
- **DRUG_DECODER**
- **Concentrations**

The first one is compulsary, the second and third ones are optional. 
The last one is not used yet.

The may format used is a comma-separated values (CSV) or tabular-separeted values (TSV). It may not be optimal but it is human readable and the file-size being small enough, this is manageable. 

The choise between CSV and TSV has not been made yet. So, we accept both format. However, the output format will always be comma separated. If a header name is ambiguous (contains already a comma), then it be be enclosed within double quotes::

Note that all data format described here can be read with one of the
**GDSCTools** function or class from the :mod:`gdsctools.readers` module as we
will show here below.


.. autosummary::

    gdsctools.readers
    gdsctools.IC50
    gdsctools.DrugDecoder
    gdsctools.GenomicFeatures



IC50
------


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
