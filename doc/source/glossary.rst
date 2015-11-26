Glossary
===========


.. https://tcga-data.nci.nih.gov/tcga/tcgaDataType.jsp


.. index:: ic50, ec50

.. glossary::
    :sorted:

    IC50
        The IC50 is the half maximal Inhibitory Concentration that is a
        measure of effectiveness of a drug or substance in inhibiting
        a specific biological or biochemical function.

        Sometimes, IC50 values are converted to the pIC50 scale that is
        :math:`-\log10 (IC50)`.

        pIC50 is usually given in terms of molar concentration (mol/L or M).
        Therefore to obtain pIC50, an IC50 should be specified in units
        of M. When IC50 is expressed in microM or nanoM, it will need to
        be converted to M before conversion to pIC50.

    EC50
        The EC50 is the half maximal Effective Concentration and refers to the
        concentration of a drug, antibody or toxicant which induces a
        response halfway between the baseline and maximum after a specified
        amount of time. It is used as a measure of drug's potency

    MSI
        Microsatellite Instability (MSI) are markers indicating the 
        presence or absence of a MSI shift, allele
        homozygosity/heterozygosity, and loss of heterozygosity (LOH) observed
        in the tumor sample for each participant

    FDR
        In genome research, it is common to examine a large number of features.
        When testing multiple hypotheses, one must guard against an 
        abundance of false positive results. A criterion for error control is 
        the false discovery rate (FDR), which is the expected proportion of 
        falsely rejected hypotheses. This error rate is equal to FWER when 
        all null hypotheses are true but is smaller otherwise. 
        Benjamini and Hochberg, proposed a step-down procedure to control 
        FDR for independent test statistics. This method is currently
        recommended in **GDSCTools** and is the default method gor multiple
        testing correction. See e.g.,
        `article <http://bioinformatics.oxfordjournals.org/content/21/6/781.full>`_


    TCGA
        The Cancer Genome Atlas (TCGA) is a project to catalogue
        genetic mutations responsible for cancer, using genome sequencing and
        bioinformatics. See  `TCGA homepage <http://cancergenome.nih.gov/>`_

    COSMIC
        COSMIC is a global resource for
        information on somatic mutations in human cancer, combining curation of
        the scientific literature with tumor resequencing data from the Cancer
        Genome Project at the Sanger Institute, U.K. See 
        the `COSMIC website <http://cancer.sanger.ac.uk/cosmic>`_ or
        `pubmed reference  <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2705836/>`_ for more details

    ODOF
        acronym for One Drug One Feature mode
    
    ODAF
        acronym for One Drug All Feature mode
    
    ADAF
        acronym for All Drug All Feature mode
   
    CSV
        acronym for Comma Separated Values, a file format. Extension must
        be `.csv`.

    TSV
        acronym for Tabular Separated Values, a file format. Extension must
        be `.tsv`.

    PANCAN
        alias for set of cancer tissues (unlike cancer-specific tissue).

    Terminal
        Under Unix-like operating systems, a terminal is a program that 
        run a shell. A unix terminal (Linux 
        or Mac) starts with a specialised shell (e.g. bash shell). Under 
        Windows, the "command prompt" available in All Programs -> Accessories
        is an entry point to a terminal for typing computer commands
       
    CLI
        A Command Line Interface (CLI) is an interface where the user
        types a command (text) and presses the return key to execute that
        command. 

    shell
        A shell is a program that provides the traditional, text-only user 
        interface for Linux and other Unix-like operating systems. 
        It is a specialised :term:`CLI` that is a command-line shell (e.g., 
        bash) where users can execute programs.

    OLS
        An ordinary least squares (OLS) or linear least squares is a
        method for estimating the unknown parameters in a linear regression
        model, with the goal of minimizing the differences between the observed
        responses in some arbitrary dataset and the responses predicted by the
        linear approximation of the data
