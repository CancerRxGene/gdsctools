


def test():
    pass



    # try some summary / info
    prog = "gdsctools"
    filename = "ANOVA_input.txt"
    params = {'prog':prog, 'filename':filename}

    cmd = "%(prog)s --input-ic50 %(filename)s --print-drug-names" % params
    cmd = "%(prog)s --input-ic50 %(filename)s --print-feature-names" % params
    cmd = "%(prog)s --input-ic50 %(filename)s --print-tissue-names" % params
    cmd = "%(prog)s --help" 
    cmd = "%(prog)s --version" 



"""
  -I INPUT_IC50, --input-ic50 INPUT_IC50
                        A file in TSV format with IC50s. First column should
                        be the COSMIC identifiers Following columns contain
                        the IC50s for a set of drugs. The header must be
                        COSMIC ID, Drug_1_IC50, Drug_2_IC50, ...
  -F INPUT_FEATURES, --input-features INPUT_FEATURES
                        A matrix of genomic features. First column is made of
                        COSMIC identifiers that should match those from the
                        IC50s matrix. Then 3 compulsary columns are required
                        named 'Sample', 'Tissue', 'MSI'. Other columns should
                        be the genomic features. There are recognised if the
                        (1) ends in _mut for mutation, or starts with loss or
                        gain for the CNA cases.
  --output-directory DIRECTORY
                        directory where to save images and HTML files.
  -d DRUG, --drug DRUG  The name of a valid drug identifier to be found in the
                        header of the IC50 matrix
  -f FEATURE, --feature FEATURE
                        The name of a valid feature to be found in the Genomic
                        Feature matrix
  -t TISSUE, --tissue TISSUE
                        The name of a specific cancer type i.e., tissue to
                        restrict the analysis to
  --include-drugs-in DRUGS [DRUGS ...]
                        todo
--include-features-in FEATURES [FEATURES ...]
                        There are many genomic features included by default.
                        You may want to select a subset of them. No effect if
                        --feature is provided
  --exclude-msi         Include msi factor in the analysis
  --summary             Print summary about the data (e.g., tissue)
  --test                Use a small IC50 data set and run the one-drug-one-
                        feature analyse with a couple of unit tests.
  --license             Print the current license
  --fast                If provided, the code will use more memory and should
                        be 10-30% faster. (1.2G for 265 drugs and 680
                        features)
"""
