"""
Use Omnibem to filter annotation data 
=========================================

This example shows how to create a new genomic feature data set
from an annotation file.
"""

#####################################################
#
from gdsctools import *
data = gdsctools_data("test_omnibem_genomic_alterations.csv.gz")
bem = OmniBEMBuilder(data)
bem.filter_by_gene_list(gdsctools_data("test_omnibem_genes.txt"))
bem.plot_number_alteration_by_tissue()


#####################################################
# Finally, create a MoBEM dataframe
mobem = bem.get_mobem()

# features
bem.filter_by_type_list(["Methylation"])
mobem = bem.get_mobem()

# Then, let us create a dataframe that is compatible with
# GenomicFeature. We just need to make sure the columns are correct
mobem[[x for x in mobem.columns if x!="SAMPLE"]]
gf = GenomicFeatures(mobem[[x for x in mobem.columns if x!="SAMPLE"]])

######################################""
# The final volcano plot
an = ANOVA(ic50_test, gf, verbose=False)
results = an.anova_all(animate=False)
results.volcano()

