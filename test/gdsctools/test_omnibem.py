from gdsctools.omnibem import OmniBEMBuilder
from gdsctools import gdsctools_data

omnibem_data = gdsctools_data("test_omnibem_genomic_alteration.csv.gz")
omnibem_genes = gdsctools_data("test_omnibem_genes.txt")


def test_omnibem():
    ob = OmniBEMBuilder(omnibem_data)
    ob.filter_by_gene_list(omnibem_genes)
    mobem = ob.get_mobem()
    assert mobem[mobem.columns[3:]].sum().sum() == 4971

    #
    ob.plot_number_alteration_by_tissue()
    ob.plot_alterations_per_cellline()
    ob.get_significant_genes()


    ob = OmniBEMBuilder(omnibem_data)
    ob.filter_by_gene_list(omnibem_genes)
    ob.filter_by_type_list(["Methylation"])
    mobem = ob.get_mobem()
    assert mobem[mobem.columns[3:]].sum().sum() == 127

    ob = OmniBEMBuilder(omnibem_data)
    ob.filter_by_gene_list(omnibem_genes)
    ob.filter_by_tissue_list(["HNSC"])
    mobem = ob.get_mobem()
    assert mobem[mobem.columns[3:]].sum().sum() == 115

    ob = OmniBEMBuilder(omnibem_data)
    ob.filter_by_gene_list(omnibem_genes)
    ob.filter_by_sample_list(["SNU-423"])
    #mobem = ob.get_mobem()
    #assert mobem[mobem.columns[3:]].sum().sum() == 4971
