from gdsctools.omnibem import OmniBEMBuilder
from gdsctools import gdsctools_data

omnibem_data = gdsctools_data("test_omnibem_genomic_alterations.csv.gz")
omnibem_genes = gdsctools_data("test_omnibem_genes.txt")


def test_omnibem():
    ob = OmniBEMBuilder(omnibem_data)
    assert len(ob) ==  56943
    ob.filter_by_gene_list(omnibem_genes)
    mobem = ob.get_mobem()
    assert mobem[mobem.columns[3:]].sum().sum() == 54061

    #
    ob.plot_number_alteration_by_tissue()
    ob.plot_alterations_per_cellline()
    ob.get_significant_genes()


    # filter by cosmic id
    ob = OmniBEMBuilder(omnibem_data)
    ob.filter_by_cosmic_list([910916])
    mobem = ob.get_mobem()
    assert mobem.shape == (1,105)
    assert mobem.ix[0,3:].sum() == 102

    ob = OmniBEMBuilder(omnibem_data)
    ob.filter_by_type_list(["Methylation"])
    mobem = ob.get_mobem()
    assert mobem[mobem.columns[3:]].sum().sum() == 12964

    ob = OmniBEMBuilder(omnibem_data)
    ob.filter_by_tissue_list(["HNSC"])
    mobem = ob.get_mobem()
    assert mobem[mobem.columns[3:]].sum().sum() == 4102

    ob = OmniBEMBuilder(omnibem_data)
    ob.filter_by_sample_list(["SNU-423"])
    mobem = ob.get_mobem()
    assert mobem[mobem.columns[3:]].sum().sum() == 63


    gf = ob.get_genomic_features()
