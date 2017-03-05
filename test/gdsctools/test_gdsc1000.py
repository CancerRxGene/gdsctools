from gdsctools import GDSC1000


def test_download(tmpdir):
    p = tmpdir.mkdir("download")
    name = "download/"

    gg = GDSC1000(data_folder_name=name)
    gg.download_data()

    gg = GDSC1000(data_folder_name=name)
    gg.load_data(annotation=False)

    # now some filtering. Let us start with alteration_type
    assert gg.genomic_df.shape == (42231, 8)
    gg.filter_by_alteration_type(["DELETION"])
    assert gg.genomic_df.shape == (1614, 8)
    gg.reset_genomic_data()
    assert gg.genomic_df.shape == (42231, 8)

    # by gene
    gg.filter_by_gene(["A2ML1"])
    assert len(gg.genomic_df)
    gg.reset_genomic_data()

    # by core genes
    # numbers labelled (for sure) were found in Liz document

    gg = GDSC1000(data_folder_name=name)
    gg.load_data() # Here, we include the annotations

    gg.filter_by_gene("Core Genes")
    assert len(gg.genomic_df['GENE'].unique()) == 310  # For sure
    assert len(gg.genomic_df)
    assert gg.get_genomic_info()["cosmic"].METHYLATION == 108 # for sure
    assert gg.get_methylation_info().ix[0,0] == 338 # for sure
    gg.get_cna_info()
    gg.reset_genomic_data()

    #  by tissues
    gg.filter_by_tissue_type(["COAD/READ"])
    assert len(gg.genomic_df)
    gg.reset_genomic_data()

    # by cell line
    gg.filter_by_cell_line()
    gg.filter_by_cell_line(["SW1116"])
    assert len(gg.genomic_df)
    gg.reset_genomic_data()

    # by cosmic id
    gg.filter_by_cosmic_id()
    gg.filter_by_cosmic_id([909746])
    assert len(gg.genomic_df)
    gg.reset_genomic_data()

    gg.filter_by_recurrence()
    gg.filter_by_recurrence(3)

    gg.make_matrix()

