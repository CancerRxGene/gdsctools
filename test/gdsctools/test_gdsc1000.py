from gdsctools.gdsc1000 import GDSC1000

try:
    from tempfile import TemporaryDirectory
except:
    from tempfile import mkdtemp
    class TemporaryDirectory(object):
        def __init__(self):
            self.name = mkdtemp()
        def cleanup(self):
            import shutil
            shutil.rmtree(self.name)

def test_download():
    fh = TemporaryDirectory()
    gg = GDSC1000(data_folder_name=fh.name)
    gg.download_data()
    gg.annotate_all()

    gg = GDSC1000(data_folder_name=fh.name)
    gg.load_data()

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
    gg.filter_by_gene("Core Genes")
    assert len(gg.genomic_df)
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

    fh.cleanup()


