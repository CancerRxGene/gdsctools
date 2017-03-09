from gdsctools.scripts import anova as pipelines
from gdsctools import  ic50_test


prog = "gdsctools"
filename = ic50_test.filename

class TestPipeline(object):

    def test_version(self):
        pipelines.main([prog, '--version'])

    def test_help(self):
        pipelines.main([prog, '--help'])
        pipelines.main([prog])

    def test_print_drug_names(self):
        pipelines.main([prog, '--input-ic50', filename,
            '--print-drug-names'])
    def test_print_tissue_names(self):
        pipelines.main([prog, '--input-ic50', filename,
            '--print-tissue-names'])
    def test_print_feature_names(self):
        pipelines.main([prog, '--input-ic50', filename,
            '--print-feature-names'])

    def test_odaf(self):
        pipelines.main([prog, '--input-ic50', filename,
            '--drug', '1047', '--no-html'])

    def test_odof(self):
        pipelines.main([prog, '--input-ic50', filename,
            '--drug', '1047', '--no-html', '--feature', 'TP53_mut'])
        pipelines.main([prog, '--input-ic50', filename,
            '--drug', '1047', '--feature', 'TP53_mut'])

    # slow takes about 30-60 seconds
    def test_adaf(self):
        pipelines.main([prog, '--input-ic50', filename,
            '--no-html'])

    def test_summary(self):
        pipelines.main([prog, '--summary'])

    def test_summary(self):
        pipelines.main([prog, '--testing'])
