from gdsctools import pipelines, ic50_test

from nose.plugins.attrib import attr



class TestPipeline(object):

    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.prog = "gdsctools"
        klass.filename = ic50_test.filename
        klass.params = {'prog': klass.prog, 'filename': klass.filename}

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""

    def setUp(self):
        """This method is run once before _each_ test method is executed"""

    def teardown(self):
        """This method is run once after _each_ test method is executed"""

    def test_version(self):
        pipelines.anova_pipeline([self.prog, '--version'])

    def test_help(self):
        pipelines.anova_pipeline([self.prog, '--help'])
        pipelines.anova_pipeline([self.prog])

    def test_print_drug_names(self):
        pipelines.anova_pipeline([self.prog, '--input-ic50', self.filename,
            '--print-drug-names'])
    def test_print_tissue_names(self):
        pipelines.anova_pipeline([self.prog, '--input-ic50', self.filename,
            '--print-tissue-names'])
    def test_print_feature_names(self):
        pipelines.anova_pipeline([self.prog, '--input-ic50', self.filename,
            '--print-feature-names'])

    def test_odaf(self):
        pipelines.anova_pipeline([self.prog, '--input-ic50', self.filename,
            '--drug', 'Drug_1047_IC50', '--no-html'])

    def test_odaf_fast(self):
        pipelines.anova_pipeline([self.prog, '--input-ic50', self.filename,
            '--drug', 'Drug_1047_IC50', '--no-html', '--fast'])

    def test_odof(self):
        pipelines.anova_pipeline([self.prog, '--input-ic50', self.filename,
            '--drug', 'Drug_1047_IC50', '--no-html', '--feature', 'TP53_mut'])
        pipelines.anova_pipeline([self.prog, '--input-ic50', self.filename,
            '--drug', 'Drug_1047_IC50', '--feature', 'TP53_mut'])
    
    def test_adof(self):
        pipelines.anova_pipeline([self.prog, '--input-ic50', self.filename,
            '--no-html', '--feature', 'TP53_mut'])
   
    @attr('skip')
    def test_adaf(self):
        pipelines.anova_pipeline([self.prog, '--input-ic50', self.filename,
            '--no-html'])


    def test_license(self):
        pipelines.anova_pipeline([self.prog, '--license'])
    def test_summary(self):
        pipelines.anova_pipeline([self.prog, '--summary'])
    def test_summary(self):
        pipelines.anova_pipeline([self.prog, '--testing'])
