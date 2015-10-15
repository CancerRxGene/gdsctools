"""Some Readers



"""
import pandas as pd
import pylab
import easydev

class Reader(object):
    """

    
    """
    def __init__(self, filename, sep="\t"):
        """

        :param input: could be a filename

        """
        self._filename = filename
        self._sep = sep

    def read_matrix_from_r(self, name):
        print("Reading matrix %s " % (name))
        self.session.run("rnames = rownames(%s)" % name)
        self.session.run("cnames = colnames(%s)" % name)
        self.session.run("data = %s" % name)

        cnames = self.session.cnames
        rnames = self.session.rnames
        data = self.session.data
        df = pd.DataFrame(data=data.copy())
        df.columns = [x.strip() for x in cnames]
        df.index = [x.strip() for x in rnames]
        return df

    def __str__(self):
        self.df.info()
        return ""


class CosmicRows(object):
    """Parent class to IC50 and GenomicFeatures"""
    def _get_cosmic(self):
        return list(self.df.index)
    cosmicIds = property(_get_cosmic, doc="return list of cosmic ids")


class IC50(Reader, CosmicRows):
    """Reader of IC50 data set

    Matrix containing the IC50 for each drug. First column must
    contain the COSMIC ID and other columns are the DRUGS

    DRUG label must be "Drug_XX_IC50"

    Extra columns will be ignored.

    """
    def __init__(self, filename='ANOVA_input.txt', sep="\t"):
        super(IC50, self).__init__(filename, sep=sep)

        if isinstance(filename, str):
            self.rawdf = pd.read_csv(self._filename, sep=self._sep)
            columns = ['COSMIC ID']
            columns += [x for x in self.rawdf.columns if x.startswith('Drug')]
            self.df = self.rawdf[columns].copy() # is copy  required ?
            self.df.set_index('COSMIC ID', inplace=True)
        else:
            self.df = filename.df.copy()

    def _get_drugs(self):
        return list(self.df.columns)
    drugIds = property(_get_drugs)


class GenomicFeatures(Reader, CosmicRows):
    def __init__(self, filename=None, sep="\t"):
        # first reset the filename to the shared data (if not provided)
        if filename is None:
            filename = easydev.get_share_file('gdsctools', 'data', 
                            'genomic_features.tsv')
        super(GenomicFeatures, self).__init__(filename)

        if isinstance(filename, str):
            self.df = pd.read_csv(self._filename, sep=self._sep)
            assert 'COSMIC ID' in self.df.columns, \
                "the features input file must contains a column named COSMIC ID"
            self.df.set_index('COSMIC ID', inplace=True)
        else:
            self.df = filename.df.copy()

        # There are several types of features e.g., mutation, CNA,
        # methylation but all are stored within the same file
        # Besides, the 3 first columns are
        #
        # Sample Name
        # Tissue Factor Value
        # MS-instability Factor Value
        #
        # Note the spaces
        # columsn are tab-delimited

    def _get_features(self):
        return list(self.df.columns)
    features = property(_get_features)




class PANCAN(Reader):
    def __init__(self, filename="PANCAN_simple_MOBEM.rdata"):
        super(PANCAN, self).__init__(filename)
        # Remove R dependencies
        from biokit.rtools import RSession
        self.session = RSession()
        self.session.run('load("%s")' %self._filename)
        self.df = self.read_matrix_from_r('MoBEM')


class Extra(Reader):
    def __init__(self, filename="djvIC50v17v002-nowWithRMSE.rdata"):
        super(Extra, self).__init__(filename)
        # Remove R dependencies
        from biokit.rtools import RSession
        self.session = RSession()
        self.session.run('load("%s")' %self._filename)

        # 3 identical matrices containing AUC, IC50 and
        self.dfAUCv17= self.read_matrix_from_r('dfAUCv17')
        self.dfIC50v17 = self.read_matrix_from_r('dfIC50v17')
        # Residual
        self.dfResv17 = self.read_matrix_from_r('dfResv17')

        # This df holds the xmid/scale parameters for each cell line
        # Can be visualised using the tools.Logistic class.
        self.dfCL= self.read_matrix_from_r('dfCL')

        # There is an extra matrix called MoBEM, which is the same as in the
        # file

    def hist_residuals(self, bins=100):
        """Plot residuals across all drugs and cell lines"""
        data = [x for x in self.dfResv17.fillna(0).values.flatten() if x != 0]
        pylab.clf()
        pylab.hist(data, bins=bins, normed=True)
        pylab.grid(True)
        pylab.xlabel('Residuals')
        pylab.ylabel(r'\#')

    def scatter(self):
        from biokit import scatter
        s = scatter.ScatterHist(self.dfCL)
        s.plot(kargs_histx={'color':'red', 'bins':20},
                kargs_scatter={'alpha':0.9, 's':100, 'c':'b'},
                kargs_histy={'color':'red', 'bins':20})

    def hist_ic50(self, bins=100):
        data = [x for x in self.dfIC50v17.fillna(0).values.flatten() if x != 0]
        pylab.clf()
        pylab.hist(data, bins=bins, normed=True)
        pylab.grid(True)
        pylab.xlabel('IC50')
        pylab.ylabel(r'\#')

    def hist_auc(self, bins=100):
        data = [x for x in self.dfAUCv17.fillna(0).values.flatten() if x != 0]
        pylab.clf()
        pylab.hist(data, bins=bins, normed=True)
        pylab.grid(True)
        pylab.xlabel('AUC')
        pylab.ylabel(r'\#')


