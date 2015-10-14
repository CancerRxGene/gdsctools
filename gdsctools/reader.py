"""Some Readers



"""
import pandas as pd
import pylab


class Reader(object):
    def __init__(self, filename):
        self.filename = filename

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

    def save(self, tag='test'):
        # There are 14 digits in the original file, so let us write
        # with 14 digits as well.
        self.ic50.to_csv("test_ic50.tsv", sep="\t", float_format="%.14f")
        self.features.to_csv("test_features.tsv", sep="\t",
                float_format="%.14f")


class IC50(Reader):
    """Reader of ANOVA data set created by Francesco

    This data seems to have all genomic feature and drug IC50 in the same matrix
    and each row correspond to a cellline.


    """
    def __init__(self, filename='ANOVA_input.txt'):
        super(IC50, self).__init__(filename)

        self._read_data()

    def _read_data(self):
        self.rawdf = pd.read_csv(self.filename, sep='\t')

        columns = ['COSMIC ID'] +[x for x in self.rawdf.columns if x.startswith('Drug')]
        self.ic50 = self.rawdf[columns].copy() # is copy  required ?
        self.ic50.set_index('COSMIC ID', inplace=True)

        columns = [x for x in self.rawdf.columns if
                x.startswith('Drug') is False]
        self.features = self.rawdf[columns]
        self.features.set_index('COSMIC ID', inplace=True)


class Feature(Reader):
    def __init__(self, filename):
        pass


class PANCAN(Reader):
    def __init__(self, filename="PANCAN_simple_MOBEM.rdata"):
        super(PANCAN, self).__init__(filename)
        # Remove R dependencies
        from biokit.rtools import RSession
        self.session = RSession()
        self.session.run('load("%s")' %self.filename)
        self.df = self.read_matrix_from_r('MoBEM')


class Extra(Reader):
    def __init__(self, filename="djvIC50v17v002-nowWithRMSE.rdata"):
        super(Extra, self).__init__(filename)
        # Remove R dependencies
        from biokit.rtools import RSession
        self.session = RSession()
        self.session.run('load("%s")' %self.filename)

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
        pylab.ylabel('#')

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
        pylab.ylabel('#')

    def hist_auc(self, bins=100):
        data = [x for x in self.dfAUCv17.fillna(0).values.flatten() if x != 0]
        pylab.clf()
        pylab.hist(data, bins=bins, normed=True)
        pylab.grid(True)
        pylab.xlabel('AUC')
        pylab.ylabel('#')


