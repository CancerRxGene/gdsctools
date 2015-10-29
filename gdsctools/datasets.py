import easydev

__all__ = ['dataset', 'ic50_test', 'genomic_features']


class Data(object):
    """A class to hold information about a dataset
    
    
    Can be used as input to :class:`ANOVA` instance.
    """
    def __init__(self):
        self.filename = None
        self.description = "No description"
        self.authors = 'GDSC consortium'

    def __str__(self):
        txt = 'location: %s\n' % self.filename
        txt += 'description: %s\n' % self.description
        txt += 'authors: %s\n' % self.authors
        return txt


def dataset(dataname):
    """Retrieve information about a dataset including location
    
    :param str: a data set's name (e.g., ic50_test)
    :return: a :class:`Data` holder

    ::

        from gdsctools.datasets import dataset
        filename = dataset(ic50_test).filename
        print dataset(i50_test).description
        # or simply
        from gdsctools.datasets import ic50_test
    
    """

    valid = ['ic50_test', 'genomic_features']
    easydev.check_param_in_list(dataname, valid)

    if dataname == 'ic50_test':
        d = Data()
        d.filename = easydev.get_share_file('gdsctools', 
                'data', 'IC50_10drugs.tsv')
        d.description = 'IC50s for 10 public drugs across cell lines'
    elif dataname == 'genomic_features':
        d = Data()
        d.filename = easydev.get_share_file('gdsctools', 
                'data', 'genomic_features.tsv')
        d.descritption = 'Set of genomic features + tissue + sample name + msi'
    return d

#: dataset with IC50s for 10 drugs
ic50_test = dataset('ic50_test')

#: dataset with genomic features for 1001 cell lines and 680 features
genomic_features = dataset('genomic_features')



