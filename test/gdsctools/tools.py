from gdsctools.readers import IC50
import easydev


def get_data(filename=None):
    if filename is None:
        filename = easydev.get_share_file('gdsctools', 'data', 
                'IC50_10drugs.tsv')
    import os
    r = IC50(filename)


