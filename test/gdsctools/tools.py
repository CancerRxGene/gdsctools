from gdsctools.readers import IC50



def get_data(filename='IC50_10drugs.tsv'):
    import os
    try:
        r = IC50(filename)
    except:
        r = IC50('test/gdsctools' + os.sep +  filename)
    return r


