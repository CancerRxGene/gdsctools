from gdsctools.drugs import ChemSpiderSearch

# Requires to use bioservices and a token for chemspider.
def _test_chemspider():

    cc = ChemSpiderSearch(['ZM447439'])
    cc.search_in_chemspider()
    cc.search_from_smile_inchembl()
    df = cc.find_chembl_ids()
    assert df.ix[0].CHEMBL_ID == 'CHEMBL202721'

