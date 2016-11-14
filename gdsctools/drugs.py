"""

Small functionalities to retrive chembl/chemspider identifiers
based on a drug name

"""
from easydev import Progress
from gdsctools import DrugDecode
import pandas as pd


__all__ = ["ChemSpiderSearch"]


class ChemSpiderSearch(object):
    """This class uses ChemSpider and ChEMBL to identify drug name

    .. warning:: this is a draft version in dev mode

    ::

        c = ChemSpiderSearch()
        c.search_in_chemspider()
        c.search_from_smile_inchembl()
        df = c.find_chembl_ids()

    It happens that most of public names can be found
    and almost none of non-public are found. As expected...

    If chemspider, chembl and pubchem are empty, search for the drug name in
    chemspider.

        CHEMSPIDER search:
            if no identifier found, the search if DROPPED
            if 1 identifier found, we keep going using the SMILE identifier
            If more than 1 identifier found, this is AMBIGUOUS.


    If chembl and pubchem, check with unichem
    If chembl, check smiles
    If chembl and chemspider, check smiles ?

    SMILES are not unique

    """
    def __init__(self, drug_decode):
        print("ChemSpiderSearch is still in progress, please do not use")
        self.dd = DrugDecode(drug_decode)
        self.dd_filled = DrugDecode(drug_decode)

        from bioservices.chemspider import ChemSpider
        from bioservices import ChEMBL
        from bioservices import UniChem

        try:
            print('Loading PubChem')
            from bioservices.pubchem import PubChem
            self.puchem = PubChem()
        except:
            # Pubchem was introduced only in dec 2015
            pass

        print('Loading ChEMBL service')
        self.chembl = ChEMBL(cache=True)

        print('Loading ChemSpider service')
        self.chemspider = ChemSpider(cache=True)

        print('Loading UniChem service')
        # in unichem db number is 22 and chembl is 1
        self.unichem = UniChem()

        print('Settings some data aliases')
        self._cs_find = self.chemspider.find
        self._cs_get = self.chemspider.GetExtendedCompoundInfo

        self.drug_ids = sorted(list(self.dd.df.index.values))
        self.drug_names = sorted(list(self.dd.df.DRUG_NAME.values))

    def filling_chembl_pubchem_using_unichem(self):
        """

        """
        N = len(self.drug_ids)
        pb = Progress(N)
        for i,this in enumerate(self.drug_ids):
            entry = self.dd.df.ix[this]
            # if no information is provided, we will need to get it 
            # from chemspider

            # From the database, when chembl is provided, it is unique
            # same for chemspider and pubchem and CAS
            select = entry[['CHEMSPIDER', 'CHEMBL', 'PUBCHEM']]
            if select.count() == 0:
                name = self.dd.df.ix[this].DRUG_NAME
                results = self._cs_find(name)
                if len(results) == 0:
                    # nothing found
                    pass
                elif len(results) == 1:
                    self.dd_filled.df.ix[this].loc['CHEMSPIDER'] = results[0]
                else:
                    # non unique
                    #chemspider = ",".join([str(x) for x in results])
                    self.dd_filled.df.ix[this].loc['CHEMSPIDER'] = results
            pb.animate(i+1)

        # Search in chemspider systematically
        for i, this in enumerate(self.drug_ids):
            entry = self.dd.df.ix[this]
            if select.count() == 1:
                res = self._cs_find(drug)

            pb.animate(i+1)

    def find_chembl_ids(self):
        """


        """
        # don't know how to search for a chembl id given the drug name...
        # so we use chemspider
        #self.search_in_chemspider()

        # but chemspider returns molecular information (not chembl id)
        # so given the smile string, we look back in chembl for valid entries
        #self.search_from_smile_inchembl()

        # finally, get the chembl identifiers
        drugs = []
        chembl_ids = []
        chemspider_ids = []
        smiles_c = []
        smiles_cs = []

        for drug in self.drug_ids:
            try:
                entry = self.results_chembl[drug]

                ids = ",".join([x['chemblId'] for x in entry])
                drugs.append(drug)
                chembl_ids.append(ids)
                ids = ",".join([str(x) for x in self.results[drug]])
            except:
                print('skipping' + drug)
                ids = ",".join([drug, '', '', '', '', ''])
            chemspider_ids.append(ids)

        for drug in self.drug_ids:
            try:
                smiles_c.append(",".join([x['smiles'] for x in
                    self.results_chembl[drug]]))
            except:
                smiles_c.append('')
            try:
                smiles_cs.append(self.results_chemspider[drug]['smiles'])
            except:
                smiles_cs.append('')

        df = pd.DataFrame([drugs, chembl_ids, chemspider_ids, smiles_c,
            smiles_cs],
                index=['DRUG_NAME','CHEMBL_ID','CHEMSPIDER_ID', 'SMILE_CHEMBL',
                    'SMILE_CHEMSPIDER'])
        df = df.T
        return df

    def get_chemspider_ids(self, drug_name):
        res = self._cs_find(drug)
        return res

    def search_in_chemspider(self):
        # Fill results attribute as a dictionary. Keys being the drug id
        # and values are list of chemspider identifiers
        #
        # SB52334 --> SB-52334
        N = len(self.dd)

        pb = Progress(N)
        self.results = {}
        results = []
        for i, index in enumerate(self.dd.df.index):
            drug = self.dd.df.index[i]
            drug_name = self.dd.df.ix[drug].DRUG_NAME
            try:
                res = self._cs_find(drug_name)
            except:
                print("This drug index (%s) / drug name (%s) was not found" %
                        (index, drug_name))
                res = []
            self.results[drug] = res
            pb.animate(i+1)
            results.append(res)
        self.dd_filled.df['CHEMSPIDER_SEARCHED'] = results

    def search_from_smile_inchembl(self):

        N = len(self.drug_ids)

        pb = Progress(N)
        self.results_chembl = {}
        self.results_chemspider = {}

        for i in range(0, N):
            drug = self.drug_ids[i]
            self.results_chembl[drug] = []

            if self.results[drug]:
                for chemspider_id in self.results[drug]:
                    chemspider_entry = self._cs_get(chemspider_id)
                    self.results_chemspider[drug] = chemspider_entry
                    smile = chemspider_entry['smiles']
                    # now search in chembl
                    res_chembl = self.chembl.get_compounds_by_SMILES(smile)
                    try:
                        res_chembl['compounds']
                        self.results_chembl[drug].extend(res_chembl['compounds'])
                    except:
                        pass

            pb.animate(i+1)
