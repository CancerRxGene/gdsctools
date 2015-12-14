from gdsctools.anova import ANOVA
from gdsctools.readers import IC50
from gdsctools.tools import extract_drug_identifiers
import pandas as pd


class IC50Cluster(IC50):
    def __init__(self, ic50, ratio_threshold=10, verbose=True):
        """

        """
        super(IC50Cluster, self).__init__(ic50)
        self.verbose = verbose
        self.ratio_threshold = ratio_threshold

        if self.verbose:
            print(self)
        self.cluster()
        if self.verbose:
            print(self)

    def _get_to_cluster(self):
        info = self.info()
        to_cluster = info[info.ratio >= 10].DRUG_ID.values
        return list(to_cluster)
    to_cluster = property(_get_to_cluster)

    def _get_mapping(self):
        from collections import defaultdict
        mapping = defaultdict(list)

        drug_ids = extract_drug_identifiers(self.df.columns)
        for drug_id, colname in zip(drug_ids, self.df.columns):
            mapping[drug_id].append(colname)
        return mapping

    def _get_duplicated(self):
        mapping = self._get_mapping()
        duplicated = [key for key in mapping.keys() if len(mapping[key])>1]
        return duplicated
    duplicated = property(_get_duplicated)


    def info(self):
        mapping = self._get_mapping()
        duplicated = self.duplicated

        max_ids = max([len(x) for x in mapping.values()])

        results = []
        for drug_id in duplicated:
            df = self.df[mapping[drug_id]]

            total = df.mean(axis=1).count()
            individuals = list(df.count().values)

            # the number of individuals may be different. In v18 one drug had 3
            # entries (all other duplicated had only 2). so, we need to add
            # dummies (NA) when required:
            individuals += [None] * (max_ids - len(individuals))

            common = sum(df.count(axis=1)>=2)
            result = [drug_id] + individuals + [total, common,
                    100*common/float(total)]
            results.append(result)

        df = pd.DataFrame(results)
        df.columns = ['DRUG_ID'] + [str(x) for x in range(1, max_ids+1)] +\
            ['total', 'common', 'ratio']
        df.sort_values(by='ratio', ascending=False, inplace=True)
        return df

    def cluster(self):
        # get list of drug identifiers to cluster
        to_cluster = self.to_cluster
        mapping = self._get_mapping()
        
        if self.verbose:
            print('Found  %s non unique drug identifiers ' %
                    len(self.duplicated))
            print('Clustering %s of them.\n' % len(self.to_cluster))
        if len(self.to_cluster) == 0:
            return 

        for identifier in to_cluster:
            drug_names = mapping[identifier]
            new_drug_name = str(identifier) + "_" + \
                    "_".join([x.split("_",1)[1] for x in  drug_names])
            # add new column with new name and mean of the columns with same
            # drug id
            self.df[new_drug_name] = self.df[drug_names].mean(axis=1)
            # Remove the individual columns 
            self.df.drop(drug_names, axis=1, inplace=True)





class GDSC(ANOVA):
    """Alias to ANOVA class with default settings 
    
    It also converts tissue names into TCGA names.
    """

    def __init__(self, ic50, genomic_features,
        drug_decode, verbose=True, low_memory=False):

        super(GDSC, self).__init__(ic50, genomic_features=genomic_features, 
            drug_decode=drug_decode, verbose=True, low_memory=low_memory)
        from gdsctools.tissues import map_to_tcga
        def tcga_names(x):
            if x in map_to_tcga.keys():
                return map_to_tcga[x]
            else:
                return x
        newnames =  self.features.df.TISSUE_FACTOR.apply(
                            lambda x: tcga_names(x))

        self.features.df['TISSUE_FACTOR'] = newnames


class DataPackagesAnalysis(ANOVA):


    def __init__(self, ic50, genomic_features,
        drug_decode, company, verbose=True, low_memory=False):

        super(DataPackagesAnalysis, self).__init__(ic50, genomic_features=genomic_features, 
            drug_decode=drug_decode, verbose=True, low_memory=low_memory)
        
        
        self.companies = sorted(self.drug_decode.df.groupby('OWNED_BY').groups.keys())


        self.company = company

        drugs_az = dd.df.query("WEBRELEASE=='Y' or OWNED_BY=='%s'" %
                    'AZ')


        # FOR v17, v18
        # Based on the GF for a given cancer type, prune IC50
        # e.g, in COREAD, IC50 is now made of 50 rows (cell lines)
        # Based on company and public compounds, keep only subset of the 
        # columns

        # in drugs, it uses DRUG_ID
        # in ic50 from v18/francesco, it uses Drug_<ID>_<something>_IC50
        # In Howard IC50, only <ID>_<something>
        ic50.df.columns = ['Drug_'+str(x)+'_IC50' for x in ic50.df.columns]

        gf = GenomicFeatures('../../gdsctools_private/v17/COREAD/DATA/INPUT/ANOVA_input.txt')

        #gf = GenomicFeatures('test/IC50_v18.csv')

        # filter cosmic ID from GF
        ic50.df = ic50.df.ix[gf.cosmicIds]

        # We should keep only those ones
        drug_ids = [x for x in drugs_az.index]


        filtered = []
        for this in ic50.df.columns:
            ID = int(this[5:].split("_")[0])
            if ID in drug_ids:
                filtered.append(True)
            else:
                filtered.append(False)

        ic50.df = ic50.df[ic50.df.columns[filtered]]
        # end up with 407 columns. In v18, only 391

        ic50.df[ic50.df.columns[ic50.df.count()>=6]]
        # In the ANOVAResults from COREAD, v18 379 unique DRUG ID, 381 unique DRUG_CODE
        # In the IC50 from coread, v18 391 unique DRUG_CODE but 381 with at
        # least 6 non NA values and 379 unique DRUG ID


        # In this file, similarly 379 unique DRUG ID but 407 entries













