# -*- python -*-
# -*- coding utf-8 -*-

#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.com>
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
from easydev import AttrDict


class TCGA(AttrDict):
    """A dictionary-like container to access to TCGA cancer types keywords
    
    .. doctest::

        >>> from gdsctools import TCGA
        >>. tt = TCGA()
        >>> tt.ACC
        'Adrenocortical Carcinoma'


    """
    def __init__(self):
        super(TCGA, self).__init__()
        # TCGA notations
        _data = {
            'ACC': 'Adrenocortical Carcinoma',
            'ALL': 'Acute Lymphoblastic Leukemia',
            'BLCA': 'Bladder Urothelial Carcinoma',
            'BRCA': 'Breast invasive Carcinoma',
            'CESC': 'Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma',
            'CLL': 'Chronic Lymphocytic Leukemia',
            'COAD': 'Colon Adenocarcinoma',
            'COREAD': 'Colon Adenocarcinoma and Rectum Adenocarcinoma',
            'DLBC': 'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma',
            'ESCA': 'Esophageal Carcinoma',
            'GBM': 'Glioblastoma multiforme',
            'HNSC': 'Head and Neck Squamous Cell Carcinoma',
            'KICH': 'Kidney Chromophobe',
            'KIRC': 'Kidney Renal Clear Cell Carcinoma',
            'KIRP': 'Kidney Renal Papillary Cell Carcinoma',
            'LAML': 'Acute Myeloid Leukemia',
            'LCML': 'Chronic Myelogenous Leukemia',
            'LCLL': 'Chronic Lymphocytic Leukemia',
        #    'LCMM': 'TOCHECK lymphoid_neoplasm_other',
            'LGG': 'Brain Lower Grade Glioma',
            'LIHC': 'Liver Hepatocellular Carcinoma',
            'LUAD': 'Lung Adenocarcinoma',
            'LUSC': 'Lung Squamous Cell Carcinoma',
            'MB':  'Medulloblastoma',
            'MESO': 'Mesothelioma',
            'MM': 'Multiple Myeloma',
            'NB': 'Neuroblastoma',
            'OV': 'Ovarian Serous Cystadenocarcinoma',
            'PAAD': 'Pancreatic Adenocarcinoma',
            'PCPG': 'Pheochromocytoma and Paraganglioma',
            'PRAD': 'Prostate Adenocarcinoma',
            'READ': 'Rectuem Adenocarcinoma',
            'SCLC': 'Small Cell Lung Cancer',
            'SKCM': 'Skin Cutaneous Melanoma',
            'STAD': 'Stomach Adenocarcinoma',
            'THCA': 'Thyroid Carcinoma',
            'UCS': 'Uterine Carcinosarcoma',
            'UCEC': 'Uterine Corpus Endometriod Carcinoma',
            }
        for k,v in _data.items():
            self[k] = v


#: TCGA keys used in GDSC1000
TCGA_GDSC1000 = ['BLCA', 'BRCA', 'COREAD', 'DLBC', 'ESCA', 'GBM',
    'HNSC', 'KIRC', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC',
    'OV', 'PAAD', 'SKCM', 'STAD', 'THCA']


# urogenital_system_other; not prostate, bladder, ovary, cervix
# could be therefore  endometrium (UCEC) 
GDSC1_2_TCGA = {
        #'aero_dig_tract': 'SKCM', # could be also HNSC (head and neck)
        'aero_dig_tract': 'ESCA',
        'bone': None,
        'Bladder': 'BLCA',
        'breast': 'BRCA',
        'cervix': 'CESC',
        'endometrium': 'UCEC',
        'kidney':'KIRC',
        #'large_intestine': 'COREAD', GDSC2
        'disgestive_system': 'COREAD',
        'lung_NSCLC': 'LUSC',
        'lung_SCLC': 'LUAD', #TODO check that one
        'nervous_system': 'GBM',
        #'nervous_system': 'LGG',
        'ovary': 'OV',
        'pancreas':'PAAD',
        'prostate': 'PRAD',
        'skin': 'SKCM', #melanoma
        'soft_tissue': None,
        'stomach': 'STAD',
        'thyroid': 'THCA',
        'biliary_tract': None,
        'leukemia': 'LAML',
        'liver':' LIHC',
        #'lymphoma': DLBC,
        #'blood': DLBC,
        'myeloma': None,
        'neuroblastoma': None,
        'testis': None,
        'urogenital_system_other': 'UCEC' #TOD check this
}

#from gdsctools import COSMICInfo
#c = COSMICInfo()
#goups = c.df.groupby('TCGA')

mapping = [
    (['leukemia', 'lymphona'], ['leukemia', 'lymphoblastic_leukemia', 
        'lymphoblastic_T_cell_leukaemia', 'lymphoid_neoplasm other', 
        'T_cell_leukemia'], 'ALL'),
    ('lymphoma', 'lymphoid_neoplasm other',  'CLL'),
    ('lymphoma' 'B_cell_lymphoma' 'DLBC'),
    ('aero_dig_tract' 'head and neck' 'HNSC'),
    ('leukemia',['acute_myeloid_leukaemia', 'leukemia'], 'LAML'),
    ('leukemia' 'chronic_myeloid_leukaemia' 'LCML'),
    ('nervous_system' 'glioma' 'LGG'),
    ('lung_NSCLC' 'lung_NSCLC_adenocarcinoma' 'LUAD'),
    ('lung_NSCLC' 'lung_NSCLC_squamous_cell_carcinoma' 'LUSC'),
    
    ('myeloma', ['haematopoietic_neoplasm other',
      'lymphoid_neoplasm other', 'myeloma'], 'MM'),
    
    
    ('kidney', 'adrenal_gland', 'ACC'),
    ('breast',  'breast',  'BRCA'),
    ('urogenital_system',  'Bladder',  'BLCA'),
    ('urogenital_system',  'cervix',  'CESC'),
    ('aero_dig_tract' 'oesophagus' 'ESCA'),
    ('nervous_system' 'glioma' 'GBM'),
    ('kidney' 'kidney' 'KIRC'),
    ('digestive_system' 'liver' 'LIHC'),
    
    ('skin', 'melanoma', 'SKCM'),
    ('digestive_system', 'stomach', 'STAD'),
    ('large_intestine',  'large_intestine',  'COREAD'),
    ('nervous_system' 'medulloblastoma' 'MB'),
    ('lung' 'mesothelioma' 'MESO'),
    ('neuroblastoma', 'neuroblastoma', 'NB'),
    ('urogenital_system', 'ovary', 'OV'),
    ('pancreas', 'pancreas', 'PAAD'),
    ('urogenital_system', 'prostate', 'PRAD'),
    ('lung_SCLC' 'lung_small_cell_carcinoma' 'SCLC'),
    ('thyroid', 'thyroid', 'THCA'),
    ('urogenital_system', 'endometrium' 'UCEC')
]

# TO check:
# 908118, 1330993
# 1240141 


GDSC1_included_in_v17 = ['myeloma', 'nervous_system', 'soft_tissue', 'bone', 'lung_NSCLC',  'skin', 'Bladder', 'cervix', 'lung_SCLC', 'lung', 'neuroblastoma',  'pancreas', 'aero_dig_tract', 'breast', 'kidney', 'leukemia',
        'ovary', 'prostate', 'large_intestine', 'lymphoma',
        'thyroid', 'stomach', 'biliary_tract', 'endometrium',
    'liver',  'urogenital_system_other', 'testis']

#GDSC1_included_in_v17 = ['myeloma', 'soft_tissue', 'bone', 'lung_NSCLC', 'cervix', 'lung_SCLC', 'lung', 'neuroblastoma',  
#        'prostate',    'biliary_tract', 'endometrium',
#     'urogenital_system_other', 'testis']

map_to_tcga = {
        'large_intestine':'COREAD',
        'lymphoma': 'DLBC',
        'Bladder': 'BLCA',
        'breast' :'BRCA',
        'aero_dig_tract': 'ESCA',
        'nervous_system': 'GBM',
        '??': 'HNSC', # stored as EScA ???
        'kidney': 'KIRC',
        'leukemia': 'LAML',
        '???': 'LGG', # stored as GBM ???
        'liver': 'LIHC',
        'lung_NSCLC_error2': 'LUAD',
        'lung_NSCLC_error1': 'LUSC',
        'ovary': 'OV',
        'pancreas': 'PAAD',
        'skin': 'SKCM',
        'stomach': 'STAD',
        'thyroid': 'THCA'
        }

# DLBC in PANCAN --> 69 but in DLBC file, only 34 ?
# ESCA in PANCAN --> 79 but in DLBC file, only 35 ?
# other discprepancies: GBM,



# Takes COSMICFetcher, get df and create groups by TCGA
# cf = COSMICFetcher()
# cf = cf.df.set_index('COSMIC_ID')
# cf.ix[groups['PRAD']]

# c= COSMICInfo()
# groups = c.df.groupby('TCGA').groups
# Does that make sense ? 



# OK: SCLCTHCA,UCEC
# ??: STAD

#UNABLE TO CLASSIFY [['lung_NSCLC' 'lung_NSCLC_not specified' 'UNABLE TO
#    CLASSIFY']
#     ['lung' 'Lung_other' 'UNABLE TO CLASSIFY']
#     ['pancreas' 'pancreas' 'UNABLE TO CLASSIFY']
#       ['urogenital_system' 'ovary' 'UNABLE TO CLASSIFY']]
#
#   ('lymphoma','B_cell_lymphoma','DLBC'),



    
    




#from easydev import swapdict
#tcga_2_gdsc = swapdict(gdsc_2_tcga, check_ambiguity=False)

