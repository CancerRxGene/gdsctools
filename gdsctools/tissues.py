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


# TCGA notations
cancer_types = {
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


included_in_v17 = ['BLCA', 'BRCA', 'COREAD', 'DLBC', 'ESCA', 'GBM',
                   'HNSC', 'KIRC', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC',
                   'OV', 'PAAD', 'SKCM', 'STAD', 'THCA']


# urogenital_system_other; not prostate, bladder, ovary, cervix
# could be therefore  endometrium (UCEC) 
gdsc_2_tcga = {
        'aero_dig_tract': 'SKCM', # could be also HNSC (head and neck)
        'Bladder': 'BLCA',
        'breast': 'BRCA',
        'cervix': 'CESC',
        'endometrium': 'UCEC',
        'kidney':'KIRC',
        'large_intestine': 'COREAD',
        'lung_NSCLC': 'LUSC',
        'lung_SCLC': 'LUAD', #TODO check that one
        'nervous_system': 'BGM',
        'ovary': 'OV',
        'pancreas':'PAAD',
        'prostate': 'PRAD',
        'skin': 'SKCM', #melanoma
        'soft_tissue': None,
        'stomach': 'STAD',
        'thyroid': 'THCA',
        'biliary_tract': None,
        'bone': None,
        'leukemia': None,
        'liver': None,
        'lymphoma': None,
        'myeloma': None,
        'neuroblastoma': None,
        'testis': None,
        'urogenital_system_other': 'UCEC' #TOD check this
}

from easydev import swapdict
tcga_2_gdsc = swapdict(gdsc_2_tcga, check_ambiguity=False)

