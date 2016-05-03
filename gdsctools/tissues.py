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
        >>> tt = TCGA()
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
        for k, v in _data.items():
            self[k] = v


# GDSC to TCGA is ambiguous. E.g, lung_NSCLC contains LUAD and LUSC labels.
# TCGA to GDSC is also ambiguous since it may be a GDSC1 label (e.g. bone but
# several GDSC2 labels), OR a GDSC2 label.  


#: TCGA keys used in GDSC1000
TCGA_GDSC1000 = ['BLCA', 'BRCA', 'COREAD', 'DLBC', 'ESCA', 'GBM',
    'HNSC', 'KIRC', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC',
    'OV', 'PAAD',  'SKCM', 'STAD', 'THCA']

# Here are some TCGA used in GDSC1000 genomic features 
# and the corresponding Tissue label used in the input file.
TCGA_2_GDSC = {
        'ACC': 'kidney',
        'ALL': 'leukemia',
        #'ALL': 'lymphona', there are ALL with lymphona and leukemia so not a one-one mapping !
        # THere were 15 leukemia and 1 lymphona so we choose here leukemia.
        'BLCA': 'Bladder',
        'BRCA':'breast',
        'CESC': 'cervix',
        'CLL': 'lymphona',
        'COREAD': 'large_instetine',
        'DLBC': 'lymphona',
        'ESCA': 'aero_dig_tract',
        'GBM': 'nervous_system',
        'HNSC': 'aero_dig_tract',
        'KIRC': 'kidney',
        'LAML': 'leukemia',
        'LGG': 'nervous_system',
        'LIHC': 'liver',
        'LUAD': 'lung_NSCLC',
        'LUSC': 'lung_NSCLC',
        'LCML': 'leukemia',
        'MB': 'nervous_system',
        'MESO': 'lung',
        'MM': 'myeloma',
        'NB': 'neuroblastoma',
        'OV': 'ovary',
        'PAAD': 'pancreas',
        'PRAD': 'prostate',
        'SCLC': 'lung_SCLC',
        'SKCM': 'skin',
        'STAD': 'stomach',
        'THCA': 'thyroid',
        'UCEC': 'endometrium',
        }


class Tissues(object):
    """List of tissues included in various analysis
    
    Contains tissues included e.g in v17,v18
    """
    def __init__(self):

        self.v17 = ['myeloma', 'nervous_system', 
            'soft_tissue', 'bone', 'lung_NSCLC',  'skin', 
            'Bladder', 'cervix', 'lung_SCLC', 'lung', 'neuroblastoma',  
            'pancreas', 'aero_dig_tract', 'breast', 'kidney', 'leukemia',
            'ovary', 'prostate', 'large_intestine', 'lymphoma',
            'thyroid', 'stomach', 'biliary_tract', 'endometrium',
            'liver',  'urogenital_system_other', 'testis']
        self.v18 = self.v17



"""
The tissue factor provided in the genomic feature file seem to mix GDSC1 and
GDSC2 labels. For instance lung_NSCLC is actually made of several TCGA labels
such as LUAD, LUSC and NAs. So for now we must use the COSMIC identifiers 
to make sure what were are speaking about. We can not assume that labels are
GDSC1 or GDSC2. 

To help us, we have the COSMICInfo class that contains a bunch of information.
with TCGA and COSMIC Identifiers together + much more.

In v18 TISSUE FACTOR       Count  In Francesco GDSC2 mapping  TCGA   In v18 web
lung_NSCLC                 111     111 GDSC1  64 LUAD, 15 LUSC, 9,  23 NA   YES (both)
leukemia                    82     85 GDSC1   28 LAML, 25 ALL, 10 LCML      LAML only
aero_dig_tract              79     79 GDSC1         42 HNSC + 35 ESCA       YES (both)
lymphoma                    69     70 GDSC1         35 DLBC 3 CLL 1 ALL     YES (DCLC)
lung_SCLC                   64     66 GDSC1         SCLC                    NO
skin                        58     58 GDSC1         55 SKCM 3 NA            YES
nervous_system              56     57 GDSC1         36 GBM 17 LGG 4 MB      YES (GBM
breast                      52     52                 51 BRCA  1NA          YES
large_intestine             50     50                 COREAD                YES
ovary                       43     43               34 OV 3 Others 6 NA     YES (OV)
bone                        39     GDSC1 39                      NA         NO
kidney                      34        34          32 KIRC 1 ACC  1 NA       YES
neuroblastoma               32        32          32 NB                     NO
pancreas                    32        32          30 PAAD + 2NA             YES
stomach                     29        29          25 STAD + 4 NA            YES
lung                        22        GDSC1 22    MESO +1 unclassified      NO
soft_tissue                 21        GDSC1 21                 NAN          NO
Bladder                     19        19                       BLCA         YES
myeloma                     18        18                       MM           NO
liver                       17        17                       LIHC         YES
thyroid                     16        16                       THCA         YES
cervix                      15       15                        CESC         YES
endometrium                 11       11                     9UCEC 2 NA      NO
prostate                     7        8                    6 PRAD 2 NA      NO
biliary_tract                5        5                        NAN          NO
urogenital_system_other      4        4                        NAN          NO
testis                       3        3                        NAN          NO

"""



