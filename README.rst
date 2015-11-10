GDSCtools 
============


.. image:: https://badge.fury.io/py/gdsctools.svg
    :target: https://pypi.python.org/pypi/gdsctools

.. image:: https://secure.travis-ci.org/CancerRxGene/gdsctools.png
    :target: http://travis-ci.org/CancerRxGene/gdsctools

.. image::  https://coveralls.io/repos/CancerRxGene/gdsctools/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/CancerRxGene/gdsctools?branch=master 

.. image:: https://badge.waffle.io/CancerRxGene/gdsctools.png?label=issues&title=issues
   :target: https://waffle.io/gdsctools/gdsctools

:Note: tested for Python 2.7, 3.3, 3.4
:Contributions: Please join https://github.com/CancerRxGene/gdsctools project

:Online documentation: `On pypi website <http://pythonhosted.org/gdsctools/>`_

:Issues and bug reports: `On github <https://github.com/CancerRxGene/gdsctools/issues>`_

Overview
-----------
Genomics of Drug Sensitivity in Cancer (GDSC) tools including pipelines related  to http://www.cancerrxgene.org/ 

Installation
---------------

::

  pip install gdsctools
  
  
QuickStart
-------------

2 inputs matrices are required: (1)  an IC50 matrix with COSMIC identifiers as rows and drugs as columns, and (2) a genomic feature matrix with COSMIC identifiers as rows and features as columns. Then, you can analyse the data with the standalone application::

    gdsctools_anova --ic50 ic50.txt --feature features.txt --analyse-all 

or as a script::

  from gdsctools import anova, reader
  ic50 = reader.IC50('ic50.txt')
  features = reader.Features('features.txt')
  an = ANOVA(ic50, features)
  an.anova_all()
  
  
More examples are provided in the documentation. You can for example select a specific drug, or a set of drugs instead of the entire screening, or perform a sub selection on features. 






.. note:: Version 1 (linear regression and ANOVA) was created based on https://github.com/francescojm/FI.GDSC.ANOVA
