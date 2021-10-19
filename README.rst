GDSCTools
============


.. image:: https://badge.fury.io/py/gdsctools.svg
    :target: https://pypi.python.org/pypi/gdsctools

.. image:: https://secure.travis-ci.org/CancerRxGene/gdsctools.png
    :target: http://travis-ci.org/CancerRxGene/gdsctools

.. image::  https://coveralls.io/repos/CancerRxGene/gdsctools/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/CancerRxGene/gdsctools?branch=master

.. image:: https://readthedocs.org/projects/gdsctools/badge/?version=master
    :target: http://gdsctools.readthedocs.io/en/master/?badge=master
    :alt: Documentation Status

|License|

:Citation: Cokelaer et al. GDSCTools for mining pharmacogenomic interactions in 
    cancer.  Bioinformatics, 2017, https://doi.org/10.1093/bioinformatics/btx744

:Note: Developed and tested for Python 3.7, 3.8, 3.9
:Contributions: Please join https://github.com/CancerRxGene/gdsctools project
:Documentation: `On ReadTheDocs <http://gdsctools.readthedocs.io/en/master>`_
:GitHub: `On github <https://github.com/CancerRxGene/gdsctools/issues>`_

Overview
-----------
Genomics of Drug Sensitivity in Cancer (GDSC) tools including pipelines related  to http://www.cancerrxgene.org/

Installation
---------------

::

  pip install gdsctools

For beginners, please visit the main documentation `Installation
<http://gdsctools.readthedocs.io/en/master/installation.html>`_ section.


QuickStart
-------------

You will need 2 input matrices:

#. an IC50 matrix with COSMIC identifiers as rows and drugs as columns,
#. a genomic feature matrix with COSMIC identifiers as rows and features as columns.

Then, you can analyse the data with the standalone application::

    gdsctools_anova --input-ic50 ic50.txt --input-features features.txt

or as a script::

  from gdsctools import anova, ic50_test
  an = ANOVA(ic50_test, features_filename)  # second arg is optional
  an.anova_all()

More examples are provided in the `documentation on ReadThedoc <http://gdsctools.readthedocs.io/en/master/index.html>`_.

Note that first versions (ANOVA analysis) were based on https://github.com/francescojm/FI.GDSC.ANOVA repository. New tools have been added (regression based on Ridge, Lasso, OmniBEM tool, ...).


.. |License| image:: https://img.shields.io/badge/license-BSD-blue.svg
   :alt: BSD License
   :target: https://github.com/CancerRxGene/gdsctools/blob/master/LICENSE


Changelog
~~~~~~~~~

Prior to v0.20, please see the online documentation (link above)

========= ====================================================================
Version   Description
========= ====================================================================
1.1.0     * Update GDSCtools to newest Python ecosytem (3.9)
          * add CI
          * remove redundant dependencies (biokit)
1.0.1     * fix the MANIFEST to include requirements.txt
1.0.0     * First release of GDSCTools with test coverage, singularity, CI
0.20.0    * Fixing dependencies. added Snakemake example
========= ====================================================================

