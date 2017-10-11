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

:Note: Developed and tested for Python 2.7, 3.5, 3.6
:Note: The GDSCTools libary works for Python 2.7 and 3.X but the standalone
       pipeline to be ran on cluster works on Python 3.X only (requires
       Snakemake).
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
