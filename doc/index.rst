GDSCTools documentation
===========================

|version|, |today|


.. image:: https://badge.fury.io/py/gdsctools.svg
    :target: https://pypi.python.org/pypi/gdsctools

.. image:: https://secure.travis-ci.org/CancerRxGene/gdsctools.png
    :target: http://travis-ci.org/CancerRxGene/gdsctools

.. image::  https://coveralls.io/repos/CancerRxGene/gdsctools/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/CancerRxGene/gdsctools?branch=master

.. image:: https://readthedocs.org/projects/gdsctools/badge/?version=master
    :target: http://gdsctools.readthedocs.io/en/master/?badge=master
    :alt: Documentation Status


:Citation: Cokelaer et al. GDSCTools for mining pharmacogenomic interactions in cancer
    Bioinformatics, 2017, https://doi.org/10.1093/bioinformatics/btx744
:Note: developed and tested for Python 3.7, 3.8, 3.9
:Contributions: Please join https://github.com/CancerRxGene/gdsctools project
:Documentation: `On ReadTheDocs <http://gdsctools.readthedocs.io/en/master>`_
:GitHub: `On github <https://github.com/CancerRxGene/gdsctools/issues>`_



**GDSCTools** is a free open-source Python library dedicated to the study of drug
responses in the context of the `GDSC (Genomics of Drug Sensitivity in Cancer) <http://www.cancerrxgene.org/>`_ project.  The main developer is `Thomas Cokelaer (Institut Pasteur) <https://research.pasteur.fr/en/member/thomas-cokelaer/>`_, and it is a
joint effort of the groups of `Mathew Garnett (Sanger Institute) <http://www.sanger.ac.uk/people/directory/garnett-mathew>`_ and `Julio Saez-Rodriguez (RWTH Aachen & EMBL-EBI) <http://www.combine.rwth-aachen.de/index.php/people/julio-saez-rodriguez.html>`_.

It contains utilities to find significant associations between drugs and genomic
features (e.g., gene mutation) based on an :ref:`ANOVA <anova_partone>` analysis. Other methods, such
as multi-factorial linear models based on :ref:`Elastic Net <multivariate_regression>` are also available.
Besides, the library should also be useful for manipulating dedicated data sets
such as :ref:`IC50 <data>` (drug response) or :ref:`MoBEM <omnibem>` (genomic features) data structures. Hence,
we hope that GDSCTools serves as basis for other scientists to develop further
methods.


.. |association| image::
    auto_examples/images/sphx_glr_plot_association_001.png
    :target: auto_examples/plot_association.html

.. |volcano| image::
    auto_examples/images/sphx_glr_plot_volcano_001.png
    :target: auto_examples/plot_volcano.html

.. |histo| image::
    auto_examples/images/sphx_glr_plot_ic50_hist_001.png
    :target: auto_examples/plot_ic50_hist.html

.. |elastic| image::
    auto_examples/images/sphx_glr_plot_elastic_tuning_001.png
    :target: auto_examples/plot_elastic_tuning.html


.. |elastic_weights| image::
    auto_examples/images/sphx_glr_plot_elastic_weights_001.png
    :target: auto_examples/plot_elastic_weights.html


.. raw:: html

    <div class="body">
    <div id="index-grid" class="section group">
    <div class="col span_1_of_3">
        <h3><a href="quickstart.html">First Steps</a></h3>
        <p>Get started with GDSCTools</p>
        <h3><a href="auto_examples/index.html">Examples</a></h3>
        <p>Visit our example gallery</p>
        <h3><a href="anova_partone.html">The ANOVA analysis</a></h3>
        <p>Browse the full documentation</p>
    </div>
    <div class="col span_2_of_3">
    <div class="jcarousel-wrapper">
    <div class="jcarousel">

* |association|
* |volcano|
* |histo|
* |elastic|
* |elastic_weights|

.. raw:: html

            </div>
        <a href="#" class="jcarousel-control-prev">&lsaquo;</a>
        <a href="#" class="jcarousel-control-next">&rsaquo;</a>
        <p class="jcarousel-pagination">
        </p>
        </div>
        </div>
        </div>
   </div>
   <div style="clear: left"></div>


.. index:: installation

**GDSCTools** is written in Python. If you are a developer and/or knows
already about the Python ecosystem and the **pip** command, just type 
the following command in a :term:`Terminal` to install **GDSCTools**::

    pip install gdsctools

add the option ``--upgrade`` to get the latest release. Conversely, if you are not
familiar with Python or the command above, please see the :ref:`Installation` section
for further details. Note also that we strongly recommend to use Anaconda to install dependencies (e.g., numpy, matplotlib); **GDSCTools** is available on bioconda channel::

    conda install gdsctools.

In the following example, we provide a short Python snippet that uses the **GDSCTools** library. You can either copy and paste the code in a file, and execute it or type the commands in an :term:`IPython` shell. With this example we investigate the associations between the :term:`IC50` of a given drug (across 52 breast cancer cell lines) and a genomic feature (here, TP53 mutation). Drugs are refer to by a unique identifier (here 1047):


.. doctest::

    from gdsctools import ANOVA, ic50_test
    gdsc = ANOVA(ic50_test)
    gdsc.set_cancer_type('breast')
    df = gdsc.anova_one_drug_one_feature(1047, 'TP53_mut', show=True)


.. plot::
    :width: 80%

    from gdsctools import ANOVA, ic50_test
    gdsc = ANOVA(ic50_test)
    gdsc.set_cancer_type('breast')
    df = gdsc.anova_one_drug_one_feature(1047, 'TP53_mut', show=True)

The :attr:`df` object returned in the last statement is a dataframe that contains information explained in :ref:`regression` section.


.. seealso:: For more examples and explanations, please visit
    the :ref:`anova_partone` section.


.. index:: warnings

The previous example may be verbose with comments and warnings. You may
set the verbose option to False and ignore warnings as follows::

    import warnings
    warnings.simplefilter("ignore","exceptions.Warning")
    gdsc = ANOVA(ic50_test, verbose=False)


.. index:: standalone

We will see more examples on how to use **GDSCTools** to perform more
systematic studies. However, let us note that **GDSCTools** also provide
a standalone application called **gdsctools_anova**, 
which can be used within a standard :term:`Terminal` (same 
output as in the previous example)::

    gdsctools_anova --input-ic50 <ic50 filename> --drug 1047
        --feature TP53_mut

If you want to have a go, please download this
:download:`IC50 example <../gdsctools/data/test_IC50.csv>`, which is required as an input.
Note that by default, **GDSCTools** loads a set of 50 genomic features and 1001 cell
lines but in general, you should provide your own genomic feature file (see :ref:`data`). The default data set contains only a small set of genomic features and can be downloaded:
:download:`GenomicFeature example <../gdsctools/data/genomic_features.tsv.gz>`, and adapted to your needs.


.. seealso:: See :ref:`standalone` section for more details about the
    standalone application and the :ref:`data` section to learn more
    about the expected input data formats.


Contents
============


.. toctree::
    :numbered:
    :maxdepth: 1

    installation.rst
    quickstart.rst
    data.rst
    anova_partone.rst
    anova_parttwo.rst
    html.rst
    data_packages.rst
    omnibem.rst
    regression.rst
    notebooks.rst
    standalone.rst
    auto_examples/index
    releases.rst
    references.rst
    developers.rst

Issues
===========

Please fill bug report in https://github.com/CancerRxGene/gdsctools/issues

Contributions
================

Please join https://github.com/CancerRxGene/gdsctools

.. toctree::
    :hidden:

    ChangeLog.rst
    faqs
    glossary

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
