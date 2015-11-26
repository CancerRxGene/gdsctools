GDSCTools documentation
===========================

|version|, |today|

**GDSCTools** is a Python library dedicated to the study of drug responses in the context of the `GDSC (Genomics of Drug Sensitivity in Cancer) <http://www.cancerrxgene.org/>`_ project. It contains utilities to find significant associations between drugs and genomic features (e.g., gene mutation), however, it should be also of interest to a wider community involved in cancer projects.

**GDSCTools** is written in Python. If you are a developer and/or knows 
already about the Python ecosystem and the **pip** command, just type the following command in a :term:`Terminal` to install **GDSCTools**::

    pip install gdsctools

If you are not familiar with this command, please see the :ref:`Installation` section for further details.

In the following example, we provide a short Python snippet on how to use **GDSCTools** (from a programmer point of view). You can either copy and paste the code in a file, and execute it or type the commands in an IPython shell. With this example we investigate the associations between the :term:`IC50` of a given drug (across 52 breast cancer cell lines) and a genomic feature (here, TP53 mutation):

.. plot::
    :width: 80%
    :include-source:

    from gdsctools import ANOVA, ic50_test
    gdsc = ANOVA(ic50_test)
    gdsc.set_cancer_type('breast')
    df = gdsc.anova_one_drug_one_feature('Drug_1047_IC50', 
        'TP53_mut', show=True)

.. seealso:: For more examples and explanations, please visit 
    the :ref:`quickstart` section.


The code above was typed within an IPython shell. Although the code remains
simple, nevertheless, we also provide a standalone application 
called **gdsctools_anova**, which can be used within a standard :term:`Terminal`::

    gdsctools_anova 
        --input-ic50 <ic50 filename> 
        --drug Drug_1047_IC50 
        --feature TP53_mut --onweb


If you want to have a go, please download this 
:download:`IC50 example <../../share/data/test_IC50.csv>`, which is required as
an input. 

.. seealso:: See :ref:`standalone` for more details about the 
    standalone application.


Contents
============


.. toctree::
    :numbered:
    :maxdepth: 1

    installation.rst
    quickstart.rst
    data.rst
    html.rst
    notebooks.rst
    standalone.rst
    settings.rst
    regression.rst
    data_packages.rst
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
