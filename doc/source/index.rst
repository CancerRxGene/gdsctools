GDSCTools documentation
===========================



**GDSCTools** is a Python library dedicated to the study of drug responses in the context of the `GDSC (Genomics of Drug Sensitivity in Cancer) <http://www.cancerrxgene.org/>`_ project. It contains utilities to find significant association between drugs and genomic features (e.g., gene mutation), however, it should be of interest to a wider community.

**GDSCTools** is written in Python. If you are a developer and knows about
Python ecosystem already, then just type this command to install the library::

    pip install gdsctools

If you are not familiar with this command, please see the :ref:`Installation` section for further details.

Here is a quick example on how to use **GDSCTools** to investigate the association between the IC50 of a given drug and a genomic feature (TP53 mutation) across 52 breast cancer cell lines.


.. plot::
    :include-source:
    :width: 80%

    from gdsctools import *
    gdsc = ANOVA(ic50_test)
    gdsc.set_cancer_type('breast')
    df = gdsc.anova_one_drug_one_feature('Drug_1047_IC50', 
        'TP53_mut', show=True)

The code above needs to be typed within an IPython shell or a notebook. We also
provide a standalone application called **gdsctools_anova**, which can be used
within a shell::

    gdsctools_anova 
        --I <ic50 filename> 
        --drug Drug_1047_IC50 
        --feature TP53_mut --onweb


See :ref:`standalone` for more details.


Contents
============


.. toctree::
    :numbered:
    :maxdepth: 1

    installation.rst
    userguide.rst
    data.rst
    standalone.rst
    references.rst
    developers.rst
    ChangeLog.rst

Issues
===========

Please fill bug report in https://github.com/CancerRxGene/gdsctools/issues

Contributions
================

Please join https://github.com/CancerRxGene/gdsctools

.. toctree::
    :hidden:

    glossary


