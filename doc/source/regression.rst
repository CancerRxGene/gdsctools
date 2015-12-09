
.. _regression:

About the regression and ANOVA analysis
==========================================



By default, the regression uses a :term:`OLS` method. 4 Factors may be taken
into account depending on the content of the 
:class:`~gdsctools.readers.GenomicFeature` data set.

Here, we will use R syntax for simplicity. Depending on the data and the
:class:`gdsctools.anova.settings` one of the following formula will be used:


.. math:: Y \sim C(tissue) + C(media) + C(MSI) + Feature


.. math:: Y \sim C(tissue) + C(MSI) + Feature


.. math:: Y \sim C(MSI) + Feature


.. math:: Y \sim Feature



