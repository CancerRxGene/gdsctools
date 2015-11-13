

.. _quickstart:

Quick Start
=============

There are currently 3 ways to use **GDSCTools**:

#. Typing commands in a **Python shell**. That method is for **developers and users** who want to have high flexibility and use the **GDSCTools** to the best advantage. This documentation provides lots of examples to get you started; the :ref:`references` should also be very informative.
#. Re-using **IPython Notebooks** available within the source code.
#. Using a **standalone application** called **gdsctools_anova**,
   which is the recommended version for **end-users**.

In this section we will focus on the first approach. This will also allow you to use IPython notebooks as explained in :ref:`notebooks` section. The third approach is presented in the :ref:`standalone` section. The standalone application can be used to reproduce the following examples and will be used to produce **data packages**. In the parlance of **GDSCTools**, a data package is the results of an analysis together with an HTML report (see :ref:`data_packages` section).


We assume now that you have **gdsctools** installed together with **IPython**.
If not, please go back to the :ref:`installation` section.

.. note:: Instead of Python, we will use IPython, which is more flexible 
    and interactive that the original version. To start IPython, type this
    command in a terminal/shell::

        ipython

    You should now see something like::

        Python 2.7.5 (default, Nov  3 2014, 14:33:39) 
        Type "copyright", "credits" or "license" for more information.

        IPython 4.0.0 -- An enhanced Interactive Python.
        ?         -> Introduction and overview of IPython's features.
        %quickref -> Quick reference.
        help      -> Python's own help system.
        object?   -> Details about 'object', use 'object??' for extra details.

        In [1]: 




.. note:: All snippets in this documentation are typed within IPython shell. 
    You may see >>> signs. They indicate a python statement typed in 
    a shell. Lines without those signs indicate the output of the previous
    statement. For instance::

        >>> a = 3
        >>> 2 + a
        5

    means the code **2 + a** should print the value 5

The IC50 input data 
-------------------------------

Before starting, we first need to get a data file to play with. The only
required data file is one that contains :term:`IC50`.  
For now, we do not need to enter into the details of the expected data
structure; it should be a CSV or TSV file as in this :download:`IC50 example <../../share/data/test_IC50.csv>` file.

.. seealso:: More details about the data format can be found in the :ref:`data` section as well as links to retrieve IC50 data sets.


Although all functionalities could be imported using::

    from gdsctools import *

we will be as explicit as possible in the following examples; we would rather use (for instance)::

    from gdsctools import IC50

This is better coding practice and has also the advantage of telling beginners
which functions are going to be used. 

Here above, we imported the :class:`~gdsctools.readers.IC50` class, that allows one to read the IC50 example file aforementioned. We will explain in details the different data sets and their formats in the :ref:`data` section. However, for now it is enough to know that it should be a CSV formatted file that contains IC50s; one value for each combination of drug and cell line. 

Note that the IC50 example is also installed with **GDSCTools** and its location obtained using::

    from gdsctools import ic50_test
    print(ic50_test.filename)

The :class:`~gdsctools.readers.IC50` class is flexible enough that you can provide the filename location or just the name **ic50_test** as in the example below:: 

    >>> from gdsctools import IC50, ic50_test
    >>> ic = IC50(ic50_test)
    >> print(ic)
    Number of drugs: 11
    Number of cell lines: 988
    Percentage of NA 0.206569746043

As you can see you can get some information about the IC50 content (e.g., 
number of drugs, percentage of NaNs) using the :ref:`print` function. See :class:`gdsctools.readers.IC50` for more details.

Getting help
---------------

At any time, you can get help about a **GDSCTools** functionality or a python function by adding question tag after a function's name::

    IC50?

With **GDSCTools**, we also provide a convenience function called :func:`~gdsctools.gsdctools_help`::

    gdsctools_help(IC50)

that should open a new tab in a browser redirecting you to the HTML help version (on `ReadTheDoc website <gdsctools.readthedocs.org>`_) of a function or class (here the :class:`IC50` class).


    

The ANOVA class
----------------
Given an IC50 data set, we can now analyse it using the main class 
called :class:`~gdsctools.anova.ANOVA`. 


Before starting, just a few words about the underlying stastistical analysis. On one hand, we have an IC50 file. It contains IC50s measured for :math:`N_d` drugs and :math:`N_c` cell lines. Each combination of drug and cell line has a unique measured IC50. On the other hand, one should also provide a data file with genomic features with the same set of :math:`$N_c$` cell lines. The other dimension being the :math:`N_f` genomic features (e.g. mutation). A default set of 680 genomic features is used so we will not give more information about that data file for now. 

.. seealso:: More details about the genomic features data format can be found in the :ref:`data` section.

For each drug, the ANOVA class computes the association (a regression analysis) between a drug and a feature leading to a p-value. This calculation is repeated across all features and all drugs. Consequently, a multiple testing correction is applied and reported in the analysis. For more information, please see the :ref:`userguide` section.

Let us now proceed to the analysis using the IC50 test example::

    from gdsctools import ANOVA, ic50_test
    gdsc = ANOVA(ic50_test)

You may provide a variety of input file here (in CSV or TSV format). 
You can now analysis the entire IC50 data across all drugs and all features but
you can actually perform three different kind of analysis::

.. index:: ODOF, ODAF, ADAF

#. compute all associations for all drugs and all features. (ADAF)
#. compute the associations between one drug and all the features (ODAF)
#. compute one association between a drug and feature (ODOF)


One Drug One Feature (ODOF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Let us start with the first case. You can analyse a given drug for 
a given genomic feature using the
:meth:`~gdsctools.anova.ANOVA.anova_one_drug_one_feature` method:

.. plot::
    :include-source:

    from gdsctools import ANOVA, ic50_test
    gdsc = ANOVA(ic50_test)
    gdsc.anova_one_drug_one_feature('Drug_999_IC50', 'TP53_mut', 
        show=True)


.. todo:: explain the analysis and the plots

.. seealso:: a link to more explanation

One Drug All Features (ODAF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a similar way, you may look at all features for a given drug:

.. plot::
    :include-source:

    from gdsctools import ANOVA, ic50_test
    gdsc = ANOVA(ic50_test)
    df = gdsc.anova_one_drug('Drug_999_IC50')  
    
    # no plots were generated in the previous statement
    from gdsctools import VolcanoANOVA
    df = gdsc.add_pvalues_correction(df)
    v = VolcanoANOVA(df)
    v.volcano_plot_all()

.. note:: When you call the ODAF method, you are actually calling
   the ODOF method for each feature. This method takes 4-10 seconds 
   per drug depending on the number of features.


.. todo:: explain the analysis and the plots
.. seealso:: a link to more explanation

All Drug All Features (ADAF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Or analyse a all drugs across all features. This takes a long depending on the
number of drugs and features (30 minutes for 250 drugs and 1000 features):

.. plot::
    :include-source:

    from gdsctools import ANOVA, ic50_test
    gdsc = ANOVA(ic50_test)
    gdsc.set_cancer_type('breast')
    results = gdsc.anova_all()

    from gdsctools import VolcanoANOVA
    v = VolcanoANOVA(results.df)
    v.volcano_plot_all()

.. note:: When you call the :meth:`anova_all` method (ADAF) you are
    actually calling the :meth:`anova_one_drug` for each drug. 
    
.. warning:: :meth:`anova_all` may take a long time to run 
    (e.g., 10 minutes, 30 minutes) depending on the number of drugs
    and features.

.. todo:: explain the analysis and the plots
.. todo:: FDR threshold to show some green/red dots
.. seealso:: a link to more explanation

HTML report
==============

You can also create a thorough HTML report 
::

    >>> from gdsctools import ANOVA, ic50_test
    >>> gdsc = ANOVA(ic50_test)
    >>> gdsc.set_cancer_type('breast')
    >>> results = gdsc.anova_all()
    >>> report = ANOVAReport(gdsc, results)



