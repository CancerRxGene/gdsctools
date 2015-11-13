

.. _quickstart:

Quick Start
=============

There are currently 3 ways to use **GDSCTools**:

#. From a **Python shell**. That method is for **developers and users** who want to have high flexibility and use the **GDSCTools** to the best advantage. This documentation provides lots of examples to get you started; the :ref:`references` should also be very informative.
#. Using **Python Notebooks** available within the source code. 
#. From a shell, using a standalone application called **gdsctools_anova**,
   which is the recommended version for **end-users**.

In this section we will focus on the first approach, which would allow you to
use the second approach is very familiar and introduced in the :ref:`notebooks` section. The third approach is presented in the :ref:`standalone` section. The standalone application can be used to reproduce the following examples and produce **data packages**. In the parlance of **GDSCTools**, a data package is the results of an analysis together with HTML report (see :ref:`data_packages` section).


We assume now that you have **gdsctools** installed together with **IPython**.
If not, please go back to the :ref:`installation` section.


.. note:: snippets here below are typed within a IPython shell. 
    You may see >>> signs. They indicate a python statement typed in 
    a shell. Lines without those signs indicate the output of the previous
    statement. For instance::

        >>> a = 3
        >>> 2 + a
        5

    means the code **2 + a** should print the value 5



IC50 reader
-------------------

Before starting, we first need to get a data file to play with. 
For now, we do not need to enter into the details of the expected data
structure; it should be a CSV or TSV file as in this :download:`IC50 example <../../share/data/test_IC50.csv>` file.

.. seealso:: more details about the data format can be found in the :ref:`data` section as well as links to retrieve IC50 data sets.


Although all functionalities could be imported using::

    from gdsctools import *

we will be try to be as explicit as possible in the following examples. So we will rather use (for instance)::

    from gdsctools import IC50

This is better coding practice and has also the advantage of telling beginners
which functions are going to be used. 

Here above, we imported the :class:`~gdsctools.readers.IC50` class, that allows one to read a data set such as the downloadable file at the top of this page. We will explain in details the different data sets and their formats in the :ref:`data` section. However, for now it is enough to know that it should be a CSV formatted file that contains IC50s; one value for each combination of drug and cell line. 

Note that the location 
of that file can be easily obtained using::

    from gdsctools import ic50_test
    print(ic50_test.filename)

This structure **ic50_test** does not contain any data per
se but the location of the file itself, which can then be read with a 
dedicated :class:`~gdsctools.readers.IC50` class::

    >>> from gdsctools import IC50, ic50_test
    >>> ic = IC50(ic50_test)
    >> print(ic)
    Number of drugs: 11
    Number of cell lines: 988
    Percentage of NA 0.206569746043

As you can see you can get some information about the IC50 content (e.g., 
number of drugs, percentage of NaNs). See :class:`gdsctools.readers.IC50` for more details.

The ANOVA class
----------------
Given an IC50 data set, we can now analyse it using the main class 
called :class:`~gdsctools.anova.ANOVA`. A default set of 680 genomic features 
is provided and we do not need to worry about it right now.

Before starting, just a few words about the underlying stastistical analysis. In a given analysis, there are :math:`N_d` drugs and :math:`N_c` cell lines. Each combination of drug and cell line has a measured IC50. A set of genomic features is provided and the corresponding :math:`$N_c$` cell lines used to get :math:`N_f` features. Then, for each drug, we compute the association (a regression analysis) between a drug and a feature leading to a p-value. This calculation is possibly repeated across all features and even all drugs. Consequently, a multiple testing correction is applied and reported in the analysis. For more information, please see the :ref:`userguide` section.

Let us now create a structure for the ANOVA analysis (the only input is the IC50)::

    from gdsctools import ANOVA, ic50_test
    gdsc = ANOVA(ic50_test)

As you can see here, we just provide the ic50_test name but one can provide the filename of a TSV or CSV file. Actually, you could provide a variety of input (dataframes, IC50 instance). The following statements are equivalent::

    from gdsctools import ANOVA, ic50_test, IC50
    gdsc = ANOVA(ic50_test)
    gdsc = ANOVA(ic50_test.filename)
    gdsc = ANOVA(IC50(ic50_test))
    gdsc = ANOVA("localfile.csv")

As briefly mentionned earlier, you can perform 3 types of analysis:

.. index:: ODOF, ODAF, ADAF

#. compute one association between a drug and feature (ODOF)
#. compute the associations between one drug and all the features (ODAF)
#. compute all associations for all drugs and all features. (ADAF)


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



