.. _html:

HTML report
==============

The output of an ANOVA analysis can be used to create an HTML report.

First, let us generate the results again (see :ref:`quickstart`).

::  
    
    >>> from gdsctools import ANOVA, ic50_test
    >>> gdsc = ANOVA(ic50_test)
    >>> gdsc.set_cancer_type('breast')
    >>> results = gdsc.anova_all() 

The :attr:`results` variable contains a dataframe :attr:`df`. This dataframe 
contains as many rows as associations of
drugs and features. See :class:`~gdsctools.anova_results.ANOVAResults` for the contents. The HTML report extracts the significant associations, and then create figures and HTML pages for each of the associations that are significant.  You can easily create HTML report as follows::

    >>> from gdsctools import ANOVAReport
    >>> report = ANOVAReport(gdsc, results)
    >>> report.create_html_pages()

The report creates a **Data Package**, which details can be found in the :ref:`data_packages` section.

Images are created for each significant associations and may take a while.

Some tunable settings are available in the :attr:`settings` (see :ref:`settings`). For instance, you can set the output directory to a user value instead of (html_gdsc_anova)::

    >>> report.settings.directory = 'BLCA'





