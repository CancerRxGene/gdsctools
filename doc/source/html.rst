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
drugs and features. See :class:`~gdsctools.anova.ANOVAResults` for the contents.  Mosf of them being non significant. The HTML report extracts the significant associations, and then create figures and HTML pages for each of them.  You can easily create HTML report as follows::

    >>> report = ANOVAReport(gdsc, results)
    >>> report.create_html_pages()

To know more about the content of the HTML report, please see the
:ref:`data_packages` section.


Some tunable settings are available in the :attr:`settings` (see :ref:`settings`). For instance, you can set the output directory::

    >>> report.settings.directory = 'BLCA'





