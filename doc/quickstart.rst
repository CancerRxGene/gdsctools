

.. _quickstart:

Quick Start
=============

You know what you are doing with **GDSCTools** already, so you have a well
formatted :term:`IC50` matrix and a genomic features binary matrix, then you can run the entire ANOVA analysis as follows::


    from gdsctools import ANOVA
    # For example, use these test files
    # from gdsctools import ic50_test as ic50_filename
    # from gdsctools import gf_v17 as genomic_feature_filename
    gdsc = ANOVA(IC50_filename, genomic_feature_filename)
    results = gdsc.anova_all()

And create an HTML report as follows:::

    from gdsctools import ANOVAReport
    report = ANOVAReport(gdsc, results) 
    report.create_html_pages()


The results variable contains all tested associations within a single 
dataframe. The report will focus on significant associations and create boxplots or volcano plots accordingly.


More details about the ANOVA analysis itself can be found in the
:ref:`anova_partone` and :ref:`anova_parttwo` sections. The data structure can
be found in the next section (:ref:`data`).


