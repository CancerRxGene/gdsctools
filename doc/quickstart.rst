

.. _quickstart:

Quick Start
=============

ANOVA analysis
-----------------

If you already know what you can do with **GDSCTools**, we assume you have a well
formatted :term:`IC50` matrix and a genomic features binary matrix. Then, 
you can run the entire ANOVA analysis as follows::


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


The **results** variable contains all tested associations within a single 
dataframe. The report will focus on significant associations and create boxplots or volcano plots accordingly.


More details about the ANOVA analysis itself can be found in the
:ref:`anova_partone` and :ref:`anova_parttwo` sections. The data structure can
be found in the next section (:ref:`data`).


Regression analysis
----------------------

Similarly, for the regression analysis, one can write a script as above::

    from gdsctools import GDSCLasso
    gdsc = GDSCLasso(IC50_filename, genomic_feature_filename)
    for drugID in lasso.drugIds:
        res = lasso.runCV(drugid, kfolds=8)
        best_model = lasso.get_model(alpha=res.alpha)
        weights = lasso.plot_weight(drugid, best_model)
        boxplots = lasso.boxplot(drugid, model=best_model, n=10, bx_vert=False)

However, we could recommend to use a worflow designed for this analysis. In a
shell, type::

    gdsctools_regression -I IC50_filename -F genomic_feature_filename 
        --method lasso -o analysis
    cd analysis

Edit the config.yaml file to change any parameters. Then, execute the pipeline::

    snakemake -s regression.rules

See :ref:`multivariate_regression` section for details.




