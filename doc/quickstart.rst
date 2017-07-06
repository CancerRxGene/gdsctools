

.. _quickstart:

Quick Start
=============

Data
----------

In the following examples and most of the examples in GDSCTools, we use 
data sets known as version 5 (v5) and version 17 (v17). Each set is made of  
a drug responses file (IC50s) and a genomic features file. The 4 files
are copies of original files to be found on the `CancerRxGene
<http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html>`_ page. 
Altough not required for the following examples, you may want to 
download them directly as follows::

    wget https://tinyurl.com/y7nn6e5h -O genomic_features_v17.csv.gz
    wget https://tinyurl.com/ybtlsrpz -O genomic_features_v5.csv.gz
    wget https://tinyurl.com/ycavjd37 -O IC50_v17.csv.gz
    wget https://tinyurl.com/yakfnqmb -O IC50_v5.csv.gz


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

And create an HTML report as follows::

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

Similarly, for the regression analysis, one can write a script as above. Here, 
we restrict the analysis to the first 4 drugs (figures are open for each
drug)::

    from gdsctools import GDSCLasso
    lasso = GDSCLasso(IC50_filename, genomic_feature_filename)
    for drugid in lasso.drugIds[0:4]:
        res = lasso.runCV(drugid, kfolds=8)
        best_model = lasso.get_model(alpha=res.alpha)
        #weights = lasso.plot_weight(drugid, best_model)
        boxplots = lasso.boxplot(drugid, model=best_model, n=10, bx_vert=False)

However, we would recommend to use a worflow designed for this analysis. If you
type this command in a shell::

    gdsctools_regression -I IC50_filename -F genomic_feature_filename 
        --method lasso -o analysis
    cd analysis

it creates a Snakemake pipeline and its configuration file. You can then edit
file named **config.yaml** and once done, execute the pipeline::

    snakemake -s regression.rules

See :ref:`multivariate_regression` section for details.




