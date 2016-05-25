from gdsctools import *
from easydev import assert_list_almost_equal
import numpy as np
import os


def test_skcm():

    an = ANOVA(gdsctools_data("test_ANOVA_input_v17_skcm.txt"), 
                gdsctools_data("test_ANOVA_input_v17_skcm.txt"))
    an.settings.pvalue_correction_method = 'qvalue'
    results = an.anova_all()

    # This create a temp directory called "skin"
    report = ANOVAReport(an, results)

    diag = report.diagnostics()

    diag = diag.set_index('text')

    assert diag.ix["Total number of ANOVA tests performed"].values == 5194
    assert diag.ix["Percentage of tests performed"].values == 81.66
    assert diag.ix["Total number of tested drugs"].values ==  265
    assert diag.ix["Total number of genomic features used"].values == 24
    assert diag.ix["Total number of screened cell lines"].values == 55
    assert diag.ix["MicroSatellite instability included as factor"].values == False
    assert diag.ix["Total number of significant associations"].values == 13
    assert diag.ix[" - sensitive"].values == 8
    assert diag.ix[" - resistant"].values == 5
    assert diag.ix["p-value significance threshold"].values == 0.001
    assert diag.ix["FDR significance threshold"].values == 25
    assert diag.ix["Range of significant p-values"].values[0] == "[9.87e-08, 0.0006358]"
    assert diag.ix["Range of significant % FDRs"].values[0] == "[0.04777 23.67]"

    assert_list_almost_equal(report.df.ix[0].values,
        np.array([1, 'BRAF_mut', 1373, None, None, 10, 35,
            -2.1750318847152079, 2.9802302648104275, -5.155262149525635,
               2.291545942648078, 3.2964327113036669, 2.1492596576572947,
               1.5638912124151449, 2.3986223028747027, 9.8695117331183039e-08,
               9.8695117331182668e-08, None, None, None, 0.047768436788292422],
            dtype=object))

    assert_list_almost_equal(
        [report.df.ix[12]['ANOVA_FEATURE_FDR']], [23.671582956786185])

    report.create_html_pages(onweb=False)

    assert os.path.exists("skin/index.html")
    assert os.path.exists("skin/associations/manova.html")
    assert os.path.exists("skin/associations/a1.html")
    assert os.path.exists("skin/associations/BRAF_mut.html")
    assert os.path.exists("skin/associations/drug_1047.html")




