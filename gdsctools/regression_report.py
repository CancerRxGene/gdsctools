# -*- python -*-
# -*- coding utf-8 -*-
#
#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.comWE HERE>
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
"""Code related to the ANOVA analysis to find associations between drug IC50s
and genomic features"""
import os
import glob
import json

from gdsctools.report import HTMLTable, ReportMain
from gdsctools.readers import DrugDecode
from gdsctools.volcano import ScatterJS

import pandas as pd
import pylab

import easydev
from colormap import cmap_builder


__all__ = ['RegressionReport']


class RegressionReport(object):
    """Class used to interpret the results and create final HTML report


    """
    def __init__(self, method, directory=".", verbose=True):
        """.. rubric:: Constructor

        :param gdsc: the instance with which you created the results to report
        :param results: the results returned by :meth:`ANOVA.anova_all`. If
            not provided, the ANOVA is run on the fly.

        """
        self.method = method
        self.directory = directory
        self.output_dir = "."
        self.verbose = verbose
        self.prefix = "gdsctools_regression_"
        self.filenames = glob.glob(self.prefix + "boxplot_*png")
        self.drugids =  [this.rstrip(".png").lstrip(self.prefix).lstrip("boxplot_") 
                        for this in self.filenames]

    def create_html_drug(self):
        for drugid in self.drugids:
            print('Creating HTML report for drug %s' % drugid)
            report = HTMLOneDrug(drugid, caller=self)
            report.create_report(onweb=False)

    def create_html_main(self, onweb=False):
        """Create HTML main document (summary)"""

        if self.verbose:
            print("Creating main HTML page in directory %s" %
                (self.directory))
        ReportMain(directory=self.directory, verbose=self.verbose, mode="summary")

        html = HTMLPageMain(caller=self)
        html.create_report(onweb=onweb)

        # report for each individual drug
        self.create_html_drug()


class HTMLOneDrug(ReportMain):
    def __init__(self, drugid, caller):
        self.drug = int(drugid)
        self.caller = caller

        filename = "drug_{0}.html".format(self.drug)
        super(HTMLOneDrug, self).__init__(
                directory=caller.output_dir,
                filename=filename, 
                template_filename='regression.html',
                init_report=False)
        self.title = 'Single Drug analysis (%s)' % self.drug
        self.params = {"drugid": self.drug}

        self.filename_template = self.caller.prefix + "%(name)s_" + "%s." % self.drug 
        results_filename = self.filename_template % {"name":"results"} + "json"

        with open(results_filename, "r") as fh:
            data = json.loads(fh.read())
        try:data["bayes"] = easydev.precision(data['bayes'], 3)
        except:pass
        try:data["alpha"] = easydev.precision(data['alpha'], 3)
        except:pass
        self.params.update(data)
        self.params['method'] = self.caller.method
        self.jinja['sections'] = []

    def _create_report(self, onweb=True):
        section = """<div>
        <b>DrugID:</b> %(drugid)s</br>
        <b>Regression method:</b> %(method)s </br>
        <b>Regression, alpha parameter used:</b> %(alpha)s</br>
        <b>Bayes factor:</b> %(bayes)s</br>
        <b>Regression coeffcient:</b> %(Rp)s</br>
        </div>
        """ % self.params
        self.jinja['sections'].append(section)

        text = {}
        text['boxplot'] = "to do"
        text['importance'] = ("Feature with non-null weights. If empty, it"
            " means no feature of interests were found")
        text['randomness'] = """Here we run the regression analysis N times and
plot the regression value (x-axis) for the real data (blue) and randomising the
variable to explain (red). The bayes and ttest metric are then computed. """
        text['weights'] = ("Feature with non-null weights. If empty, it"
            " means no feature of interests were found")

        for this in ["boxplot", "randomness", "importance", "weights"]:
            self.params['name'] = this
            self.params['text'] = text[this]
            filename = self.caller.prefix + "%(name)s_%(drugid)s.png" % self.params
            self.params["filename"] = filename
            if os.path.exists(filename):
                section = """<div>
                <h2>%(name)s results</h2>
                <p>%(text)s</p>
                <img src="%(filename)s">
                """ % self.params
                self.jinja['sections'].append(section)
            else:
                print("%s not found. Skipped" % filename)


class HTMLPageMain(ReportMain):
    def __init__(self, caller):
        sepjoin = os.sep.join
        super(HTMLPageMain, self).__init__(
                directory=caller.directory,
                template_filename="regression.html",
                filename="index.html", mode="summary")

        self.caller = caller
        self.jinja['analysis_domain'] = "PANCAN"
        self.jinja['title'] = "Regression analysis summary"
        self.jinja['sections'] = []
        #self.jinja["collaborator"] = report.company

    def _create_report(self, onweb=True):
        section = """<div>
        <b>Regression method:</b> %s </br>
        </div><hr>
        """ % self.caller.method
        self.jinja['sections'].append(section)

        filename = self.caller.prefix + "results.csv"
        df = pd.read_csv(filename)
        df['ttest (-log10)'] = -pylab.log10(df['ttest'])
        # prevents inf to fail in the HTMLTable
        table = HTMLTable(df)
        table.add_bgcolor('bayes')
        table.add_bgcolor('Rp')
        table.df['drugid'] = ['<a href="drug_%s.html">%s</a>' % (x,x)
                          for x in table.df['drugid']]

        html = ("<div><p>This table contains links to all drugs. The Rp columns"
                " contains the coefficient of correlation found with the"
                " method for the give alpha parameter. The ln_alpha column"
                " is just the -log10(alpha) value. The bayes and ttest columns"
                " gives an idea of the significance of the correlation as "
                " compared to a null distribution.</p>")
        html += table.to_html(index=False) +"</div>"
        self.jinja['sections'].append(html)

        filename = self.caller.prefix + "scatter_plot.png"
        html = "<hr><div>"
        html += "<img src=%s></img>" % filename
        html += "/<div>"

        self.jinja['sections'].append(html)

        # The scatter plot. First the javascript in the header
        self._set_scatter()
        # and the section itself
        html = """
         <div class="wrap">
          <div class="content">
              <center>
                  <canvas id='canvasVolcano' width='800' height='540'></canvas>
              </center>
          </div>
          <div class="clear">&nbsp;</div>
        </div>
        """
        self.jinja["sections"].append(html)

    def _set_scatter(self):
        filename = self.caller.prefix + "results.csv"
        df = pd.read_csv(filename)
        df["log_ttest"] = -pylab.log10(df["ttest"])
        df["markersize"] = 20
        js = ScatterJS(df, x="Rp", y="log_ttest", color="bayes", size="markersize" )
        js.xlabel = "Regression coefficient"
        js.ylabel = "log10(ttest pvalue)"
        self.jinja['volcano_jsdata'] = js.get_html()
        

