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
"""Code related to the Regression analysis to find associations between drug IC50s
and genomic features"""
import os
import glob
import json

from gdsctools.report import HTMLTable, ReportMain
from gdsctools.readers import DrugDecode
from gdsctools.volcano import ScatterJS

import colorlog as logger

import pandas as pd
import pylab

import easydev
from colormap import cmap_builder


__all__ = ['RegressionReport']


class RegressionReport(object):
    """Class used to interpret the results and create final HTML report


    """
    def __init__(self, method, directory=".", verbose=True, image_dir="images",
            data_dir="data", config={"boxplot_n": "?"}):
        """.. rubric:: Constructor

        :param method: Method used in the regression analysis (lasso,
            elasticnet, ridge)
        :param results: 

        """
        self.image_dir = image_dir
        self.method = method
        self.config = config

        self.directory = directory
        self.output_dir = "."
        self.verbose = verbose
        self.prefix = "gdsctools_regression_"
        self.prefix_images = image_dir + os.sep +"gdsctools_regression_"
        self.prefix_data = data_dir + os.sep +"gdsctools_regression_"
        self.filenames = glob.glob(self.prefix_images + "boxplot_*png")
        self.drugids =  [this.rstrip(".png").lstrip(self.prefix_images).lstrip("boxplot_") 
                        for this in self.filenames]

    def create_html_drug(self):
        """report for each individual drug"""
        for drugid in self.drugids:
            logger.info('Creating HTML report for drug %s' % drugid)
            report = HTMLOneDrug(drugid, caller=self)
            report.create_report(onweb=False)

    def create_html_main(self, onweb=False):
        """Create HTML main document (summary)"""

        if self.verbose:
            logger.info("Creating main HTML page in directory %s" %
                (self.directory))
        ReportMain(directory=self.directory, verbose=self.verbose, mode="summary")

        html = HTMLPageMain(caller=self)
        html.create_report(onweb=onweb)


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

        filename_template = self.caller.prefix_data + "%(name)s_" + "%s." % self.drug 
        results_filename = filename_template % {"name":"results"} + "json"

        with open(results_filename, "r") as fh:
            data = json.loads(fh.read())
        try:data["bayes"] = easydev.precision(data['bayes'], 3)
        except:pass
        try:data["alpha"] = easydev.precision(data['alpha'], 5)
        except:pass
        try:data["Rp"] = easydev.precision(data['Rp'], 4)
        except:pass
        self.params.update(data)
        self.params['method'] = self.caller.method
        self.jinja['sections'] = []
        self.jinja['goback'] = True

    def _create_report(self, onweb=True):
        section = """<div>
        <b>DrugID:</b> %(drugid)s</br>
        <b>Regression method:</b> %(method)s </br>
        <b>Regression, alpha parameter used:</b> %(alpha)s</br>
        <b>Bayes factor:</b> %(bayes)s</br>
        <b>Coefficient of regression (pearson):</b> %(Rp)s</br>
        </div>
        """ % self.params
        self.jinja['sections'].append(section)

        text = {}
        text['boxplot'] = ("This boxplot shows the %s most important features "
            "(based on the weights of the regression).")
        text['boxplot'] %= self.caller.config["boxplot_n"]

        text['importance'] = ("Feature with non-null weights. If empty, it"
            " means no feature of interests were found")

        text['randomness'] = ("Here we run the regression analysis %s times and "
            "plot the regression value (x-axis) for the real data (blue) "
            "and randomising the variable to explain (red). " )
        text['randomness'] %= self.caller.config['randomness']

        text['weights'] = ("Feature with non-null weights. If empty, it"
            " means no feature of interests were found")

        for this in ["boxplot", "randomness", "importance", "weights"]:
            self.params['name'] = this
            self.params['text'] = text[this]
            filename = self.caller.prefix_images + "%(name)s_%(drugid)s.png" % self.params
            self.params["filename"] = filename
            self.params['title'] = this.title()
            if os.path.exists(filename):
                section = """<div>
                <h2>%(title)s results</h2>
                <p>%(text)s</p>
                <img src="%(filename)s">
                """ % self.params
                self.jinja['sections'].append(section)
            else:
                logger.warning("%s not found. Skipped" % filename)


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

        # The top section with standard information
        section = """<div>
        <b>Regression method:</b> %s </br>
        </div><hr>
        """ % self.caller.method
        self.jinja['sections'].append(section)


        # The main CSV tables with bayes factor and links to each drug ID
        filename = self.caller.prefix_data + "results.csv"

        df = pd.read_csv(filename)
        df['ttest (-log10)'] = -pylab.log10(df['ttest'])
        # prevents inf to fail in the HTMLTable
        table = HTMLTable(df)
        table.add_bgcolor('bayes')
        table.add_bgcolor('Rp')
        table.df['drugid'] = ['<a href="drug_%s.html">%s</a>' % (x, x)
                          for x in table.df['drugid']]

        html = ("<div><p>This table contains links to all drugs (first column)."
                " The Rp column contains the coefficient of correlation"
                " (pearson) found with the regression method for the alpha"
                " parameter provided in column 3. The alpha value is the optimised"
                " value obtained using a cross validation (see below)."
                " The ln_alpha column is just the -log10(alpha) value. The bayes"
                " factor gives an idea of the significance of the correlation as "
                " compared to a null distribution. See"
                ' <a href="http://gdsctools.readthedocs.io/en/master/references.html">'
                'gdstools documentation.</a> for details.'
                "<br>"
                " Note also that the optimisation of the alpha parameter is"
                " performed using a cross validation and depends on a few"
                " parameters such as the range of alpha values, number of "
                " cross validation, ....</p>")
        html += table.to_html(index=False)

        pattern = '<div>%s <p>Download the CSV <a href="%s">file</a></p></div><hr>' 
        pattern = pattern % (html, filename)
        html = pattern

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
        filename = self.caller.prefix_data + "results.csv"
        df = pd.read_csv(filename)
        df["markersize"] = 20
        js = ScatterJS(df, x="Rp", y="bayes", color="bayes", size="markersize")
        js.xlabel = "Coefficient correlation (pearson)"
        js.ylabel = "Bayes factor"
        self.jinja['volcano_jsdata'] = js.get_html()

