# coding=utf-8
# -*- python -*-
#
#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.com>
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
"""Base classes to create HTML reports easily"""
import os
import shutil
import glob
sepjoin = os.sep.join

from gdsctools import gdsctools_data

import easydev
import pandas as pd
from jinja2.environment import Environment
from jinja2 import FileSystemLoader

from colormap import rgb2hex, cmap_builder
# note that the sorttable javascript is from
# `http://www.kryogenix.org/code/browser/sorttable/
# with an X11 license


from reports import HTMLTable

__all__ = ['ReportMain']


# Should re-user reports package
from reports import Report


class ReportMain(Report):
    """A base class to create HTML pages

    This :class:`Report` class holds the CSS and HTML layout and will ease
    the creation of new reports and HTML pages. For instance, it will add
    a footer and header automatically, save files in the proper directory,
    create the directory if it is missing, copy CSS and JS files in the
    directory automatically.

    ::

        from gdsctools import Report
        r = Report()
        r.add_section('Example with some text', 'Example' )
        r.report(onweb=True)

    .. note:: **For developers** the original CSS and JS files are stored in
        the share/data directory.

    The idea is that you create sections (text +  title) that you add little by
    little in your HTML documents. Then, you create the report. The report will
    add a footer, header, table of contents before the sections. The
    **text** of a section can contain any HTML document.

    """

    def __init__(self, filename='index.html', directory='report',
                 overwrite=True, verbose=True, template_filename='index.html', 
                 mode=None, init_report=True):
        """.. rubric:: Constructor

        :param filename: default to **index.html**
        :param directory: defaults to **report**
        :param overwrite: default to True
        :param verbose: default to True
        :param dependencies: add the dependencies table at the end of the
            document if True.
        :param str mode: if none, report have a structure that contains the
            following directories: OUTPUT, INPUT, js, css, images, code.
            Otherwise, if mode is set to 'summary', only the following
            directories are created: js, css, images 

        """
        gdsctools_path = easydev.get_package_location('gdsctools')
        extra_css_path = sepjoin([gdsctools_path, "gdsctools", "data", "css"])
        extra_js_path = sepjoin([gdsctools_path, "gdsctools", "data", "javascript"])

        extra_css_list = glob.glob(extra_css_path + os.sep + "*css")
        extra_js_list = glob.glob(extra_js_path + os.sep + "*js")

        searchpath = sepjoin([gdsctools_path, "gdsctools", "data", "templates"])

        super(ReportMain, self).__init__(searchpath, filename=filename,
                template_filename=template_filename, directory=directory,
                extra_css_list=extra_css_list, extra_js_list=extra_js_list,
                init_report=init_report)

        self.jinja['dependencies'] = self.get_table_dependencies("gdsctools").to_html()
        self.jinja['analysis'] = 'anova'
        from gdsctools import version
        self.jinja['version'] = version
        self.jinja['title'] = 'ANOVA analysis summary'
        self.jinja["analysis_domain"] = "PANCAN"
        self.jinja['resource_path'] = "."

        self._directory = directory
        self._filename = filename

        if mode is None:
            self._to_create = ['OUTPUT', 'INPUT', 'images', 'css',
                    'js', 'code', 'associations']
        elif mode == 'summary':
            self._to_create = ['images', 'css', 'js',]
        if init_report:
            self._init_report()

    def show(self):
        """Opens a tab in a browser to see the document"""
        from easydev.browser import browse as bs
        bs(self.abspath)

    def _init_report(self):
        super(ReportMain, self)._init_report()

        for filename in ['EBI_logo.png', 'sanger-logo.png']:
            target = os.sep.join([self.directory, 'images', filename ])
            if os.path.isfile(target) is False:
                filename = gdsctools_data("images" + os.sep + filename)
                shutil.copy(filename, target)

    def create_report(self, onweb=True):
        self._create_report()
        self.write()
        if onweb is True:
            self.onweb()

