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

__all__ = ['ReportMAIN']


# Should re-user reports package
from reports import Report


class ReportMAIN(Report):
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
                 overwrite=True, verbose=True,
                template_filename='index.html', mode=None):
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
            directories are created: js, css, images, code

        """
        #: name of the analysis added in the title
        self.analysis = 'anova'
        self.pkgname = 'gdsctools'

        from gdsctools import version
        #: version added in the sub title
        self.version = version

        self._directory = directory
        self._filename = filename

        # This contains the sections and their names when
        # method add_section is used
        self.sections = []
        self.section_names = []

        #: flag to add dependencies
        self.add_dependencies = False

        self.title = 'ANOVA analysis summary'
        self.analysis_type = "PANCAN"

        # For jinja2 inheritance, we need to use the environment
        # to indicate where are the parents' templates
        template_directory = gdsctools_data('templates')

        self.env = Environment()
        self.env.loader = FileSystemLoader(template_directory)

        # use template provided inside gdsctools
        self.template = self.env.get_template(template_filename)

        self.jinja = {
                'time_now': self.get_time_now(),
                "analysis": self.analysis,
                "version": self.version,
                "title": self.title,
                "analysis_domain": self.analysis_type,
                'dependencies': self.get_table_dependencies().to_html(),
                }

        if mode is None:
            self._to_create = ['OUTPUT', 'INPUT', 'images', 'css',
                    'js', 'code']
        elif mode == 'summary':
            self._to_create = ['images', 'css', 'js',]

        self._init_report()

    def _get_filename(self):
        return self._filename
    def _set_filename(self, filename):
        self._filename = filename
    filename = property(_get_filename, _set_filename,
        doc="The filename of the HTML document")

    def _get_directory(self):
        return self._directory
    def _set_directory(self, directory):
        self._directory = directory
    directory = property(_get_directory, _set_directory,
            doc="The directory where to save the HTML document")

    def _get_abspath(self):
        return self.directory + os.sep + self.filename
    abspath = property(_get_abspath,
            doc="The absolute path of the document (read only)")

    def show(self):
        """Opens a tab in a browser to see the document"""
        from easydev.browser import browse as bs
        bs(self.abspath)

    def _init_report(self):
        """create the report directory and return the directory name"""
        self.sections = []
        self.section_names = []
        # if the directory already exists, print a warning

        try:
            if os.path.isdir(self.directory) is False:
                print("Created directory {}".format(self.directory))
                os.mkdir(self.directory)

            # list of directories created in the constructor
            for this in self._to_create:
                try:
                    os.mkdir(self.directory + os.sep + this)
                except:
                    pass # already created ?
        except Exception:
            pass
        finally:
            for filename in ['gdsc.css', 'github-gist.css']:
                target = os.sep.join([self.directory, 'css', filename ])
                if os.path.isfile(target) is False:
                    filename = gdsctools_data(filename)
                    shutil.copy(filename, target)

            for filename in ['sorttable.js', 'highlight.pack.js']:
                target = os.sep.join([self.directory, 'js', filename ])
                if os.path.isfile(target) is False:
                    filename = gdsctools_data(filename)
                    shutil.copy(filename, target)

            for filename in ['EBI_logo.png', 'sanger-logo.png']:
                target = os.sep.join([self.directory, 'images', filename ])
                if os.path.isfile(target) is False:
                    filename = gdsctools_data("images" + os.sep + filename)
                    shutil.copy(filename, target)


    def to_html(self):
        self.jinja['time_now'] = self.get_time_now()
        return self.template.render(self.jinja)

    def write(self):
        with open(self.abspath, "w") as fh:
            data = self.to_html()
            fh.write(data)

    def onweb(self):
        """Open the HTML document in a browser"""
        from easydev import onweb
        onweb(self.abspath)

    def create_report(self, onweb=True):
        try:
            # some parent zill have that method implemented
            self._create_report()
        except:
            pass
        self.write()
        if onweb is True:
            self.onweb()

    def get_time_now(self):
        """Returns a time stamp"""
        import datetime
        import getpass
        username = getpass.getuser()
        # this is not working on some systems: os.environ["USERNAME"]
        timenow = str(datetime.datetime.now())
        timenow = timenow.split('.')[0]
        msg = '<div class="date">Created on ' + timenow
        msg += " by " + username +'</div>'
        return msg

    def get_table_dependencies(self):
        """Returns dependencies of the pipeline as an HTML/XML table

        The dependencies are the python dependencies as returned by
        pkg_resource module.

        """
        dependencies = easydev.get_dependencies(self.pkgname)
        # TODO: Could re-use new method in HTMLTable for adding href
        # but needs some extra work in the add_href method.
        names = [x.project_name for x in dependencies]
        versions = [x.version for x in dependencies]
        links = ["""https://pypi.python.org/pypi/%s""" % p for p in names]
        df = pd.DataFrame({
            'package': ["""<a href="%s">%s</a>""" % (links[i], p)
                for i, p in enumerate(names)],
            'version': versions})
        table = HTMLTable(df, name="dependencies", escape=False)
        table.sort('package')
        return table
