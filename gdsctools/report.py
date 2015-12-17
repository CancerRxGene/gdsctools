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

import easydev
import pandas as pd
from jinja2.environment import Environment
from jinja2 import FileSystemLoader

from colormap import rgb2hex, cmap_builder
# note that the sorttable javascript is from
# `http://www.kryogenix.org/code/browser/sorttable/
# with an X11 license

__all__ = ['HTMLTable', 'ReportMAIN']


class HTMLTable(object):
    """Handler to export dataframe into HTML table.

    Dataframe in Pandas already have a to_html method to export the dataframe
    into a HTML formatted table. However, we provide here a few handy features:

        * Takes each cell in a given column and creates an HTML
          reference in each cell. See :meth:`add_href` method.
        * add an HTML background into cells (numeric content) of
          a given column using different methods (e.g., normalise).
          See :meth:`add_bgcolor`

    ::

        import pandas as pd
        df = pd.DataFrame({'A':[1,2,10], 'B':[1,10,2]})
        from gdsctools import HTMLTable
        html = HTMLTable(df)

    .. note:: similar project exists such as prettytable but could not do
        exactly what we wanted at the time gdsctools was developed.

    .. note:: Could be moved to biokit or easydev package.

    """
    def __init__(self, df, name=None, **kargs):
        """.. rubric:: Constructor


        :param dataframe df: a pandas dataframe to transform into a table
        :param str name: not used yet

        There is an :attr:`pd_options` attribute to reduce the max column
        width or the precision of the numerical values.

        """
        self.df = df.copy() # because we will change its contents possibly
        self.name = name
        self.pd_options = {
                'max_colwidth': -1,
                'precision': 2}

    def to_html(self, index=False, escape=False, header=True,
            collapse_table=True, **kargs):
        """Return HTML version of the table

        This is a wrapper of the to_html method of the pandas dataframe.

        :param bool index: do not include the index
        :param bool escape: do not escape special characters
        :param bool header: include header
        :param bool collapse_table: long tables are shorten with a scroll bar
        :param kargs: any parameter accepted by
            :meth:`pandas.DataFrame.to_html`

        """
        _buffer = {}
        for k, v in self.pd_options.items():
            # save the current option
            _buffer[k] = pd.get_option(k)
            # set with user value
            pd.set_option(k, v)

        # class sortable is to use the sorttable javascript
        # note that the class has one t and the javascript library has 2
        # as in the original version of sorttable.js
        table = self.df.to_html(escape=escape, header=header, index=index,
                classes='sortable', **kargs)

        # get back to default options
        for k, v in _buffer.items():
            pd.set_option(k, v)
        if len(self.df) > 20 and collapse_table is True:
            return '<div class="table_outer">' + table+"</div>"
        else:
            return '<div class="table">' + table+"</div>"

    def add_bgcolor(self, colname, cmap='copper', mode='absmax',
            threshold=2):
        """Change column content into HTML paragraph with background color

        :param colname:
        :param cmap: a colormap (matplotlib) or created using
            colormap package (from pypi).
        :param mode: type of normalisation in 'absmax', 'max', 'clip'
            (see details below)
        :param threshold: used if mode is set to 'clip'

        Colormap have values between 0 and 1 so we need to normalised the data
        between 0 and 1. There are 3 mode to normalise the data so far.

        If mode is set to 'absmax', negatives and positives values are
        expected to be found in a range from -inf to inf. Values are
        scaled in between [0,1] X' = (X / M +1) /2. where m is the absolute
        maximum. Ideally a colormap should be made of 3 colors, the first
        color used for negative values, the second for zeros and third color
        for positive values.

        If mode is set to 'clip', values are clipped to a max value (parameter
        *threshold* and values are normalised by that same threshold.

        If mode is set to 'max', values are normalised by the max.

        """
        try:
            # if no cmap provided, it may be just a known cmap name
            cmap = cmap_builder(cmap)
        except:
            pass

        data = self.df[colname].values

        if len(data) == 0:
            return
        if mode == 'clip':
            data = [min(x, threshold)/float(threshold) for x in data]
        elif mode == 'absmax':
            m = abs(data.min())
            M = abs(data.max())
            M = max([m, M])
            data = (data / M + 1)/2.
        elif mode == 'max':
            data = data / float(data.max())

        # the expected RGB values for a given data point
        rgbcolors = [cmap(x)[0:3] for x in data]
        hexcolors = [rgb2hex(*x, normalised=True) for x in rgbcolors]

        # need to read original data again
        data = self.df[colname].values
        # need to set precision since this is going to be a text not a number
        # so pandas will not use the precision for those cases:
        data = [easydev.precision(x, self.pd_options['precision'])
                for x in data]
        html_formatter = '<p style="background-color:{0}">{1}</p>'
        self.df[colname] = [html_formatter.format(x, y)
                for x, y in zip(hexcolors, data)]

    def add_href(self, colname, url=None, newtab=False, suffix=None):
        """

        default behaviour: takes column content and put into::

            <a href={content}.html>content</a>

        This is used to link to local files. If url is provided, you typically
        want to link to an external url where the content is an identifier::

            <a href={url}{content}>content</a>

        Note that in the first case, *.html* is appended but not in the second
        case, which means cell's content should already have the .html
        Also in the second case, a new tab is open whereas in the first case
        the url is open in the current tab.

        .. note:: this api may change in the future.

        """
        if url is not None:
            if suffix is None:
                suffix = ''
            if newtab is False:
                formatter = '<a  alt="{1}" href="{0}{1}{2}">{1}</a>'
            else:
                formatter = '<a target="_blank" alt={1} href="{0}{1}{2}">{1}</a>'
            self.df[colname] = self.df[colname].apply(lambda x:
                    formatter.format(url, x, suffix))
        else:
            if suffix is None:
                suffix = '.html'

            if newtab is False:
                formatter = '<a alt="{1}" href="{0}{2}">{1}</a>'
            else:
                formatter = '<a target="_blank" alt="{1}" href="{0}{2}">{1}</a>'
            self.df[colname] = self.df[colname].apply(lambda x:
                formatter.format(x,x, suffix))

    def sort(self, name, ascending=True):
        # for different pandas implementations
        try:
            self.df.sort_values(by=name, inplace=True, ascending=ascending)
        except:
            self.df.sort(columns=name, inplace=True, ascending=ascending)


class ReportMAIN(object):
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
        gdsctools_path = easydev.get_shared_directory_path('gdsctools')
        template_directory = os.sep.join([gdsctools_path, 'data', 'templates'])

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
                    filename = easydev.get_share_file("gdsctools", "data",
                        filename)
                    shutil.copy(filename, target)
            
            for filename in ['sorttable.js', 'highlight.pack.js']:
                target = os.sep.join([self.directory, 'js', filename ])
                if os.path.isfile(target) is False:
                    filename = easydev.get_share_file("gdsctools", "data",
                        filename)
                    shutil.copy(filename, target)

            for filename in ['EBI_logo.png', 'sanger-logo.png']:
                target = os.sep.join([self.directory, 'images', filename ])
                if os.path.isfile(target) is False:
                    dire = 'data' + os.sep + 'images'
                    filename = easydev.get_share_file("gdsctools", dire,
                        filename)
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

    def _create_report(self):
        pass
        # this shoudl populate the self.params and create the figures.

    def create_report(self, onweb=True):
        self._create_report()
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
