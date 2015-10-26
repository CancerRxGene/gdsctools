# -*- python -*-
#
#  This file is part of GDSCtools software
#
#  Copyright (c) 2013-2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt/cellnopt
#
##############################################################################
import os
import shutil

import easydev
import pandas as pd


import bs4

class HTMLTable(object):
    def __init__(self, df, name, **kargs):
        self.df = df.copy() # because we will change its contents possibly
        self.name = name
        self.pd_options = {
                'max_colwidth': -1,
                'precision': 2}

    def to_html(self, index=False, escape=False, header=True, **kargs):
        _buffer = {}
        for k, v in self.pd_options.items():
            # save the current option
            _buffer[k] = pd.get_option(k)
            # set with user value
            pd.set_option(k, v)
    
        table = self.df.to_html(escape=escape, header=header, index=index, 
                **kargs)

        # get back to default options
        for k, v in _buffer.items():
            pd.set_option(k, v)
        return table

    def add_bgcolor(self, colname, cmap, mode='absmax', threshold=None):
        """

        add a background color style by adding <p> tags and bgcolor
        This is apply on one column. 
        The color are set according to the colormap provided (cmap)
        and a normalisation of the data defined by the mode. cmap values
        are between 0 and 1 so, let us normalise the data in that range as
        well. Then, we can easily get the hex values. If the mode is absmax, 
        the max is the abs max and data is scaled between 0 and 1.
        If you have only positive values then, data is between 0.5 and 1.
        The clip mode means that data are positives and 
        """
        from colormap import rgb2hex 

        data = self.df[colname].values
        if len(data) == 0:
            return
        if mode == 'clip':
            data = [min(x, 2)/2. for x in data]
        elif mode == 'absmax':
            m = abs(data.min())
            M = abs(data.max())
            M = max([m,M])
            data = (data / M + 1)/2.
        elif mode == 'max':
            data = data/float(data.max())

        # the expected RGB values for a given data point 
        rgbcolors = [cmap(x)[0:3] for x in data]
        hexcolors = [rgb2hex(*x, normalised=True) for x in rgbcolors]

        # need to read original data again
        data = self.df[colname].values
        # need to set precision since this is going to be a text not a number
        # so pandas will not use the precision for those cases:
        data = [easydev.precision(x, self.pd_options['precision']) for x in data]
        html_formatter = '<p style="background-color:{0}">{1}</p>'
        self.df[colname] = [html_formatter.format(x,y) 
                for x,y in zip(hexcolors, data)]

    def add_href(self, colname, urls=None):
        if urls is not None:
            raise NotImplementedError
        else:
            self.df[colname] = self.df[colname].apply(lambda x:
              '<a href="{0}.html">{1}</a>'.format(x,x))


class Report(object):

    def __init__(self, filename='index.html', directory='report',
                 overwrite=True, verbose=True, dependencies=True):
        self.analysis = 'anova'
        from gdsctools import version
        self.version = version
        self.report_directory = directory
        self.filename = filename
        self.sections = []
        self.section_names = []
        self.add_dependencies = False
        self.pretoc = None

    def show(self):
        from browse import browse as bs
        bs(self.report_directory + os.sep + self.filename)

    def close_body(self):
        return "</body>"

    def close_html(self):
        return "</html>"

    def get_footer(self):
        return self.close_body() + "\n" + self.close_html()

    def _init_report(self, directory=None):
        """create the report directroy and return the directory name"""
        self.sections = []
        self.section_names = []

        if directory is None:
            directory = self.report_directory

        # if the directory already exists, print a warning
        try:
            os.mkdir(directory)
            print("Created directory {}".format(directory))
        except Exception:
            pass
            # already exists
            #txt = "Existing directory {}. Files may be overwritten".format(self.report_directory)
            #if self._overwrite_report is True:
            #    self.warning("Directory %s exists already. Files may be overwritten" % directory)
            #else:
            #    raise IOError('Directory %s exists already. Set
            #    _overwrite_report to True or delete the directory' %
            #    directory)=

        for filename in ["dana.css"]:
            filename = easydev.get_share_file("gdsctools", "data", 
                    filename)
            shutil.copy(filename, directory)
        return directory

    def get_header(self):
        """a possible common header ? """
        params = {"analysis": self.analysis,
                   "version": self.version}
        str_ =  """
 <!DOCTYPE html>
 <html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">
     <head>
     <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
     <title>GDSCtools report</title>
     <link rel="stylesheet" href="dana.css" type="text/css" />
     <script type='text/javascript' src='tools.js'></script>
 </head>

 <body>
  <div class="document" id="unset">

     <h1 class="title">ANOVA analysis summary</h1>
     <h2 class="subtitle", id="unset2">Report created with <br>gdsctools</br> (version %(version)s)</h2>
     <p>See <a href="https://www.github.com/CancerRxGene/gdsctools">GDSCtools github page</a> for downloads and documentation.</p>
     <br>
     Go back to <a href="index.html">main page</a>.<br>
     
     """ % params
        return str_


    def get_time_now(self):
        import datetime
        import getpass
        username = getpass.getuser()
        # this is not working on some systems: os.environ["USERNAME"]
        msg = '<div class="date">Created on ' + str(datetime.datetime.now())
        msg +=  " by " + username +'</div>'
        return msg

    def get_table_dependencies(self):
        """Returns dependencies of the pipeline into a HTML/XML table

        dependencies are the python dependencies as returned by pkg_resource.
        additionally, r dependencies added in :attr:`dependencies` are also added.

        """

        dependencies = easydev.get_dependencies('gdsctools')
        names = [x.project_name for x in dependencies]
        versions = [x.version for x in dependencies]
        links = ["""https://pypi.python.org/pypi/%s"""%p for p in names]

        df = pd.DataFrame({
            'package': ["""<a href="%s">%s</a>"""%(links[i],p)
                for i,p in enumerate(names)],
            'version':versions})

        table = HTMLTable(df, name="dependencies", escape=False)

        return table

    def add_pretoc(self,content):
        self.pretoc = content

    def add_rawhtml(self, content):
        self.sections.append(content)
        #self.section_names.append(None)

    def add_section(self, content, title, references=[], position=None):

        reftxt = self._create_references(references)
        section = """<div class="section" id="%(id)s">
        <h2> <a class="toc-backref" href="#id%(index)s">%(title)s</a></h2>

        %(references)s\n
        %(content)s
    </div>
        """ % {'title':title, 'references':reftxt,'content':content,
               'id': title.replace(" ", "_"), 'index':len(self.sections)+1}
        # check that it is correct
        if position is not None:
            self.sections.insert(position, section)
            self.section_names.insert(position, title)
        else:
            self.sections.append(section)
            self.section_names.append(title)

    def get_toc(self):
        toc = """<div class="contents local topic" id="contents">
        <ul class="simple">"""
        for i, name in enumerate(self.section_names):
            if name is None:
                continue
            toc += """<li>
%(i)s - <a class="reference internal" href="%(href)s" id="%(id)s">  %(name)s</a>
</li>""" % {'i':i+1, 'name':name, 'href':"#"+name.replace(" ", "_"), 'id':'id%s' % str(i+1)}
        toc += """</ul>\n</div>"""
        return toc

    def _create_parameters(self):
        raise NotImplementedError

    def _create_references(self, references):
        if len(references) == 0:
            return ""

        txt = """
        <div class="admonition-documentation admonition">
            <p class="first admonition-title">Documentation</p>
            <ul class="last simple">
            <li>"""

        for ref in references:
            txt += """      <a class="reference external" href=%(url)s>%(title)s</a>""" %  {'url':ref[0],'title':ref[1]}
        txt += """
            </li>
            </ul>
        </div>"""
        return txt

    def write(self ):
        fh =  open(self.report_directory + os.sep + self.filename, "w")

        contents = self.get_header()


        # Get toc should be done here and no more sections should be added
        self.add_section(self.get_toc(), "Contents", position=0)
        # pretoc section must be created after the TOC, so that the position
        # is 0 and it is not taken into account in the TOC.
        if self.pretoc is not None:
            self.add_section(self.pretoc, "", position=0)

        if self.add_dependencies:
            self.add_section(self.get_table_dependencies().to_html(),
                'Dependencies')

        for i, section in enumerate(self.sections):
            if i == 0 or self.section_names[i] == 'Contents':
                contents += section
            else:
                contents += section.replace("<h2>", "<h2> %s - " %i, 1)

        contents += "<hr>" + self.get_time_now()
        contents += self.get_footer()

        contents = bs4.BeautifulSoup(contents).prettify()
        fh.write(contents)
        fh.close()


    def report(self, browse=True):
        self._create_report()
        self.write()
        if browse:
            self.show()

    def _create_report(self):
        raise NotImplementedError






