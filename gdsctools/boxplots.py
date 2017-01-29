# -*- python -*-
# -*- coding utf-8 -*-

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
import pandas as pd
import scipy
import pylab

from gdsctools import boxswarm

import jinja2
from jinja2 import Environment, PackageLoader
from easydev import get_package_location

import colorlog as logger

__all__ = ['BoxPlots', 'BoxPlotsJS']


class BoxPlots(object):
    """Box plot for a given association of drug versus genomic feature


    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import ANOVA, ic50_test
        from gdsctools.boxplots import BoxPlots

        gdsc = ANOVA(ic50_test)

        # Perform the entire analysis
        odof = gdsc._get_one_drug_one_feature_data(1047, 'TP53_mut')

        # Plot volcano plot of pvalues versus signed effect size
        bx = BoxPlots(odof)
        bx.boxplot_association()


    If the **gdsc** analysis has the MSI factor and tissue factor on, then
    additional plots can  be created using :meth:`boxplot_pancan`.


    Note that :attr:`odof` in the example above is a dictionary
    with the following keys:

    - drug_name
    - feature_name
    - masked_tissue: a dataframe with cosmic ids as index and 1 column of
      tissues names.
    - Y: list with the IC50s
    - masked_features: a dataframe with cosmic ids as index and 1 column of
      masked feature  (1/0)
    - masked_msi: same as masked_features
    - negatives: subset of the IC50s corresponding to positive feature
    - positives: subst of the IC50s corresponding to negative feature

    .. seealso:: :class:`gdsctools.boxswarm.BoxSwarm`

    """
    def __init__(self, odof, fontsize=20, savefig=False, directory='.'):
        """.. rubric:: Constructor

        """
        #: dictionary as returned by ANOVA._get_one_drug_one_feature_data
        self.odof = odof

        #: fontsize for the plots
        self.fontsize = fontsize

        #: boolean to save figure
        self.savefig = savefig

        #: directory where to save the figure.
        self.directory = directory

        #: linewidth used in the plots
        self.lw = 3

        self.drug = self.odof.drug_id
        self.feature = self.odof.feature_name.replace("_", " ")


    def boxplot_pancan(self, mode, fignum=1, title_prefix=''):
        """Create boxplot related to the MSI factor or Tissue factor

        :param mode: either set to **msi** or **tissue**

        """
        assert mode in ['tissue', 'msi', "media"]

        results = self._get_boxplot_data(mode)
        if results is None:
            logger.warning("No tissue with at least 2 pos and 2 neg found (no image created).")
            return

        fig = pylab.figure(fignum)
        oldsize = fig.get_size_inches()

        pylab.clf()  # or close ?
        data, names, significance = results
        N = len(names)
        if N <= 2: # msi or 2 tissues
            fontsize = self.fontsize
        elif N<=14:
            fontsize = max(4, int(self.fontsize - (N-2.)/(self.fontsize-4.)))
        else:
            fontsize = max(4, int(self.fontsize/1.4))

        bb = boxswarm.BoxSwarm(data, names, fontsize=fontsize)

        bb.xlabel = r'%s log(IC50)' % self.drug
        if mode == 'tissue':
            bb.title = 'FEATURE/Cancer-type interactions'
        elif mode == 'msi':
            bb.title = 'FEATURE/MS-instability interactions'
        elif mode == "media":
            bb.title = 'FEATURE/Media interactions'
        ax = bb.plot(vert=False)
        # get info from left axis
        common_ylim = ax.get_ylim()
        common_ticks = ax.get_yticks()

        self.ax = ax.twinx()
        self.ax.set_ylim(common_ylim)
        self.ax.set_yticks(common_ticks)
        self.ax.set_yticklabels([str(len(this))+" " for this in data],
                fontsize=fontsize/1.4)
        try:
            pylab.tight_layout()
        except:
            pass

        if self.savefig is True:
            filename = self.directory + os.sep
            filename += 'ODOF_{}_DRUG_{}____{}'.format(mode,
                    self.drug, self.feature)
            fig.set_size_inches(14, 16)
            pylab.savefig(filename + '.png', bbox_inches='tight')
            fig.set_size_inches(oldsize)
            fig.canvas.draw()
            #pylab.savefig(filename + '.svg', bbox_inches='tight')

    def boxplot_association(self, fignum=1):
        """Boxplot of the association (negative versus positive)

        :param fignum: number of the figure
        """
        pylab.figure(fignum)
        pylab.clf()

        # the plot itself
        boxswarm.boxswarm(
                {'pos': self.odof.positives, 'neg': self.odof.negatives},
                lw=self.lw, fontsize=self.fontsize)


        if self.odof.drug_name:
            drug_name = "%s (%s)" % (self.drug, self.odof.drug_name)
        else:
            drug_name = self.drug

        pylab.title('Individual association\n {0} versus {1}'.format(drug_name,
            self.feature), fontsize=self.fontsize)
        pylab.ylabel("{0} logIC50".format(self.drug),
                fontsize=self.fontsize)

        try:pylab.tight_layout()
        except: pass

        if self.savefig is True:
            filename = self.directory + os.sep
            filename += 'ODOF_all_DRUG_{}____{}'.format(self.drug,
                    self.feature)
            pylab.savefig(filename + '.png', bbox_inches='tight')

    def _get_boxplot_data(self, mode='tissue'):
        # should be called by anova_one_drug_one_feature
        # since masked_tissue, masked_ic50 attributes must
        # be populated.
        assert mode in ['tissue', 'msi', "media"]

        # Let us use Pandas, this will be easier
        try:
            df = pd.DataFrame({
                 'ic50': self.odof.Y,
                 'media': self.odof.masked_media,
                 'feature': self.odof.masked_features,
                 'msi': self.odof.masked_msi.values,
                'tissue': self.odof.masked_tissue.values})
        except:
            df = pd.DataFrame({
                 'ic50': self.odof.Y,
                 'tissue': self.odof.masked_tissue,
                 'feature': self.odof.masked_features,
                 'msi': self.odof.masked_msi.values})

        if mode == 'tissue':
            if 'msi' in df.columns:
                df.drop('msi', inplace=True, axis=1)
            if "media" in df.columns:
                df.drop('media', inplace=True, axis=1)
        elif mode == 'msi':
            if 'tissue' in df.columns:
                df.drop('tissue', inplace=True, axis=1)
            if "media" in df.columns:
                df.drop('media', inplace=True, axis=1)
        elif mode == 'media':
            if 'tissue' in df.columns:
                df.drop('tissue', inplace=True, axis=1)
            if 'msi' in df.columns:
                df.drop('msi', inplace=True, axis=1)

        groups = df.groupby(['feature', mode])

        # counts items in each category and fill with NA
        counts = groups.count().unstack().fillna(0)


        # if positive or negative for a combo, is not>=2, drop it
        # pandas 0.16.2
        cc = (counts >= 2).all()
        # create a groups structure
        categories = list(cc.unstack().columns[cc])

        """
        # Seems to be fixed (May 2016)
        try:
            # pandas 0.16.2
            cc = (counts >= 2).all()
            # create a groups structure
            categories = list(cc.unstack().columns[cc])
        except:
            # pandas 0.13 for the doc only
            categories = list(df.tissue.unique())[0:10]
        """

        groups = df.query(mode + ' in @categories',
                engine='python').groupby([mode, 'feature'])

        # TODO; move all this if block into a method
        # figure out the delta between pos and neg
        means = groups.mean().unstack(mode)

        if len(means):  # need 2 values
            delta = means.ix[0] - means.ix[1]
            try:
                # new pandas v0.17
                delta.sort_values(inplace=True)
            except Exception as err:
                # conda for py3.3 is 0.16.2 where sort_values sdoes not exists
                # Useful also for readthedocs
                delta.sort(inplace=True)

            significance = {}
            data = []
            names = []
            for category in delta.ix['ic50'].index:
                prefix_query = mode+"==@category"
                neg = df.query(prefix_query+' and feature==0',
                        engine='python')['ic50']
                pos = df.query(prefix_query+' and feature==1',
                        engine='python')['ic50']
                # HERE in the original code, equal_var is False. why ?
                res = scipy.stats.ttest_ind(neg, pos, equal_var=False)

                #print res[1], self.odof.ttest
                # this should be computed outside ??
                significance[category] = res[1] # p-values
                data.append(neg.values)
                data.append(pos.values)
                if mode == 'tissue':
                    name = category
                elif mode == 'msi':
                    if category == 0:
                        name = 'MSI-stable'
                    elif category == 1:
                        name = 'MSI-unstable'
                elif mode == "media":
                    name = category

                for this in [0.05, 0.01, 0.001]:
                    if significance[category] < this:
                        name = '*' + name
                names.append(name + ' neg')
                names.append(name + ' pos')
            return (data, names, significance)
        else:
            return None


class BoxPlotsJS(BoxPlots):
    def __init__(self, odof, fontsize=20, savefig=False, directory='.'):
        super(BoxPlotsJS, self).__init__(odof, fontsize=fontsize,
                savefig=savefig, directory=directory)

    def get_html_association(self):
        env = Environment()
        env.loader = jinja2.FileSystemLoader(
                       get_package_location("gdsctools")
                       + "/gdsctools/data/templates/")
        template = env.get_template("boxplot_association.html")

        jinja = {}
        N = len(self.odof.negatives) + len(self.odof.positives)
        jinja["title"] = "Individual Association"

        if self.odof.drug_name:
            drug_name = "%s (%s)" % (self.drug, self.odof.drug_name)
            jinja["subtitle"] = "%s versus %s" % (drug_name, self.feature)
        else:
            jinja["subtitle"] = "%s versus %s" % (self.drug, self.feature)
        jinja['factor'] = ["neg"] * len(self.odof.negatives) + ["pos"] * len(self.odof.positives)
        jinja["data"] = self.odof.negatives.tolist() + self.odof.positives.tolist()
        jinja["smps"] = [str(this) for this in self.odof.indices]
        jinja["subject"] = [str(this) for this in self.odof.indices]

        jinja['ylabel'] = '"logIC50"'

        html = template.render(jinja)
        return html

    def get_html_media(self):
        env = Environment()
        env.loader = jinja2.FileSystemLoader(
                       get_package_location("gdsctools")
                       + "/gdsctools/data/templates/")
        template = env.get_template("boxplot_media.html")

        data = self._get_boxplot_data("media")
        if data is None:
            return ""
        # Show from bottom to top
        labels = data[1][::-1]
        data = data[0][::-1]

        jinja = {}
        N = len(self.odof.negatives) + len(self.odof.positives)
        jinja["title"] = "FEATURE/Media interactions"
        jinja["subtitle"] = "%s versus %s" % (self.drug, self.feature)
        factor = []
        for i, thisdata in enumerate(data):
            factor.extend( [labels[i]] * len(thisdata))
        jinja['sign'] = [x.split()[1] for  x in factor]
        jinja['status'] = [x.split()[0] for  x in factor]
        jinja["data"] = list(pylab.flatten([list(this) for this in data]))
        jinja['xlabel'] = '"logIC50"'

        html = template.render(jinja)
        return html


    def get_html_msi(self):
        env = Environment()
        env.loader = jinja2.FileSystemLoader(
                       get_package_location("gdsctools")
                       + "/gdsctools/data/templates/")
        template = env.get_template("boxplot_msi.html")

        data = self._get_boxplot_data("msi")
        if data is None:
            return ""
        # Show from bottom to top
        labels = data[1][::-1]
        data = data[0][::-1]

        jinja = {}
        N = len(self.odof.negatives) + len(self.odof.positives)
        jinja["title"] = "FEATURE/MS-instability interactions"
        jinja["subtitle"] = "%s versus %s" % (self.drug, self.feature)
        factor = []
        for i, thisdata in enumerate(data):
            factor.extend( [labels[i]] * len(thisdata))
        jinja['sign'] = [x.split()[1] for  x in factor]
        jinja['status'] = [x.split()[0] for  x in factor]
        jinja["data"] = list(pylab.flatten([list(this) for this in data]))
        jinja['xlabel'] = '"logIC50"'

        html = template.render(jinja)
        return html

    def get_html_tissue(self):
        env = Environment()
        env.loader = jinja2.FileSystemLoader(
                       get_package_location("gdsctools")
                       + "/gdsctools/data/templates/")
        template = env.get_template("boxplot_tissue.html")

        data = self._get_boxplot_data("tissue")
        if data is None:
            return ""
        # Show from bottom to top
        labels = data[1][::-1]
        data = data[0][::-1]

        jinja = {}
        N = len(self.odof.negatives) + len(self.odof.positives)
        jinja["title"] = "FEATURE/MS-instability interactions"
        jinja["subtitle"] = "%s versus %s" % (self.drug, self.feature)
        factor = []
        for i, thisdata in enumerate(data):
            factor.extend( [labels[i]] * len(thisdata))
        jinja['sign'] = [x.split()[1] for  x in factor]
        jinja['status'] = [x.split()[0] for  x in factor]
        jinja["data"] = list(pylab.flatten([list(this) for this in data]))
        jinja['xlabel'] = '"logIC50"'

        if len(labels)/2 >= 10:
            jinja["minTextSize"] = 10

        html = template.render(jinja)
        return html

    


