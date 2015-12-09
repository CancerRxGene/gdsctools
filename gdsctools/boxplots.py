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

from easydev import Logging

from gdsctools import boxswarm


__all__ = ['BoxPlots']


class BoxPlots(Logging):
    """Box plot for a given association of drug versus genomic feature


    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import ANOVA, ic50_test
        from gdsctools.boxplots import BoxPlots

        gdsc = ANOVA(ic50_test)

        # Perform the entire analysis
        odof = gdsc._get_one_drug_one_feature_data('Drug_1047_IC50', 
            'TP53_mut')

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
    def __init__(self, odof, fontsize=20, savefig=False, directory='.',
            verbose='INFO'):
        """.. rubric:: Constructor

        """
        super(BoxPlots, self).__init__(level=verbose)
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

    def boxplot_pancan(self, mode, fignum=1, title_prefix=''):
        """Create boxplot related to the MSI factor or Tissue factor

        :param mode: either set to **msi** or **tissue**

        """
        assert mode in ['tissue', 'msi']
        drug_name = self.odof.drug_name.replace("_", " ")

        results = self._get_boxplot_data(mode)
        if results is None:
            self.logging.info("INFO: no tissue with at least 2 pos and 2 neg found. " + "No image created.")
            return

        fig = pylab.figure(fignum)
        oldsize = fig.get_size_inches()

        pylab.clf()  # or close ?
        data, names, significance = results
        N = len(names)
        if N<=2: # msi or 2 tissues
            fontsize = self.fontsize 
        else:
            fontsize = max(4, int(self.fontsize - (N-2.)/(self.fontsize-4.)))

        bb = boxswarm.BoxSwarm(data, names, fontsize=fontsize)

        bb.xlabel = r'%s log(IC50)' % drug_name
        if mode == 'tissue':
            bb.title = 'FEATURE/Cancer-type interactions'
        else:
            bb.title = 'FEATURE/MS-instability interactions'
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
            filename += 'ODOF_{}_{}____{}'.format(mode,
                    self.odof.drug_name, self.odof.feature_name)
            fig.set_size_inches(12, 14)
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
        # aliases
        drug_name = self.odof.drug_name.replace("_", " ")
        feature_name = self.odof.feature_name.replace("_", " ")

        # the plot itself
        boxswarm.boxswarm(
                {'pos': self.odof.positives, 'neg': self.odof.negatives},
                lw=self.lw, fontsize=self.fontsize)

        pylab.title('Individual association\n {0} versus {1}'.format(drug_name,
            feature_name), fontsize=self.fontsize)
        pylab.ylabel("{0} logIC50".format(drug_name),
                fontsize=self.fontsize)

        try:pylab.tight_layout()
        except: pass

        if self.savefig is True:
            filename = self.directory + os.sep
            filename += 'ODOF_all_{}____{}'.format(self.odof.drug_name,
                    self.odof.feature_name)
            pylab.savefig(filename + '.png', bbox_inches='tight')

    def _get_boxplot_data(self, mode='tissue'):
        # should be called by anova_one_drug_one_feature
        # since masked_tissue, masked_ic50 attributes must
        # be populated.
        assert mode in ['tissue', 'msi']

        # Let us use Pandas, this will be easier
        df = pd.DataFrame(
            {'tissue': self.odof.masked_tissue.values,
             'ic50': self.odof.Y,
             'feature': self.odof.masked_features,
             'msi': self.odof.masked_msi.values})

        if mode == 'tissue':
            df.drop('msi', inplace=True, axis=1)
        elif mode == 'msi':
            df.drop('tissue', inplace=True, axis=1)

        groups = df.groupby(['feature', mode])

        # counts items in each category and fill with NA
        counts = groups.count().unstack().fillna(0)

        # if positive or negative for a combo, is not>=2, drop it
        try:
            # pandas 0.16.2
            cc = (counts >= 2).all()
            # create a groups structure 
            categories = list(cc.unstack().columns[cc])
        except:
            # pandas 0.13 for the doc only
            categories = list(df.tissue.unique())[0:10]

        groups = df.query(mode + ' in @categories',
                engine='python').groupby([mode, 'feature'])

        # TODO; move all this if block into a method
        # figure out the delta between pos and neg
        means = groups.mean().unstack(mode)


        if len(means):
            delta = means.ix[0] - means.ix[1]
            try:
                # new pandas v0.17
                delta.sort_values(inplace=True)
            except:
                # sort_values not in anaconda for py3.3
                # FIXME does not work in readthedocs.
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

                for this in [0.05, 0.01, 0.001]:
                    if significance[category] < this:
                        name = '*' + name
                names.append(name + ' neg')
                names.append(name + ' pos')
            return (data, names, significance)
        else:
            return None
