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
    """ANOVA analysis of the IC50 vs Feature matrices


    """
    def __init__(self, odof, fontsize=20, savefig=False, directory='.',
            verbose='INFO'):
        """.. rubric:: Constructor

        """
        super(BoxPlots, self).__init__(level=verbose)
        self.odof = odof
        self.fontsize = fontsize
        self.savefig = savefig
        self.directory = directory

    def boxplot_pancan(self, mode, fignum=1, title_prefix=''):
        assert mode in ['tissue', 'msi']
        drug_name = self.odof.drug_name.replace("_", "\_")
        #feature_name = odof.feature_name.replace("_", "\_")

        results = self._get_boxplot_data(mode)
        if results is None:
            self.logging.info("INFO: no tissue with at least 2 pos and 2 neg found. " + "No image created.")
            return

        pylab.figure(fignum)
        pylab.clf()  # or close ?
        data, names, significance = results

        bb = boxswarm.BoxSwarm(data, names)
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
        self.ax.set_yticklabels([len(this) for this in data])

        #pylab.tight_layout()
        if self.savefig is True:
            filename = self.directory + os.sep
            filename += 'ODOF_{}_{}____{}'.format(mode,
                    self.odof.drug_name, self.odof.feature_name)
            pylab.savefig(filename + '.png', bbox_inches='tight')
            #pylab.savefig(filename + '.svg', bbox_inches='tight')

    def boxplot_association(self, fignum=1):
        data = self.odof
        pylab.figure(fignum)
        pylab.clf()
        # aliases
        drug_name = data.drug_name.replace("_", "\_")
        feature_name = data.feature_name.replace("_", "\_")
        fontsize = self.fontsize

        # the plot itself
        boxswarm.boxswarm({'pos': data.positives, 'neg': data.negatives},
                lw=3, fontsize=self.fontsize)

        pylab.title('Individual association\n {0} versus {1}'.format(drug_name,
            feature_name), fontsize=fontsize)
        pylab.ylabel("{0} logIC50".format(drug_name),
                fontsize=fontsize)

        pylab.tight_layout()
        if self.savefig is True:
            filename = self.directory + os.sep
            filename += 'ODOF_all_{}____{}'.format(data.drug_name,
                    data.feature_name)
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
        cc = (counts >= 2).all()

        # create a groups structure 
        categories = list(cc.unstack().columns[cc])
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

