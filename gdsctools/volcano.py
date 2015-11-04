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
"""Volcano plot utilities"""
import pandas as pd
import pylab
import numpy as np
import easydev

from easydev import Progress, AttrDict
from gdsctools.tools import Savefig


__all__ = ['VolcanoANOVA']


class VolcanoANOVA(Savefig):
    """Utilities related to volcano plots


    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import ANOVA, ic50_test, VolcanoANOVA
        an = ANOVA(ic50_test)

        an = anova.ANOVA(ic50_test)
        an.set_cancer_type('breast')
        df = an.anova_all()
        v = VolcanoANOVA(df)


    .. note:: **for developers**: there is experimental code to
        connect the volcano plot using mpld3 package and javascript.
    """
    def __init__(self, data, sep="\t", settings=None):
        """.. rubric:: Constructor

        :param data: a dataframe are returned by e.g., :meth:`ANOVA.anova_all`
        :param settings: an instance of :class:`ANOVASettings`

        """
        super(VolcanoANOVA, self).__init__()

        if isinstance(data, str):
            self.df = pd.read_csv(data, sep=sep)
        else:
            self.df = data.copy()

        # this is redundant coul reuse the input ??
        if settings is None:
            from gdsctools.settings import ANOVASettings
            self.settings = ANOVASettings()
        else:
            self.settings = AttrDict(**settings)
        self.colname_drugid = 'DRUG_ID'

        self.drugs = set(self.df[self.colname_drugid]) # set on values
        self.varname_pvalue = 'FEATURE_ANOVA_pval'
        self.varname_qvalue = 'ANOVA_FEATURE_FDR_%'

        # intensive calls made once for all
        self.groups_by_drugs = self.df.groupby(self.colname_drugid).groups
        self.groups_by_features = self.df.groupby('FEATURE').groups

    def selector(self, df, Nbest=1000, Nrandom=1000, inplace=False):
        """Select first N best rows and N random ones

        Used to select a representative set of rows and the best rows (in
        a sorted dataframe). This is used to create the volcano plots quickly.

        """
        if len(df)<Nbest:
            return df
        Nmax =  Nbest + Nrandom
        N  = len(df)
        if N > Nbest:
            x = range(Nbest, N)
            pylab.shuffle(x)
            n2pick = min(N, Nmax) - Nbest
            indices = range(0, Nbest) + x[0:n2pick]
        else:
            indices = range(0,Nbest)
        df = df.ix[indices]
        if inplace is True:
            self.df = df
        else:
            return df

    def volcano_plot_all_drugs(self):
        """Volcano plot for each drug and savefig

        :param df: output of :meth:`anova_all`
        """
        drugs = list(self.df[self.colname_drugid].unique())
        pb = Progress(len(drugs), 1)
        pylab.ioff()
        for i, drug in enumerate(drugs):
            self.volcano_plot_one_drug(drug)
            if self.settings.savefig is True:
                self.savefig("volcano_%s.png" % drug)
            pb.animate(i+1)
        pylab.ion()

    def volcano_plot_all_features(self):
        """Volcano plot for each feature

        :param df: output of :meth:`anova_all`

        Takes about 10 minutes for 265 drugs and 677 features
        """
        features = list(self.df['FEATURE'].unique())
        print('Creating image for each feature (using all drugs)')
        pb = Progress(len(features), 1)
        for i, feature in enumerate(features):
            self.volcano_plot_one_feature(feature)
            if self.settings.savefig is True:
                self.savefig("volcano_%s.png" % feature)
            pb.animate(i+1)

    def volcano_plot_all(self):
        # no annotations for all features.
        # this is slow, we can drop non relevant data
        data = self._get_volcano_sub_data('ALL')
        data['annotation'] = ['' for x in range(len(data))]

        self.volcano_plot(data, title='all drugs')
        if self.settings.savefig is True:
            self.savefig("volcano_all.png")

    def _get_fdr_from_pvalue_interp(self, pvalue):
        pvalue += 1e-15
        qvals = self.df[self.varname_qvalue]
        pvals = self.df[self.varname_pvalue]
        ya = qvals[pvals < pvalue].max()
        yb = qvals[pvals > pvalue].min()
        xa = pvals[pvals < pvalue].max()
        xb = pvals[pvals > pvalue].min()
        dx = xb - xa
        dy = yb - ya
        yc = ya + dy * (pvalue - xa) / dx
        return yc

    def _get_pvalue_from_fdr(self, fdr):
        qvals = self.df[self.varname_qvalue]
        pvals = self.df[self.varname_pvalue]
        if isinstance(fdr, list):
            pvalues = [pvals[qvals < this].max() for this in fdr]
            return pvalues
        else:
            return pvals[qvals < fdr].max()

    def _get_pvalue_from_fdr_interp(self, fdr):
        # same as get_pvalue_from_fdr but with a linear inerpolation
        fdr += 1e-15
        qvals = self.df[self.varname_qvalue]
        pvals = self.df[self.varname_pvalue]
        ya = pvals[qvals < fdr].max()
        yb = pvals[qvals > fdr].min()
        xa = qvals[qvals < fdr].max()
        xb = qvals[qvals > fdr].min()
        dx = xb - xa
        dy = yb - ya
        xc = fdr
        yc = ya + dy * (xc - xa) / dx
        return yc

    def _get_volcano_global_data(self):
        # using all data
        minN = self.df['N_FEATURE_pos'].min()
        maxN = self.df['N_FEATURE_pos'].max()
        pvalues = self._get_pvalue_from_fdr(self.settings.FDR_threshold)
        return {'minN': minN, 'maxN': maxN,
                'pvalues': (self.settings.FDR_threshold, pvalues)}

    def volcano_plot_one_feature(self, feature):
        """Volcano plot for one feature (all drugs)"""
        data = self._get_volcano_sub_data('FEATURE', feature)
        self.volcano_plot(data, title=feature)

    def _get_volcano_sub_data(self, mode, target=None):
        # using data related to the given drug

        # groups created in the constructor once for all
        if mode == 'DRUG_ID':
            subdf = self.df.ix[self.groups_by_drugs[target]]
            texts = subdf['FEATURE']
        elif mode == 'FEATURE':
            subdf = self.df.ix[self.groups_by_features[target]]
            texts = subdf['DRUG_ID']
        elif mode == 'ALL':
            # nothing to do, get all data
            subdf = self.df
            texts = subdf['FEATURE'] # TODO + drug
        else:
            raise ValueError("mode parameter must be in [FEATURE, DRUG_ID, ALL]")

        # replaced by groups created in the constructor
        #subdf = self.df[self.df[mode] == target]

        deltas = subdf['FEATURE_deltaMEAN_IC50']
        effects = subdf['FEATURE_IC50_effect_size']
        signed_effects = list(np.sign(deltas) * effects)
        qvals = list(subdf[self.varname_qvalue])
        pvals = list(subdf[self.varname_pvalue])

        colors = []
        annotations = []

        data = pd.DataFrame(index=range(len(qvals)))
        data['pvalue'] = pvals
        data['signed_effect'] = signed_effects
        data['Feature'] = list(subdf['FEATURE'])
        data['Drug'] = list(subdf['DRUG_ID'])
        data['text'] = texts.values
        data['FDR'] = subdf['ANOVA_FEATURE_FDR_%']
        annotations = []

        # just an alias
        FDR_threshold = self.settings.FDR_threshold
        if self.settings.analysis_type == 'PANCAN':
            for sign, qval, pval in zip(signed_effects, qvals, pvals):
                if sign <= -self.settings.effect_threshold and \
                        qval <= FDR_threshold:
                    colors.append('green')
                    annotations.append(True)
                elif sign >= self.settings.effect_threshold and \
                        qval <= FDR_threshold:
                    colors.append('red')
                    annotations.append(True)
                else:
                    colors.append('black')
                    annotations.append(False)
        else:
            for delta, qval, pval in zip(deltas, qvals, pvals):
                if pval <= self.settings.pvalue_threshold and \
                        qval <= FDR_threshold and delta < 0:
                    colors.append('green')
                    annotations.append(True)
                elif pval <= self.settings.pvalue_threshold and \
                        qval <= FDR_threshold and delta > 0:
                    colors.append('red')
                    annotations.append(True)
                else:
                    colors.append('black')
                    annotations.append(False)

        # here we normalise wrt the drug. In R code, normalised
        # my max across all data (minN, maxN)
        markersize = subdf['N_FEATURE_pos'] / subdf['N_FEATURE_pos'].max()
        markersize = list(markersize*800)
        markersize = [x if x > 50 else 50 for x in markersize]

        data['color'] = colors
        data['annotation'] = annotations
        data['markersize'] = markersize
        return data

    def volcano_plot_one_drug(self, drug_id):
        """Volcano plot for one drug and all genomic features

        :param df: output of :meth:`anova_all`

        """
        assert drug_id in self.drugs, 'unknown drug name'
        # needs to run :meth:`anova_all` first
        # using all data, get the FDR limits
        data = self._get_volcano_sub_data(self.colname_drugid, drug_id)
        self.volcano_plot(data, title=drug_id)

    def volcano_plot(self, data, title=''):
        """Plots signed effects versus pvalues

        THis is for internal usage but maybe useful elsewhere.

        """
        colors = list(data['color'].values)
        pvalues = data['pvalue'].values
        signed_effects = data['signed_effect'].values
        markersize = data['markersize'].values

        Y = -np.log10(list(pvalues)) # should be cast to list ?
        # This is horrendously slow as compared to plot()
        # However, plot cannot take different sizes/colors
        # Takes about 0.36 s per call and half the time
        # is spent in the scatter() call

        pylab.close(1)
        fig = pylab.figure(num=1)
        ax = pylab.axes(axisbg='#EEEEEE')
        ax = fig.gca()
        X = [easydev.precision(x, digit=2) for x in signed_effects]
        Y = [easydev.precision(y, digit=2) for y in Y]

        # black markers will have the same size
        # and will not be labelled
        scatter = ax.scatter(X, Y, s=markersize,
                alpha=0.3, c=colors,
                linewidth=0, picker=True)
        #pylab.plot(list(signed_effects), Y, markersize=5,
        #            alpha=0.4, c='grey',
        #            linewidth=0)

        m = abs(signed_effects.min())
        M = abs(signed_effects.max())
        pylab.xlabel("Signed effect size", fontsize=self.settings.fontsize)
        pylab.ylabel('-log10(pvalues)', fontsize=self.settings.fontsize)
        l = max([m, M]) * 1.1
        pylab.xlim([-l, l])
        ax.grid(color='white', linestyle='solid')

        #self.stats = self._get_volcano_global_data()

        fdr = self.settings.FDR_threshold
        pvalue = self._get_pvalue_from_fdr(fdr)
        ax.axhline(-np.log10(pvalue), linestyle='--',
            color='red', alpha=1, label="FDR %s " %  fdr + " \%")
        #pvalue = self._get_pvalue_from_fdr_interp(fdr)
        #ax.axhline(-np.log10(pvalue), linestyle='--',
        #    color='red', alpha=1, label="FDR %s " %  fdr + " \%")

        pvalue = self._get_pvalue_from_fdr(10)
        #pvalue = self._get_pvalue_from_fdr_interp(10)
        ax.axhline(-np.log10(pvalue), linestyle='-.',
            color='red', alpha=1, label="FDR 10 \%")

        #pvalue = self._get_pvalue_from_fdr_interp(1)
        pvalue = self._get_pvalue_from_fdr(1)
        ax.axhline(-np.log10(pvalue), linestyle=':',
            color='red', alpha=1, label="FDR 1 \%")

        pvalue = self._get_pvalue_from_fdr_interp(0.01)
        pvalue = self._get_pvalue_from_fdr(0.01)
        ax.axhline(-np.log10(pvalue), linestyle='--',
            color='black', alpha=1, label="FDR 0.01 \%")

        pylab.ylim([0, pylab.ylim()[1]*1.2]) # times 1.2 to put the legend

        ax.axvline(0, color='gray', alpha=0.5)
        axl = pylab.legend(loc='upper left')
        axl.set_zorder(-1) # in case there is a circle behind the legend.

        #self.ax = ax
        #self.axx = ax.twinx()
        #self.common_ticks = ax.get_yticks()
        #self.common_ylim = ax.get_ylim()
        #pvals = self.df[self.varname_pvalue]
        #y1 = pvals.min()
        #y2 = pvals.max()
        #fdr1 = self._get_fdr_from_pvalue_interp(y1)
        #fdr2 = self._get_fdr_from_pvalue_interp(y2-2e-15) # make sure it exists
        #self.axx.set_ylim([fdr2, fdr1])
        #self.axx.set_ylabel('FDR \%', fontsize=self.settings.fontsize)

        # For the static version
        pylab.title("%s" % title.replace("_","\_"))
        labels = []
        """for index, row in data.iterrows():
            text = data[['feature', 'drug',
                    'signed_effect','pvalue']].ix[0].to_frame()
            if row.annotation == True:
                x, y= row.signed_effect, row.pvalue
                ax.text(x, -pylab.log10(y), row.text.replace("_","\_"))
            labels.append(text)
        def onpick(event):
            ind = event.ind
            print('on pick scatter:', ind, np.take(X, ind)[0],
                    np.take(Y, ind)[0])
        #fig.canvas.mpl_connect('pick_event', onpick)
        """

        # for the JS version
        # TODO: for the first 1 to 2000 entries ?
        import mpld3
        labels = []
        for i, row in data[['Drug','Feature', 'FDR']].iterrows():
            label = row.to_frame()
            label.columns = ['Row {0}'.format(i)]
            # .to_html() is unicode; so make leading 'u' go away with str()
            labels.append(str(label.to_html(header=False)))

        css = """
        table{  font-size:0.8em;  }
        th {  color: #ffffff;  background-color: #aaaaaa;  }
        td { color: blue; background-color: #cccccc; }
        """

        tooltip = mpld3.plugins.PointHTMLTooltip(scatter, labels=labels,
                css=css)
        mpld3.plugins.connect(fig, tooltip)

        #mpld3.display()

        self.scatter = scatter
        self.current_fig = fig
        #fh = open('test.html', 'w'); mpld3.save_html(v.current_fig, fh);
        #fh.close()
