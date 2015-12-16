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
"""Volcano plot utilities

The :class:`VolcanoANOVA` is used in the creation of the report but
may be used by a usr in a Python shell.

"""
import os

import pandas as pd
import pylab
import numpy as np
import easydev


from easydev import Progress, AttrDict
from gdsctools.tools import Savefig


__all__ = ['VolcanoANOVA']


class VolcanoANOVA(object):
    """Utilities related to volcano plots

    This class is used in :mod:`gdsctools.anova` but can also
    be used independently as in the example below.

    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import ANOVA, ic50_test, VolcanoANOVA
        an = ANOVA(ic50_test)

        # retrisct analysis to a tissue to speed up computation
        an.set_cancer_type('lung_NSCLC')

        # Perform the entire analysis
        results = an.anova_all()

        # Plot volcano plot of pvalues versus signed effect size
        v = VolcanoANOVA(results)
        v.volcano_plot_all()

    .. note:: Within an IPython shell, you should be able to click
        on a circle and the title will be updated
        with the name of the drug/feature and FDR value.

    .. note:: (**for developers**) A javascript version is also
        created on the fly using mpld3 library. It is used in the
        creation of the HTML report but one can use it as well in an
        ipython notebook::

            # The **v** instance is as created in the example above
            # Then, type the following code to create an HTML with the
            # javascript plot embedded.
            import mpld3
            htmljs = mpld3.fig_to_html(v.current_fig)
            fh = open('volcano_doc.html', 'w')
            fh.write(htmljs)
            fh.close()

    There are 5 methods to plot volcano plots depending on what you want to see

    - :meth:`volcano_plot_all` as above plots all associations
    - :meth:`volcano_plot_all_drugs` creates a volcano plot for each drug and
      save it into a PNG file. This method calls :meth:`volcano_plot_one_drug`.
    - :meth:`volcano_plot_all_features` creates a volcano plot for each feature
      and save it into a PNG file. This method calls
      :meth:`volcano_plot_one_feature`.

    """
    _colname_pvalue = 'ANOVA_FEATURE_pval'
    _colname_qvalue = 'ANOVA_FEATURE_FDR'
    _colname_drugid = 'DRUG_ID'
    _colname_feature = 'FEATURE'
    _colname_deltas = 'FEATURE_delta_MEAN_IC50'
    _colname_effect_size = 'FEATURE_IC50_effect_size'
    _colname_N_feature_pos = 'N_FEATURE_pos'

    def __init__(self, data, sep="\t", settings=None):
        """.. rubric:: Constructor

        :param data: an :class:`~gdsctools.anova.ANOVAResults` instance
            or a dataframe with the proper columns names (see below)
        :param settings: an instance of
            :class:`~gdsctools.settings.ANOVASettings`

        Expected column names to be found if a filename is provided::

            ANOVA_FEATURE_pval
            ANOVA_FEATURE_FDR
            FEATURE_delta_MEAN_IC50
            FEATURE_IC50_effect_size
            N_FEATURE_pos
            N_FEATURE_pos
            FEATURE
            DRUG_ID

        If the plotting is too slow, you can use the :meth:`selector` to prune
        the results (most of the data are noise and overlap on the middle
        bottom  area of the plot with little information.

        """
        # a copy since we do may change the data
        try:
            # an ANOVAResults contains a df attribute
            self.df = data.df.copy()
        except:
            # probably a dataframe
            self.df = data.copy()

        # this is redundant could reuse the input ??
        if settings is None:
            from gdsctools.settings import ANOVASettings
            self.settings = ANOVASettings()
        else:
            self.settings = AttrDict(**settings)

        self.figtools = Savefig()
        self.figtools.directory = self.settings.directory

        self.drugs = set(self.df[self._colname_drugid])
        self.features = set(self.df[self._colname_feature])

        # intensive calls made once for all
        self.groups_by_drugs = self.df.groupby(self._colname_drugid).groups
        self.groups_by_features = self.df.groupby(self._colname_feature).groups

    def selector(self, df, Nbest=1000, Nrandom=1000, inplace=False):
        """Select only the first N best rows and N random ones

        Sometimes, there are tens of thousands of associations and future
        analysis will include more features and drugs. Plotting volcano plots
        should therefore be fast and scalable. Here, we provide a naive
        way of speeding up the plotting by selecting only a subset of the data
        made of Nbest+Nrandom associations.

        :param df: the input dataframe with ANOVAResults
        :param int Nbest: how many of the most significant association
            should be kept
        :param int Nrandom: on top of the Nbest significant association,
            set how many other randomly chosen associations are to be kept.
        :return: pruned dataframe

        """
        if len(df) < Nbest:
            return df
        Nmax = Nbest + Nrandom
        N = len(df)
        if N > Nbest:
            x = range(Nbest, N)
            pylab.shuffle(x)
            n2pick = min(N, Nmax) - Nbest
            indices = range(0, Nbest) + x[0:n2pick]
        else:
            indices = range(0, Nbest)
        # indices in the index may not be order
        indices = [df.index[xx] for xx in indices]
        df = df.ix[indices]
        if inplace is True:
            self.df = df
        else:
            return df

    def volcano_plot_all_drugs(self):
        """Create a volcano plot for each drug and save in PNG files

        Each filename is set to **volcano_<drug identifier>.png**
        """
        drugs = list(self.df[self._colname_drugid].unique())
        pb = Progress(len(drugs), 1)
        for i, drug in enumerate(drugs):
            self.volcano_plot_one_drug(drug)
            self.savefig_and_js("volcano_%s.png" % drug, size_inches=(10, 10))
            pb.animate(i+1)

            # This prevent memory leak.
            self.current_fig.canvas.mpl_disconnect(self.cid)
            try:
                import mpld3
                mpld3.plugins.clear(self.current_fig)
            except:
                pass

    def volcano_plot_all_features(self):
        """Create a volcano plot for each feature and save in PNG files

        Each filename is set to **volcano_<feature name>.png**
        """
        features = list(self.df[self._colname_feature].unique())
        print('Creating image for each feature (using all drugs)')
        pb = Progress(len(features), 1)
        for i, feature in enumerate(features):
            self.volcano_plot_one_feature(feature)
            self.savefig_and_js("volcano_%s.png" % feature, 
                    size_inches=(10, 10))
            pb.animate(i+1)

            # This prevent memory leak.
            self.current_fig.canvas.mpl_disconnect(self.cid)
            try:
                import mpld3
                mpld3.plugins.clear(self.current_fig)
            except:
                pass

    def volcano_plot_all(self):
        """Create an overall volcano plot for all associations

        This method saves the picture in a PNG file named **volcano_all.png**.
        """
        # no annotations for all features.
        # this is slow, we can drop non relevant data
        data = self._get_volcano_sub_data('ALL')
        data['annotation'] = ['' for x in range(len(data))]

        self._volcano_plot(data, title='all drugs all features')

    def _get_fdr_from_pvalue_interp(self, pvalue):
        """Here, FDR are computed using an interpolation"""
        pvalue += 1e-15
        qvals = self.df[self._colname_qvalue]
        pvals = self.df[self._colname_pvalue]
        ya = qvals[pvals < pvalue].max()
        yb = qvals[pvals > pvalue].min()
        xa = pvals[pvals < pvalue].max()
        xb = pvals[pvals > pvalue].min()
        dx = xb - xa
        dy = yb - ya
        yc = ya + dy * (pvalue - xa) / dx
        return yc

    def _get_pvalue_from_fdr(self, fdr):
        """Get pvalue for a given FDR threshold

        This is equivalent to v17 of the R version but is not very precise
        we should use _get_pvalue_from_fdr_interp instead but needs to be
        tested.

        """
        qvals = self.df[self._colname_qvalue]
        pvals = self.df[self._colname_pvalue]
        if isinstance(fdr, list):
            pvalues = [pvals[qvals < this].max() for this in fdr]
            return pvalues
        else:
            return pvals[qvals < fdr].max()

    def _get_pvalue_from_fdr_interp(self, fdr):
        # same as get_pvalue_from_fdr but with a linear inerpolation
        fdr += 1e-15
        qvals = self.df[self._colname_qvalue]
        pvals = self.df[self._colname_pvalue]
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
        colname = self._colname_N_feature_pos
        minN = self.df[colname].min()
        maxN = self.df[colname].max()
        pvalues = self._get_pvalue_from_fdr(self.settings.FDR_threshold)
        return {'minN': minN, 'maxN': maxN,
                'pvalues': (self.settings.FDR_threshold, pvalues)}

    def _get_volcano_sub_data(self, mode, target=None):
        # Return data needed for each plot
        # TODO could be simplified but works for now

        # groups created in the constructor once for all
        if mode == self._colname_drugid:
            subdf = self.df.ix[self.groups_by_drugs[target]]
            texts = subdf[self._colname_feature]
        elif mode == 'FEATURE':
            subdf = self.df.ix[self.groups_by_features[target]]
            texts = subdf[self._colname_drugid]
        elif mode == 'ALL':
            # nothing to do, get all data
            subdf = self.df
            texts = subdf[self._colname_feature] # TODO + drug
        else:
            raise ValueError("mode parameter must be in [FEATURE, %s, ALL]" %
                    (self._colname_drugid))

        # replaced by groups created in the constructor
        #subdf = self.df[self.df[mode] == target]
        deltas = subdf[self._colname_deltas]
        effects = subdf[self._colname_effect_size]
        signed_effects = list(np.sign(deltas) * effects)
        qvals = list(subdf[self._colname_qvalue])
        pvals = list(subdf[self._colname_pvalue])
        #assocs = list(subdf['ASSOC_ID'])

        colors = []
        annotations = []

        data = pd.DataFrame(index=range(len(qvals)))
        data['pvalue'] = pvals
        data['signed_effect'] = signed_effects
        data['Feature'] = list(subdf[self._colname_feature])
        data['Drug'] = list(subdf[self._colname_drugid])
        data['text'] = texts.values
        #data['Assoc'] = assocs
        ## !! here, we need to use .values since the pandas dataframe
        # index goes from 1 to N but the origignal indices in subdf
        # may not be from 1 to N but random between 1 and M>>N
        data['FDR'] = subdf[self._colname_qvalue].values
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
        colname = self._colname_N_feature_pos
        markersize = subdf[colname] / subdf[colname].max()
        markersize = list(markersize * 800)
        markersize = [x if x > 80 else 80 for x in markersize]

        data['color'] = colors
        data['annotation'] = annotations
        data['markersize'] = markersize
        return data

    def volcano_plot_one_feature(self, feature):
        """Volcano plot for one feature (all drugs)

        :param feature: a valid feature name to be found in the results
        """
        assert feature in self.features, 'unknown feature name'
        # FEATURE is the mode's name, not a column's name
        data = self._get_volcano_sub_data('FEATURE', feature)
        self._volcano_plot(data, title=feature)

    def volcano_plot_one_drug(self, drug_id):
        """Volcano plot for one drug (all genomic features)

        :param drug_id: a valid drug identifier to be found in the results
        """
        assert drug_id in self.drugs, 'unknown drug name'
        data = self._get_volcano_sub_data(self._colname_drugid, drug_id)
        self._volcano_plot(data, title=drug_id)

    def _volcano_plot(self, data, title=''):
        """Main volcano plot function called by other methods
        such as volcano_plot_all"""
        # This functio is a bit complicated because it does create a few tricky
        # plots

        # It creates a volcano plot, which is the easy part
        # Then, it creates tooltips for the user interface in an IPython
        # shell using a callback to 'onpick' function coded here below
        # finally, it creates a Javascript connection using mpld3 that
        # will allow the creation of a JS version of the plot.

        # !! There is a memory leak in this function due to matplotlib
        # This is not easy to track down.

        # You have to call clf() to make sure the content is erase.
        # One reason for the memory leak is that it is called in the
        # Report to loop over all drugs and then all featuers.
        # To see the memory leak, you will need to call the
        # volcano_plot_all_drugs function (or volcano_plot_all_features).
        colors = list(data['color'].values)
        pvalues = data['pvalue'].values
        signed_effects = data['signed_effect'].values
        markersize = data['markersize'].values

        Y = -np.log10(list(pvalues)) # should be cast to list ?

        num = 1
        #pylab.close(num)
        fig = pylab.figure(num=1)
        fig.clf()
        ax = fig.add_subplot(111)
        ax.set_axis_bgcolor('#EEEEEE')
        ax.cla()
        X = [easydev.precision(x, digit=2) for x in signed_effects]
        Y = [easydev.precision(y, digit=2) for y in Y]

        # Using scatter() is slow as compared to plot()
        # However, plot cannot take different sizes/colors
        scatter = ax.scatter(X, Y, s=markersize,
                alpha=0.3, c=colors, linewidth=1, picker=True)
        scatter.set_zorder(11)

        m = abs(signed_effects.min())
        M = abs(signed_effects.max())
        pylab.xlabel("Signed effect size", fontsize=self.settings.fontsize)
        pylab.ylabel('-log10(pvalues)', fontsize=self.settings.fontsize)
        l = max([m, M]) * 1.1
        pylab.xlim([-l, l])
        ax.grid(color='white', linestyle='solid')

        # some aliases
        fdr = self.settings.FDR_threshold
        fdrs = sorted(self.settings.volcano_additional_FDR_lines)
        fdrs = fdrs[::-1] # reverse sorting

        styles = ['--', ':', '-.']
        if self.settings.volcano_FDR_interpolation is True:
            get_pvalue_from_fdr = self._get_pvalue_from_fdr_interp
        else:
            get_pvalue_from_fdr = self._get_pvalue_from_fdr

        pvalue = get_pvalue_from_fdr(fdr)
        ax.axhline(-np.log10(pvalue), linestyle='--', lw=2,
            color='red', alpha=1, label="FDR %s " %  fdr + " \%")

        for i, this in enumerate(fdrs):
            if this < self.df[self._colname_qvalue].min() or\
                this > self.df[self._colname_qvalue].max():
                    continue
            pvalue = get_pvalue_from_fdr(this)
            ax.axhline(-np.log10(pvalue), linestyle=styles[i],
                color='red', alpha=1, label="FDR %s " % this +" \%")

        pylab.ylim([0, pylab.ylim()[1]*1.2]) # times 1.2 to put the legend

        ax.axvline(0, color='gray', alpha=0.8, lw=2)
        axl = pylab.legend(loc='best')
        axl.set_zorder(10) # in case there is a circle behind the legend.

        #self.ax = ax
        #self.axx = ax.twinx()
        #self.common_ticks = ax.get_yticks()
        #self.common_ylim = ax.get_ylim()
        #pvals = self.df[self._colname_pvalue]
        #y1 = pvals.min()
        #y2 = pvals.max()
        #fdr1 = self._get_fdr_from_pvalue_interp(y1)
        #fdr2 = self._get_fdr_from_pvalue_interp(y2-2e-15) # make sure it exists
        #self.axx.set_ylim([fdr2, fdr1])
        #self.axx.set_ylabel('FDR \%', fontsize=self.settings.fontsize)

        # For the static version
        title_handler = pylab.title("%s" % title.replace("_","  "),
                fontsize=self.settings.fontsize/1.2)
        labels = []

        # This code allows the ipython user to click on the matplotlib figure
        # to get information about the drug and feature of a given circles.
        def onpick(event):
            ind = event.ind[0]
            try:
                title = data.ix[ind]['Drug'] + " / " + data.ix[ind].Feature
                title += "\nFDR=" + str(data.ix[ind]['FDR'])
                title_handler.set_text(title.replace("_","  "))
            except:
                print('Failed to create new title on click')
            print(data.ix[ind].T)
            fig.canvas.draw()

        # keep track on the id for further memory release
        # For more info search for "matplotlib memory leak mpl_connect"
        self.cid = fig.canvas.mpl_connect('pick_event', onpick)

        # for the JS version
        # TODO: for the first 1 to 2000 entries ?
        labels = []


        self.data = data
        for i, row in data[['Drug', 'Feature', 'FDR']].iterrows():

            template = """
<table border="1" class="dataframe">
  <tbody>
    <tr>
      <th>Drug</th>
      <td>%(Drug)s</td>
    </tr>
    <tr>
      <th>Feature</th>
      <td>%(Feature)s</td>
    </tr>
    <tr>
      <th>FDR</th>
      <td>%(FDR)s</td>
    </tr>
  </tbody>
</table>""" % row.to_dict()
            labels.append(template)

            # this is more elegant but slower
            #label = row.to_frame()
            #label.columns = ['Row {0}'.format(i)]
            #labels.append(str(label.to_html(header=False)))
        css = """
        svg.mpld3-figure { border: 2px black solid;margin:10px;}
        table{  font-size:0.8em;  }
        th {  color: #ffffff;  background-color: #aaaaaa;  }
        td { color: blue; background-color: #cccccc; }"""

        try:
            tooltip = mpld3.plugins.PointHTMLTooltip(scatter, labels=labels,
                css=css)
            mpld3.plugins.connect(fig, tooltip)
        except:
            pass
        self.scatter = scatter
        self.current_fig = fig
        # not sure is this is required. could be a memory leak here
        import gc
        gc.collect()

    def mpld3_to_html(self):
        """This require to call a plotting figure before hand"""
        from gdsctools import gdsctools_data
        # This copy the full path and therefore HTML cannot
        # be moved in another directory. to be fixed.
        js_path1 = gdsctools_data('d3.v3.min.js', where='javascript')
        js_path2 = gdsctools_data('mpld3.v0.2.js', where='javascript')
        try:
            # mpld3 is great but there are a couple of issues
            # 1 - legend zorder is not used so dots may be below the legend,
            #     hence we set the framealpha =0.5
            # 2 - % character even though there well interpreted in matploltib
            #     using \%, they are not once parsed by mpld3. So, here
            #     we remove the \ character
            axl = pylab.legend(loc='best', framealpha=0.8, borderpad=1)
            axl.set_zorder(10) # in case there is a circle behind the legend.
            texts = [this.get_text() for this in axl.get_texts()]

            for i, text in enumerate(texts):
                text = text.replace("\\%", "%")
                text += "  "
                axl.get_texts()[i].set_text(text)
            import mpld3
            js_path1 = "https://raw.githubusercontent.com/CancerRxGene/" + \
                    "gdsctools/master/share/data/javascript/d3.v3.min.js"
            js_path2 = "https://raw.githubusercontent.com/CancerRxGene/" + \
                    "gdsctools/master/share/data/javascript/mpld3.v0.2.js"
            htmljs = mpld3.fig_to_html(self.current_fig)
            #                d3_url=js_path1,
            #                mpld3_url=js_path2)
        except:
            htmljs = ""
        return """<div class="jsimage"> """ + htmljs + "</div>"

    def savefig_and_js(self, filename, size_inches=(10, 10)):
        # Save the PNG first. The savefig automatically set the size
        # to a defined set and back to original figsize.
        self.figtools.savefig(filename + '.png', size_inches=size_inches)

        # now the javascript. 
        fig = self.current_fig
        oldsize = fig.get_size_inches()
        fig.set_size_inches(size_inches)
        htmljs = self.mpld3_to_html()
        fh = open(self.settings.directory + os.sep + filename + ".html", "w")
        fh.write(htmljs)
        fh.close()
        fig.set_size_inches(*oldsize)


