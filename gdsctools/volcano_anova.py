import pandas as pd
import pylab
import numpy as np
import easydev

from easydev import Progress, AttrDict
from cno.misc.profiler import do_profile


__all__ = ['VolcanoANOVA']

from gdsctools.tools import Savefig

class VolcanoANOVA(Savefig):
    def __init__(self, data, sep="\t"):
        super(VolcanoANOVA, self).__init__()

        if isinstance(data, str):
            self.df = pd.read_csv(data, sep=sep)
        else:
            self.df = data.copy()

        # this is redundant coul reuse the input ??
        self.settings = {
            'pval_threshold': np.inf,
            'fdr_threshold': 25,
            'effect_threshold': 0,
            'fontsize': 20,
            'analysis_type': 'PANCAN',
            'savefig': True,
            'directory': '.'}
        self.settings = AttrDict(**self.settings)
        self.drugs = set(self.df['Drug id']) # set on values
        self.varname_pvalue = 'FEATURE_ANOVA_pval'
        self.varname_qvalue = 'ANOVA FEATURE FDR %'

        self.stats = self._get_volcano_global_data()

        # intensive calls made once for all
        self.groups_by_drugs = self.df.groupby('Drug id').groups
        self.groups_by_features = self.df.groupby('FEATURE').groups

    def volcano_plot_all_drugs(self):
        """Volcano plot for each drug and savefig

        :param df: output of :meth:`anova_all`
        """
        drugs = list(self.df['Drug id'].unique())
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

        self.volcano_plot(data, title='all')
        if self.settings.savefig is True:
            self.savefig("volcano_all.png")

    def _get_volcano_global_data(self):

        # using all data
        minN = self.df['N_FEATURE_pos'].min()
        maxN = self.df['N_FEATURE_pos'].max()
        qvals = self.df[self.varname_qvalue]
        pvals = self.df[self.varname_pvalue]
        fdrlim = pvals[qvals < self.settings.fdr_threshold].max()
        fdrlim1 = pvals[qvals < 10].max()
        fdrlim2 = pvals[qvals < 1].max()
        fdrlim3 = pvals[qvals < 0.01].max()

        return {'minN': minN, 'maxN': maxN,
                'fdrs': {
                    self.settings.fdr_threshold: fdrlim,
                    0.01: fdrlim3,
                    1: fdrlim2,
                    10: fdrlim1}
                }

    def volcano_plot_one_feature(self, feature):
        """Volcano plot for one feature (all drugs)


        """
        data = self._get_volcano_sub_data('FEATURE', feature)
        self.volcano_plot(data, title=feature)

    #@do_profile()
    def _get_volcano_sub_data(self, mode, target=None):
        # using data related to the given drug

        # groups created in the constructor once for all
        if mode == 'Drug id':
            subdf = self.df.ix[self.groups_by_drugs[target]]
            texts = subdf['FEATURE']
        elif mode == 'FEATURE':
            subdf = self.df.ix[self.groups_by_features[target]]
            texts = subdf['Drug id']
        elif mode == 'ALL':
            # nothing to do, get all data
            subdf = self.df
            texts = subdf['FEATURE'] # TODO + drug
        else:
            raise ValueError("mode parameter must be in [FEATURE, Drug id, ALL]")

        # replaced by groups created in the constructor
        #subdf = self.df[self.df[mode] == target]

        deltas = subdf['FEATURE_deltaMEAN_IC50']
        effects = subdf['FEATURE_IC50_effect_size']
        assoc_ids = list(subdf.index)
        signed_effects = list(np.sign(deltas) * effects)
        qvals = list(subdf[self.varname_qvalue])
        pvals = list(subdf[self.varname_pvalue])

        colors = []
        annotations = []

        data = pd.DataFrame(index=range(len(qvals)))
        data['pvalue'] = pvals
        data['signed_effect'] = signed_effects
        data['feature'] = list(subdf['FEATURE'])
        data['drug'] = list(subdf['Drug id'])
        data['text'] = texts.values
        annotations = []

        if self.settings.analysis_type == 'PANCAN':
            for sign, qval, pval in zip(signed_effects, qvals, pvals):
                if sign <= -self.settings.effect_threshold and \
                        qval <= self.settings.fdr_threshold:
                    colors.append('green')
                    annotations.append(True)
                elif sign >= self.settings.effect_threshold and \
                        qval <= self.settings.fdr_threshold:
                    colors.append('red')
                    annotations.append(True)
                else:
                    colors.append('black')
                    annotations.append(False)
        else:
            for delta, qval, pval in zip(deltas, qvals, pvals):
                if pval <= self.settings.pvalue_threshold and \
                        qval <= self.settings.fdr_threshold and delta<0:
                    colors.append('green')
                    annotations.append(True)
                elif pval <= self.settings.pvalue_threshold and \
                        qval <= self.settings.fdr_threshold and delta>0:
                    colors.append('red')
                    annotations.append(True)
                else:
                    colors.append('black')
                    annotations.append(False)


        # here we normalise wrt the drug. In R code, normalised
        # my max across all data (minN, maxN)
        markersize = subdf['N_FEATURE_pos'] / subdf['N_FEATURE_pos'].max()
        markersize = list(markersize*800)
        markersize = [x if x>50 else 50 for x in markersize]

        data['color'] = colors
        data['annotation'] = annotations
        data['markersize'] = markersize
        return data

    #@do_profile()
    def volcano_plot_one_drug(self, drug_id):
        """Volcano plot for one drug and all genomic features

        :param df: output of :meth:`anova_all`

        """
        assert drug_id in self.drugs, 'unknown drug name'
        # needs to run :meth:`anova_all` first
        # using all data, get the FDR limits
        data = self._get_volcano_sub_data('Drug id', drug_id)
        self.volcano_plot(data, title=drug_id)

    #@do_profile()
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
        fig = pylab.figure(num=1)
        pylab.clf()
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

        #print(fdrlim, fdrlim1, fdrlim2, fdrlim3)
        fdrlim = self.stats['fdrs'][self.settings.fdr_threshold]
        fdrlim1 = self.stats['fdrs'][10]
        fdrlim2 = self.stats['fdrs'][1]
        fdrlim3 = self.stats['fdrs'][0.01]

        ax.axhline(-np.log10(fdrlim), linestyle='--',
            color='red', alpha=1,
            label="FDR %s " %  self.settings.fdr_threshold + " \%")
        ax.axhline(-np.log10(fdrlim1), linestyle='-.',
            color='red', alpha=1, label="FDR 10 \%")
        ax.axhline(-np.log10(fdrlim2), linestyle=':',
            color='red', alpha=1, label="FDR 1 \%")
        ax.axhline(-np.log10(fdrlim3), linestyle='--',
            color='black', alpha=1, label="FDR 0.01 \%")

        pylab.ylim([0, pylab.ylim()[1]*1.2]) # times 1.2 to put the legend

        ax.axvline(0, color='gray', alpha=0.5)
        axl = pylab.legend(loc='upper left')
        axl.set_zorder(-1) # in case there is a circle behind the legend.

        # For the static version
        pylab.title("%s" % title.replace("_","\_"))
        labels = []
        for index, row in data.iterrows():
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
        """fig.canvas.mpl_connect('pick_event', onpick)

        # for the JS version
        import mpld3
        labels = []
        for i, row in data[['drug','feature']].iterrows():
            label = row.to_frame()
            label.columns = ['Row {0}'.format(i)]
            # .to_html() is unicode; so make leading 'u' go away with str()
            labels.append(str(label.to_html(header=False)))

        tooltip = mpld3.plugins.PointHTMLTooltip(scatter, labels=labels)
        mpld3.plugins.connect(fig, tooltip)
        """
        #mpld3.display()

        self.scatter = scatter
        self.current_fig = fig
        #fh = open('test.html', 'w'); mpld3.save_html(v.current_fig, fh);
        #fh.close()
        

