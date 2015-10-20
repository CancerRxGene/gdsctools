import pandas as pd
import pylab
import numpy as np


from easydev import Progress, AttrDict


__all__ = ['VolcanoANOVA']


class VolcanoANOVA(object):
    def __init__(self, data, sep="\t"):
        if isinstance(data, str):
            self.df = pd.read_csv(data, sep=sep)
        else:
            self.df = data.copy()
        self.settings = {
            'pval_threshold': np.inf,
            'FDR_threshold': 20,
            'effect_threshold': 0,
            'fontsize': 20,
            'analysis_type': 'PANCAN',
            'savefig': True}
        self.settings = AttrDict(**self.settings)

        self.varname_pvalue = 'FEATURE_ANOVA_pval'
        self.varname_qvalue = 'ANOVA FEATURE FDR %'

    def volcano_plot_all_drugs(self):
        """Volcano plot for each drug and savefig

        :param df: output of :meth:`anova_all`
        """
        drugs = list(self.df['Drug id'].unique())
        pb = Progress(len(drugs), 1)
        for i, drug in enumerate(drugs):
            print drug
            self.volcano_plot_one_drug(drug)
            if self.settings.savefig is True:
                pylab.savefig("volcano_%s.png" % drug)
            pb.animate(i+1)

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
                pylab.savefig("volcano_%s.png" % feature)
            pb.animate(i+1)

    def _get_volcano_global_data(self):

        # using all data
        minN = self.df['N_FEATURE_pos'].min()
        maxN = self.df['N_FEATURE_pos'].max()
        qvals = self.df[self.varname_qvalue]
        pvals = self.df[self.varname_pvalue]
        fdrlim = pvals[qvals < self.settings.FDR_threshold].max()
        fdrlim1 = pvals[qvals < 10].max()
        fdrlim2 = pvals[qvals < 1].max()
        fdrlim3 = pvals[qvals < 0.01].max()

        return {'minN': minN, 'maxN': maxN,
                'fdrs': {
                    self.settings.FDR_threshold: fdrlim,
                    0.01: fdrlim3,
                    1: fdrlim2,
                    10: fdrlim1}
                }

    def volcano_plot_one_feature(self, feature):
        """Volcano plot for one feature (all drugs)


        """
        stats = self._get_volcano_global_data()
        data = self._get_volcano_sub_data('FEATURE', feature)

        self.volcano_plot(data.signed_effects, data.pvals, data.markersize,
                data.colors, data.annotations, stats,
                title=feature)

    def _get_volcano_sub_data(self, mode, target):
        # using data related to the given drug

        if mode == 'Drug id':
            other = 'FEATURE'
        elif mode == 'FEATURE':
            other = 'Drug id'
        else:
            raise ValueError("mode parameter must be 'FEATURE' or 'Drug id'")

        subdf = self.df[self.df[mode] == target].copy()

        deltas = subdf['FEATURE_deltaMEAN_IC50']
        effects = subdf['FEATURE_IC50_effect_size']
        signed_effects = np.sign(deltas) * effects

        qvals = list(subdf[self.varname_qvalue])
        pvals = list(subdf[self.varname_pvalue])
        features = subdf[other]

        colors = []
        annotations = []
        labels = features
        if self.settings.analysis_type == 'PANCAN':
            for sign, qval, pval, label in zip(signed_effects, qvals,
                    pvals, labels):
                if sign <= -self.settings.effect_threshold and \
                        qval <= self.settings.FDR_threshold:
                    colors.append('green')
                    annotations.append((sign, pval, label))
                elif sign >= self.settings.effect_threshold and \
                        qval <= self.settings.FDR_threshold:
                    colors.append('red')
                    annotations.append((sign, pval, label))
                else:
                    colors.append('grey')
        else:
            for delta, qval, pval, label in zip(deltas, qvals,
                    pvals, labels):
                if pval <= self.settings.pvalue_threshold and \
                        qval <= self.settings.FDR_threshold and delta<0:
                    colors.append('green')
                    annotations.append((delta, pval, label))
                elif pval <= self.settings.pvalue_threshold and \
                        qval <= self.settings.FDR_threshold and delta>0:
                    colors.append('red')
                    annotations.append((delta, pval, label))
                else:
                    colors.append('grey')

        # here we normalise wrt the drug. In R code, normalised
        # my max across all data (minN, maxN)
        markersize = subdf['N_FEATURE_pos'] / subdf['N_FEATURE_pos'].max()
        markersize = list(markersize*800)
        markersize = [x if x>50 else 50 for x in markersize]

        dd = {'colors':colors, 'annotations':annotations,
                'signed_effects':signed_effects,
                'markersize':markersize, 'qvals':qvals, 'pvals':pvals}
        dd = AttrDict(**dd)
        return dd


    def volcano_plot_one_drug(self,drug_id):
        """Volcano plot for one drug and all genomic features

        :param df: output of :meth:`anova_all`

        """
        assert drug_id in self.df['Drug id'].values, 'unknown drug name'
        # needs to run :meth:`anova_all` first
        # using all data, get the FDR limits
        stats = self._get_volcano_global_data()
        data = self._get_volcano_sub_data('Drug id', drug_id)

        self.volcano_plot(data.signed_effects, data.pvals, data.markersize,
                data.colors, data.annotations, stats,
                title=drug_id)

    def volcano_plot(self, signed_effects, pvals, markersize,
            colors, annotations, stats, title=''):
        """Plots signed effects versus pvalues

        THis is for internal usage but maybe useful elsewhere.

        """
        Y = -np.log10(list(pvals)) # somehow should be cast to list ?

        pylab.clf()

        pylab.scatter(list(signed_effects), Y, s=markersize,
                alpha=0.4, c=colors,
                linewidth=0)

        m = abs(signed_effects.min())
        M = abs(signed_effects.max())
        pylab.xlabel("Signed effect size", fontsize=self.settings.fontsize)
        pylab.ylabel('-log10(pvalues)', fontsize=self.settings.fontsize)
        l = max([m, M]) * 1.1
        pylab.xlim([-l, l])

        #print(fdrlim, fdrlim1, fdrlim2, fdrlim3)
        fdrlim = stats['fdrs'][self.settings.FDR_threshold]
        fdrlim1 = stats['fdrs'][10]
        fdrlim2 = stats['fdrs'][1]
        fdrlim3 = stats['fdrs'][0.01]
        pylab.axhline(-np.log10(fdrlim), linestyle='--',
            color='gray', alpha=1,
            label="FDR %s " %  self.settings.FDR_threshold + " \%")
        pylab.axhline(-np.log10(fdrlim1), linestyle='-.',
            color='gray', alpha=1, label="FDR 10 \%")
        pylab.axhline(-np.log10(fdrlim2), linestyle=':',
            color='gray', alpha=1, label="FDR 1 \%")
        pylab.axhline(-np.log10(fdrlim3), linestyle='--',
            color='black', alpha=1, label="FDR 0.01 \%")

        pylab.ylim([0, pylab.ylim()[1]*1.2]) # times 1.2 to put the legend

        pylab.axvline(0, color='gray', alpha=0.5)
        ax = pylab.legend(loc='upper left')
        ax.set_zorder(-1) # in case there is a circle behind the legend.
        pylab.title("%s" % title.replace("_","\_"))
        for this in annotations:
            x,y,text = this
            pylab.text(x,-pylab.log10(y),text.replace("_", "\_"))

