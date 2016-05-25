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
import pylab
import numpy as np


__all__ = ['boxswarm', 'BoxSwarm']


def boxswarm(data, names=None, vert=True, widths=0.5, **kwargs):
    """Plot boxplot with all points as circles.

    This function is a wrapper of :class:`BoxSwarm`

    :param data: a dataframe. Each column is a data set from which a boxplot
        is created.
    :param names:
    :param vert: orientation of the boxplots
    :param widths: widths of the boxes
    :param kargs: any argument accepted by :class:`BoxSwarm`

    See :class:`BoxSwarm` documentation for details

    """
    b = BoxSwarm(data, names=names, **kwargs)
    b.plot(vert=vert, widths=widths, **kwargs)
    return b


class BoxSwarm(object):
    """Simple beeswarm plot (boxplot + dots for each data point)


    .. plot::
        :include-source:
        :width: 80%

        from pylab import randn
        from gdsctools.boxswarm import BoxSwarm
        b = BoxSwarm({'a':randn(100), 'b':randn(20)+2})
        b.plot(vert=False)


    .. note:: could use pybeeswarm, which is a proper implementation
        of beeswarm.

    """
    def __init__(self, data, names=None, fontsize=20, hold=False,
            title='', lw=2, colors=['lightgrey', 'blue']):
        """.. rubric:: Constructor

        :param: a list of list (not same size) or a dictionary of lists

        :param data:
        :param names:
        :param fontsize:
        :param hold:
        :param title:
        :param lw: width of lines
        :param colors: loop over the list of colors provided to fill boxplots
        :param **kargs:


        """
        # if a list, we create a dictionary internally
        try:
            # a dataframe ?
            self.data = data.to_dict('list')
            if names is None: # no order, let us sort alphabetically
                self.names = sorted(self.data.keys())
        except:
            # a dictionary ? nothing to do
            if isinstance(data, dict):
                self.data = data
                if names is None: # no order, let us sort alphabetically
                    self.names = sorted(self.data.keys())
            else:
                # probably a list of list or arrays or array without names
                if names is None:
                    self.names = range(0, len(data))
                else:
                    assert len(names) == len(data)
                    self.names = names

                self.data = dict([(name, d)
                    for name, d in zip(self.names, data)])

        self.ylabel = ''
        self.xlabel = ''
        self.fontsize = fontsize
        self.colors = colors
        self.hold = hold
        self.title = title
        self.lw = lw
        self.markersize = 6

    def beeswarm(self, data, position, ratio=2.):
        r"""Naive plotting of the data points

        We assume gaussian distribution so we expect fewers dots
        far from the mean/median. We'd like those dots to be close to the
        axes. conversely, we expect lots of dots centered around the mean, in
        which case, we'd like them to be spread in the box. We uniformly
        distribute position using

        .. math::

            X = X + \dfrac{ U()-0.5 }{ratio} \times factor

        but the factor is based on an arctan function:

        .. math::

            factor = 1 - \arctan( \dfrac{X - \mu }{\pi/2})

        The farther the data is from the mean :math:`\mu`,
        the closest it is to the axes that goes through the box.

        """
        N = len(data)
        m = np.median(data)
        sd = np.std(data)
        # arctan function to have a tapering window
        factor = 1. - np.abs(np.arctan((data-m)/sd)/1.570796)  # pi/2

        newdata = position + (pylab.random(N) - 0.5)/float(ratio) * factor
        return newdata

    def plot(self, vert=True, alpha=0.4, widths=0.5, **kwargs):
        """Plot the boxplots and dots


        """
        self.widths = widths
        if self.hold is False:
            pylab.clf()

        ordered_data = [self.data[key] for key in self.names]

        for i, vector in enumerate(ordered_data):
            N = len(vector)

            color = self.colors[i%len(self.colors)]
            if vert is True:
                X, Y = self.beeswarm(vector, i+1), vector
            else:
                X, Y = vector, self.beeswarm(vector, i+1)

            pylab.plot(X, Y,
                'o', markersize=self.markersize, markerfacecolor=color,
                markeredgewidth=1, alpha=alpha)

        #show means but not outliers
        try:
            d = pylab.boxplot(ordered_data, widths=self.widths,
               vert=vert, patch_artist=True,
                positions=range(1, len(ordered_data)+1),
            showmeans=True, showfliers=False)
        except:
            # ReadTheDocs uses matplotlib 1.3.1 for now, so
            # need this without showmeans parameter
            d = pylab.boxplot(ordered_data, widths=self.widths,
               vert=vert, patch_artist=True,
                positions=range(1, len(ordered_data)+1))

        # for further tuning if needed.
        self.tuning = d
        # This is now in matplotlib 1.4.3 (dots instead of lines
        # though)

        # additional line for the 1 std
        means = [pylab.mean(data) for data in ordered_data]
        stds = [pylab.std(data) for data in ordered_data]
        for i, this in enumerate(means):
            if vert is True:
                x1 = (i+1) - widths/2. / 1.5
                x2 = (i+1) + widths/2. / 1.5
                X = pylab.array([x1, x2])
                y = this + stds[i]
                pylab.plot(X, [y, y], lw=2, color='purple')
                y = this - stds[i]
                pylab.plot(X, [y, y], lw=2, color='purple')
            else:
                y1 = (i+1) - widths/2. / 1.5
                y2 = (i+1) + widths/2. / 1.5
                Y = pylab.array([y1, y2])
                x = this + stds[i]
                pylab.plot([x, x], Y, lw=2, color='purple')
                x = this - stds[i]
                pylab.plot([x, x], Y, lw=2, color='purple')

        for i, this in enumerate(d['boxes']):
            this.set_color('k')
            this.set_linewidth(self.lw)
            color = self.colors[i%len(self.colors)]
            this.set_facecolor(color)
            this.set_alpha(0.3) # 0.4 is less than the alpha of the dots to ...
            # ... so as to see the dots inside the boxes
            this.set_zorder(10) # this moves the box on top of all dots
        for this in d['caps']:
            this.set_linewidth(self.lw)
        for this in d['whiskers']:
            this.set_linewidth(self.lw)
        for this in d['medians']:
            this.set_linewidth(self.lw)

        # we will extend the limits by 5%
        m = min([min(this) for this in self.data.values()])
        M = max([max(this) for this in self.data.values()])
        extend = 0.05
        R = (M-m) * extend
        X, Y = range(1, len(self.names)+1), self.names
        Y = [y.replace("_", " ") for y in Y]
        if vert is True:
            pylab.ylabel(self.ylabel, fontsize=self.fontsize)
            pylab.xticks(X, Y, fontsize=self.fontsize, rotation=90)
            pylab.ylabel(self.xlabel, fontsize=self.fontsize)
            pylab.yticks(pylab.yticks()[0], fontsize=self.fontsize)
            pylab.ylim([m-R, M+R])
        else:
            pylab.xlabel(self.xlabel, fontsize=self.fontsize)
            if len(X) > 20:
                pylab.yticks(X, Y, fontsize=self.fontsize/1.6, rotation=00)
            else:
                pylab.yticks(X, Y, fontsize=self.fontsize, rotation=00)

            pylab.ylabel(self.ylabel, fontsize=self.fontsize)
            pylab.xticks(pylab.xticks()[0], fontsize=self.fontsize)
            pylab.xlim([m-R, M+R])

        pylab.title(self.title, fontsize=self.fontsize*1.25)
        pylab.grid()
        try:pylab.tight_layout()
        except:pass

        return pylab.gca()
