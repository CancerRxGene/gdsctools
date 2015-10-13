import pandas as pd
from easydev import Progress, AttrDict
import pylab





def boxswarm(data, names=None, vert=True, widths=0.5, **kwargs):
    """Plot boxplot with all points as circles



    See :class:`BoxSwarm`documentation for details


    """
    b = BoxSwarm(data, names=names)
    b.plot(vert=vert, widths=widths, **kwargs)
    return b


class BoxSwarm(object):
    def __init__(self, data, names=None):
        """


        :param: a list of list (not same size) or a dictionary of lists



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

                self.data = dict([(name, d) for name, d in zip(self.names, data)])

        self.ylabel = ''
        self.xlabel = ''
        self.fontsize = 20
        self.colors = ['lightgrey', 'blue']
        self.hold = False
        self.title = ''

    def beeswarm(self, data, ratio=2.):
        N = len(data)
        newdata = data + (pylab.random(N) - 0.5)/float(ratio)
        return newdata


    def plot(self, vert=True, widths=0.5, **kwargs):
        """


        """
        self.widths = widths
        if self.hold is False:
            pylab.clf()

        ordered_data = [self.data[key] for key in self.names]

        for i, vector in enumerate(ordered_data):
            N = len(vector)

            color = self.colors[i%len(self.colors)]
            if vert is True:
                X, Y = self.beeswarm([i+1] * N), vector
            else:
                X, Y = vector, self.beeswarm([i+1] * N)
            pylab.plot(X, Y,
                'o', markersize=15, markerfacecolor=color,
                markeredgewidth=0, alpha=0.5)

        means = [pylab.mean(data) for data in ordered_data]
        stds = [pylab.std(data) for data in ordered_data]

        d = pylab.boxplot(ordered_data, widths=self.widths,
                vert=vert, patch_artist=True)

        # for further tuning if needed.
        self.tuning = d
        for i, this in enumerate(means):
            if vert is True:
                x1 = (i+1) - widths/2. / 1.5
                x2 = (i+1) + widths/2. / 1.5
                X = pylab.array([x1,x2])
                pylab.plot(X, [this, this], lw=2, color='r')
            else:
                y1 = (i+1) - widths/2. / 1.5
                y2 = (i+1) + widths/2. / 1.5
                Y = pylab.array([y1,y2])
                pylab.plot([this, this], Y, lw=2, color='r')

        for i, this in enumerate(d['boxes']):
            this.set_color('k')
            this.set_linewidth(2)
            color = self.colors[i%len(self.colors)]
            this.set_facecolor(color)
            this.set_alpha(0.4) # 0.4 is less than the alpha of the dots to ...
            # ... so as to see the dots inside the boxes
            this.set_zorder(10) # this moves the box on top of all dots
        for this in d['caps']:
            this.set_linewidth(2)
        for this in d['whiskers']:
            this.set_linewidth(2)
        for this in d['medians']:
            this.set_linewidth(2)
        for this in d['fliers']:
            #delete ?
            this.set_linewidth(0)

        N = len(self.names)
        if vert is True:
            pylab.ylabel(self.ylabel, fontsize=self.fontsize)
            pylab.xticks(range(1, N+1), self.names, fontsize=self.fontsize,
                rotation=90)
            pylab.ylabel(self.xlabel, fontsize=self.fontsize)
            pylab.yticks(pylab.yticks()[0], fontsize=self.fontsize)
            m = min([min(this) for this in self.data.values()])
            M = max([max(this) for this in self.data.values()])
            R = (M-m)* 0.05
            pylab.ylim([m-R,M+R])
        else:
            pylab.xlabel(self.xlabel, fontsize=self.fontsize)
            pylab.yticks(range(1,N+1), self.names, fontsize=self.fontsize,
                rotation=00)
            pylab.ylabel(self.ylabel, fontsize=self.fontsize)
            pylab.xticks(pylab.xticks()[0], fontsize=self.fontsize)
            m = min([min(this) for this in self.data.values()])
            M = max([max(this) for this in self.data.values()])
            R = (M-m) * 0.05
            pylab.xlim([m-R,M+R])
        pylab.title(self.title, fontsize=self.fontsize)
        pylab.grid()
        pylab.tight_layout()

