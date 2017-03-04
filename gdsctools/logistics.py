# -*- python -*-
# -*- coding utf-8 -*-
#
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
"""Sets of miscellaneous tools"""
import math
import warnings

import numpy as np
import pylab

__all__ = ['Logistic', 'LogisticMatchedFiltering']


def nextpower(x):
    """Returns next power of 2 times 2

    ::

        >>> nextpower(9)
        32
        >>> nextpower(16)
        32
        >>> nextpower(17)
        64

    """
    return int(2**(math.ceil(pylab.log2(x)) + 1 ))


class Logistic(object):
    """Simple logistic class to see the curve implied by xmid/scale parameters

    .. plot::
        :include-source:

        from gdsctools.logistics import Logistic
        from pylab import legend
        tl = Logistic(2, 1)
        tl.plot()
        tl.scale = 4
        tl.plot(hold=True)
        legend(['scale=1', 'scale=4'])

    The X values are set automatically given the number of data points
    :attr:`N` and the minimum and maximum X values. By default a sensible
    values for the minimun and maximum x values are guessed based on
    :attr:`scale` parameter but one can set the :attr:`xmin` and :attr:`xmax`

    .. plot::
        :include-source:
        :width: 80%

        from gdsctools.logistics import Logistic
        from pylab import legend
        tl = Logistic(2,1)
        tl.plot()
        tl.xmin= -4
        tl.xmax = 4
        tl.N = 40
        tl.plot(hold=True)
        legend(['default', 'user defined'])

    """
    def __init__(self, xmid, scale, Asym=1, N=9, increase=False):
        r""".. rubric:: Constructor

        :param xmid: the first logistic function parameter
        :param scale: the second logistic function parameter
        :param Asym: Amplitude of the function
        :param N: number of data point
        :param increase: increasing or decreasing function.
        :param xmin: starting range of x-values
        :param xmax: ending range of x-values

        if :attr:`increase` is False:

        .. math:: L(x) = \frac{Asym}{1 + exp ((xmid-X)/scale)}

        if :attr:`increase` is True

        .. math:: L(x) = 1 - \frac{Asym}{1 + exp ((xmid-X)/scale)}

        Changing :attr:`xmin`, :attr:`xmax` or :attr:`N` does change
        the content of X. You can change :attr:`X` directly.

        """
        self.xmid = xmid
        self._scale = None
        self.scale = scale
        self.Asym = Asym
        self._N = N
        self.increase = increase
        self._xmin = None
        self._xmax = None
        self._update_x()

    def _get_scale(self):
        return self._scale
    def _set_scale(self, scale):
        self._scale = scale
    scale = property(_get_scale, _set_scale)

    def _update_x(self):
        """Return a sensible linear space of X values """
        dX = abs(self.scale)
        if self.xmin  is None:
            self._xmin = self.xmid - dX * self.N

        if self.xmax  is None:
            self._xmax = self.xmid + dX * self.N

        self._X = np.linspace(self.xmin, self.xmax, self.N)

    def plot(self, hold=False):
        """Plot the logistic function"""
        if hold is False:
            pylab.clf()
        pylab.plot(self.X, self.Y, 'o-')
        pylab.grid(True)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        pylab.axvline(self.xmid, color='r', linestyle='--')

    def _get_xmin(self):
        return self._xmin
    def _set_xmin(self, xmin):
        self._xmin = xmin
        self._update_x()
    xmin = property(_get_xmin, _set_xmin,
            doc="set/get the minimum x range")

    def _get_xmax(self):
        return self._xmax
    def _set_xmax(self, xmax):
        self._xmax = xmax
        self._update_x()
    xmax = property(_get_xmax, _set_xmax,
            doc="set/get the maximum x range")

    def _get_N(self):
        return self._N
    def _set_N(self,N):
        self._N = N
        self._update_x()
    N = property(_get_N, _set_N,
            doc="set/get the number points")

    def _get_x(self):
        return self._X
    def _set_x(self, X):
        self._xmin = min(X)
        self._xmax = max(X)
        self._N = len(X)
        self._X = np.array(X)
        # no need to call _update
    X = property(_get_x, _set_x,
            doc="get/set of the x-values")

    def _get_y(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if self.increase is True:
                Y = self.Asym / (1. + np.exp((self.xmid - self.X)/self.scale))
            else:
                Y = 1. - self.Asym / (1. + np.exp((self.xmid - self.X)/self.scale))
        return Y
    Y = property(_get_y, doc='Getter for Y-values')

    def __str__(self):
        txt = "xmid: %s\n" % self.xmid
        txt += "scale: %s\n" % self.scale
        txt += "xmin: %(xmin)s\n" % {'xmin':self.xmin}
        txt += "xmax: %(xmax)s" % {'xmax':self.xmax}
        return txt


class LogisticMatchedFiltering(object):
    """Experimental class to identify parameters of a noisy Logistic function


    This class implements two methods to identify the parameters (xmid and
    scale) of a logistic function.

    Note One constraint is to define the x-values.

    .. plot::
        :include-source:
        :width: 80%

        from gdsctools import logistics
        import pylab

        mf = logistics.LogisticMatchedFiltering(1,2)
        mf.scan(pylab.linspace(-3,3, 10), pylab.linspace(0,5,10))

    """
    def __init__(self, xmid, scale, N=9):
        """.. rubric:: constructor

        :param xmid:
        :param scale:
        :param N:

        """
        # for book keeping and plotting
        self.signal = Logistic(xmid=xmid, scale=scale, N=N)

        # will be used
        self.data = self.signal.Y

        # will be used for matchefd filtering
        self.template = Logistic(xmid=xmid, scale=scale, N=N)

        self.X = pylab.linspace(-8,8, N)

    def _get_N(self):
        return self.signal.N
    def _set_N(self, N):
        self.signal.N = N
        self.template.N = N
        self.data = self.signal.Y
    N = property(_get_N, _set_N)

    def _get_xmid(self):
        return self.signal.xmid
    def _set_xmid(self, xmid):
        self.signal.xmid = xmid
        self.template.xmid = xmid
        self.data = self.signal.Y
    xmid = property(_get_xmid, _set_xmid,
        doc="get/set the xmid value of the logistic signal")

    def _get_scale(self):
        return self.signal.scale
    def _set_scale(self, scale):
        self.signal.scale = scale
        self.template.scale = scale
        self.data.xmid = xmid
    scale = property(_get_scale, _set_scale,
        doc="get/set the scale value of the logistic signal")

    def _get_X(self):
        return self.signal.X
    def _set_X(self, X):
        self.signal.X = X
        self.template.X = X
        self.data = self.signal.Y
    X = property(_get_X, _set_X,
        doc="get/set the X values the logistic signal")

    def set_noisy_data(self, sigma=0.1):
        self.data = self.signal.Y + sigma * pylab.randn(self.N)
        self.data
        self.data[self.data>1] = 1
        self.data[self.data<0] = 0

    def get_snr(self, template, show=True):
        #losc.ligo.org/tutorial_optimal

        # first take fft of data and template
        # we double vector size to next power of 2
        N = nextpower(self.N)

        data = np.concatenate([self.data, [0] * (N-self.N+1)], axis=0)
        data_fft = np.fft.fft(data)

        #For this to work, we need the template and the data to be the
        # same length. So, we'll zero-pad the template before we take the FFT:
        template = np.concatenate([template, [0] * (N-self.N+1)], axis=0)
        template_fft = np.fft.fft(template)

        # we need an estimate of the noise power in each FFT bin.
        # we can assume white noise.

        # -- Calculate the PSD of the data
        #power_data, freq_psd = plt.psd(data[12*fs:],
        #   Fs=fs, NFFT=fs, visible=False)

        # -- Interpolate to get the PSD values at the needed frequencies
        fs = 1
        datafreq = np.fft.fftfreq(data.size)*fs
        #power_vec = np.interp(datafreq, freq_psd, power_data)
        power_vec  = [1] * (N +1)   ## 32 +f0
        # we multiply the Fourier Space template and data, and divide
        # by the noise power in each frequency bin. Taking the Inverse
        # Fourier Transform (IFFT) of the filter output puts it back
        # in the time domain, so the result will be plotted as a function
        # of time off-set between the template and the data:

        # -- Calculate the matched filter output
        optimal = data_fft * template_fft.conjugate() / power_vec
        optimal_time = 2*np.fft.ifft(optimal)

        # Finally, we can normalize the matched filter output so that
        # we expect a value of 1 at times of just noise. Then, the peak
        # of the matched filter output will tell us the signal-to-noise
        # ratio (SNR) of the signal.

        # -- Normalize the matched filter output
        df = np.abs(datafreq[1] - datafreq[0])
        sigmasq = 2*(template_fft * template_fft.conjugate() / power_vec).sum() * df
        sigma = np.sqrt(np.abs(sigmasq))
        SNR = abs(optimal_time) / (sigma)

        # -- Plot the result
        if show:pylab.plot(SNR)

        return SNR

    def scan(self, xmid_range, scale_range, method='mf', show=True):

        N = len(xmid_range)
        results = np.zeros((N, N))
        coords = None
        best = 0

        self.template.X = self.signal.X

        for i, xmid in enumerate(xmid_range):
            for j, scale in enumerate(scale_range):
                self.template.xmid = xmid
                self.template.scale = scale

                if method == 'corr':
                    M = np.corrcoef(self.template.Y, self.data)[0,1]
                    results[i, j] = M
                elif method == 'mf':
                    M = self.get_snr(self.template.Y, show=False).max()
                    results[i, j] = M

                if M > best:
                    best = M
                    coords = (xmid, scale)

        res = self.optimise(coords)

        Ex = abs(self.xmid-res.x[0])/self.xmid*100
        Ey = abs(self.scale-res.x[1])/self.scale*100

        if show is True:
            self.plot(results, coords, best, xmid_range, scale_range)

            pylab.figure(1)
            pylab.plot(self.xmid, self.scale, 'ok', markersize=10)
            pylab.plot(res.x[0], res.x[1], 'g', marker='D')

            pylab.title('Ex = {:.4} ; Ey={:.4}'.format(Ex, Ey))

        results = {'results':results, 'coords':coords, 'best':best,
                'Ex':Ex, 'Ey':Ey}

        return results

    def plot(self, results, coords, best, xmid_range, scale_range):
        pylab.figure(1)
        pylab.clf()
        pylab.contourf(results.transpose(), 20,
                extent=[xmid_range[0], xmid_range[-1],
                scale_range[0], scale_range[-1]], aspect='auto')
        pylab.plot(coords[0], coords[1], 'ko')
        pylab.colorbar()

        pylab.figure(2)
        pylab.clf()

        pylab.plot(self.signal.X, self.signal.Y, '--k', label='input')
        pylab.plot(self.signal.X, self.data, 'o', label='input')

        self.fitted = Logistic(coords[0], coords[1], N=self.N)
        self.fitted.X = self.X

        pylab.plot(self.fitted.X, self.fitted.Y, 'o-', label='input')
        pylab.legend(['signal','signal+noise','fitted'])

        pylab.grid()

    def objective_function(self, xmid, scale):
        temp = Logistic(xmid, scale)
        temp.X = self.signal.X
        M = self.get_snr(temp.Y, show=False).max()
        return M

    def optimise(self, guess=[1, 1]):
        from scipy.optimize import minimize
        def func(inputs):
            xmid = inputs[0]
            scale = inputs[1]
            M = self.objective_function(xmid, scale)
            return -M

        res = minimize(func, guess, method='nelder-mead',
            options={'xatol': 1e-8, 'disp': False})

        return res
