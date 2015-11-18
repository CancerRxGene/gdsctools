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
import numpy as np
import pylab
import os


class Savefig(object):
    """A simple class to save matploltib figures in the proper place

    .. note:: For developers only
    """
    def __init__(self, verbose=False):
        self.verbose = verbose
        self._directory = None
        #: directory where to save figures
        self.directory = '.'

    def _get_directory(self):
        return self._directory
    def _set_directory(self, directory):
        self._directory = directory
        try:
            if os.path.isdir(directory) is False:
                os.mkdir(self.directory)
                if self.verbose:
                    print("Created directory {}".format(directory))
        except Exception:
            if self.verbose:
                print("Could not create the directory")
    directory = property(_get_directory, _set_directory, doc="")

    def savefig(self, name, size_inches=None,**kargs):
        """Save a matplotlib figure

        :param str filename: where to save the figure.
        :param **kargs: accepts all parameters known by pylab.savefig
        """
        # if not provided, don't do anything.
        if name is None:
            return

        try:
            directory = self.settings.directory
        except:
            directory = self.directory

        filename = directory + os.sep + name
        fig = pylab.gcf()
        if size_inches is not None:
            fig.set_size_inches(size_inches)
        if self.verbose:
            print("saving file in %s" % filename)
        pylab.savefig(filename, **kargs)


class Logistic(object):
    """Simple logistic class to see the curve implied by xmid/scale parameters

    .. plot::
        :include-source:

        from gdsctools.tools import Logistic
        from pylab import legend
        tl = Logistic(2, 1)
        tl.plot()
        tl.scale = 4
        tl.plot(hold=True)
        legend(['scale=1', 'scale=4'])

    """
    def __init__(self, xmid, scale, Asym=1, N=9):
        r""".. rubric:: Constructor

        :param xmid:
        :param scale:
        :param Asym:

        .. math:: L(x) = \frac{Asym}{1 + exp ((xmid-X)/scale)}

        """
        self.xmid = xmid
        self.scale = scale
        self.Asym = Asym
        self._N = N

    def get_x(self):
        """Return a sensible linear space of X values """
        N = self._N
        dX = abs(self.scale)
        X = np.linspace(self.xmid - dX * self._N,
                self.xmid + dX * self._N, N)
        return  X

    def get_y(self, X=None):
        """Get the Y values given X and the 2 logistic function parameters"""
        if X is None:
            X = self.get_x()
        else:
            X = np.array(X)

        Y = self.Asym / (1.+ np.exp((self.xmid - X)/self.scale))
        return Y

    def plot(self, X=None, hold=False):
        """Plot the logistic function"""
        Y = self.get_y(X=X)
        import pylab
        if X is None:
            X = self.get_x()
        if hold is False:
            pylab.clf();
        pylab.plot(X, Y, 'o-')
        pylab.grid(True)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        pylab.axvline(self.xmid, color='r', linestyle='--')

    #def add_noise(self, sigma=0.1):
    #    noise = pylab.randn(self.N)
    #    self.data += sigma*noise

    def _get_N(self):
        return self._N
    def _set_N(self,N):
        self._N =N
    N = property(_get_N, _set_N)


class LogisticMatchedFiltering(object):

    def __init__(self, xmid, scale, xx=None,N=9, xmax=8):

        self.N = N
        self.xx = pylab.linspace(-xmax, xmax, N)
        self.data = Logistic(xmid=xmid, scale=scale, N=N).get_y(self.xx)
        self.xmid = xmid
        self.scale = scale

    def noisy_data(self, sigma=0.1):
        self.data = Logistic(xmid=self.xmid, scale=self.scale,
                N=self.N).get_y(self.xx)
        self.data += sigma * pylab.randn(self.N)

    def get_snr(self, template, show=True):
        #losc.ligo.org/tutorial_optimal

        # first take fft of data and template
        # we double vector size to next power of 2
        N = 32
        data = self.data.copy()
        data = np.concatenate([data, [0] * (32-self.N+1)], axis=0)
        data_fft = np.fft.fft(data)

        #For this to work, we need the template and the data to be the
        # same length. So, we'll zero-pad the template before we take the FFT:
        template = np.concatenate([template, [0] * (32-self.N+1)], axis=0)
        template_fft = np.fft.fft(template)

        # we need an estimate of the noise power in each FFT bin.
        # we can assume white noise.

        # -- Calculate the PSD of the data
        #power_data, freq_psd = plt.psd(data[12*fs:], Fs=fs, NFFT=fs, visible=False)

        # -- Interpolate to get the PSD values at the needed frequencies
        fs = 1
        datafreq = np.fft.fftfreq(data.size)*fs
        #power_vec = np.interp(datafreq, freq_psd, power_data)
        power_vec  = [1] * 33   ## 32 +f0
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

    def scan(self, xmid_range, scale_range):

        N = len(xmid_range)
        results = np.zeros((N,N))
        coords = None
        best = 0
        for i,xmid in enumerate(xmid_range):
            for j,scale in enumerate(scale_range):
                temp = Logistic(xmid, scale).get_y(X=self.xx)
                M = self.get_snr(temp, show=False).max()
                results[i,j] = M
                if M>best:
                    best=M
                    coords = (xmid,scale)

        self.plot(results, coords, best, xmid_range, scale_range)

        res = self.optimise(coords)
        pylab.figure(1)
        pylab.plot(res.x[0], res.x[1], 'go')
        
        pylab.title('Ex = {:.4} ; Ey={:.4}'.format(
            abs(self.xmid-res.x[0])/self.xmid*100,
            abs(self.scale-res.x[1])/self.scale*100))
        return results, coords, best

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
        oridata = Logistic(xmid=self.xmid, scale=self.scale,
                N=self.N).get_y(self.xx)
        pylab.plot(self.xx, oridata,'--k', label='input')
        pylab.plot(self.xx, self.data,'o', label='input')

        fitted = Logistic(coords[0], coords[1], N=self.N)
        pylab.plot(self.xx, fitted.get_y(self.xx),'o-', label='input')

        pylab.grid()

    def objective_function(self, xmid, scale):
        temp = Logistic(xmid, scale).get_y(X=self.xx)
        M = self.get_snr(temp, show=False).max()
        return M

    def optimise(self, guess=[1,1]):
        from scipy.optimize import minimize
        def func(inputs):
            xmid = inputs[0]
            scale = inputs[1]
            M = self.objective_function(xmid, scale)
            return -M


        res = minimize(func, guess, method='nelder-mead',
            options={'xtol': 1e-8, 'disp': True})

        
        return res













