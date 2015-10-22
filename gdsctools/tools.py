import numpy as np
import pylab
import os

class Savefig(object):
    def __init__(self):
        # This should be a parent class. 
        # The child should have a directory attribute 
        # or settings.directory. If not, default to local directory.
        self.directory = '.'

    def savefig(self, name, **kargs):
        try:
            directory = self.settings.directory
        except:
            directory = self.directory

        filename = directory + os.sep + name
        print('Saving %s' % filename)
        pylab.savefig(filename, **kargs)


class Logistic(object):
    """Simple logistic class to see the curve implied by xmid/scale parameters

    .. plot::
        :include-source:

        from gdsctools.tools import Logistic
        tl = Logistic(2, 1)
        tl.plot()
        tl.scale = 4
        tl.plot(hold=True)
        legend(['scale=1', 'scale=4'])



    """
    def __init__(self, xmid, scale, Asym=1):
        self.xmid = xmid
        self.scale = scale
        self.Asym = Asym
        self._N = 10

    def getX(self, N=100):

        dX = abs(self.scale) 
        X = np.linspace(self.xmid - dX * self._N, 
                self.xmid + dX * self._N, N)
        return  X

    def getY(self, X=None, N=100):
        if X is None:
            X = self.getX()
        else:
            X = np.array(X)

        Y = self.Asym / (1.+ np.exp((self.xmid - X)/self.scale))
        return Y

    def plot(self, X=None, N=100, hold=False):
        Y = self.getY(X=X, N=N)
        import pylab
        if X is None:
            X = self.getX()
        if hold is False:
            pylab.clf();
        pylab.plot(X, Y, 'o-')
        pylab.grid(True)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        pylab.axvline(self.xmid, color='r', linestyle='--')
