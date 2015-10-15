from gdsctools.tools import Logistic




def test_logistic():
    import pylab
    tl = Logistic(2, 1)
    tl.plot()
    tl.scale = 4
    tl.plot(hold=True)
    pylab.legend(['scale=1', 'scale=4'])

