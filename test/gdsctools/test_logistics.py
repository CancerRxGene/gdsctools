from gdsctools.logistics import Logistic, LogisticMatchedFiltering
import pylab


def test_logistic():
    import pylab
    tl = Logistic(2, 1)
    tl.plot()
    tl.scale = 4
    tl.plot(hold=True)
    pylab.legend(['scale=1', 'scale=4'])

    # settter
    tl.xmin= -4
    tl.xmax = 4
    tl.N = 40
    tl.X = [-1,0,1]
    assert tl.N == 3

def test_logistic_mf():
    mf = LogisticMatchedFiltering(1,2)
    mf.scan(pylab.linspace(-3,3, 10), pylab.linspace(0,5,10))


