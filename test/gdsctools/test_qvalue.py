from gdsctools import qvalue
import numpy as np
from easydev import assert_list_almost_equal

def test_qvalue():


    pvalues = np.array([0.8, 0.5, 0.5] + list(10**-np.linspace(1,10,9)))

    # This gives same answer as in R qvalue  library

    # pvalues = c(c(0.8,0.5,0.5), 10**(-seq(1,10,2)))
    # qvalue(pvalues)

    # Note that here we enforce the pi0 value to be the same as in R
    # If not, this particular fails because pi0 is negative...
    # qvalue not very robust.
    qv = qvalue.QValue(list(pvalues), pi0=0.1109898)
    qv = qvalue.QValue(pvalues, pi0=0.1109898)

    assert_list_almost_equal(qv.qvalue(),
            np.array([  8.87918400e-02,   6.05398909e-02,     6.05398909e-02,
                1.47986400e-02,   1.24845912e-03,   1.06995688e-04,
                9.36080212e-06,   8.42353356e-07,  7.89594880e-08,
                7.89483504e-09,   8.88043662e-10,  1.33187760e-10]))


    try:
        qv = qvalue.QValue(pvalues)
        assert False
    except:
        assert True
    try:
        qv = qvalue.QValue(pvalues, pi0=0.1109898, lambdas=[0,.5,0.7])
        assert False
    except:
        assert True
    try:
        qv = qvalue.QValue(pvalues, pi0=0.1109898, lambdas=[0,.5,0.7,.9,500])
        assert False
    except:
        assert True
    try:
        qv = qvalue.QValue(pvalues, pi0=0.1109898, lambdas=[-10,.5,0.7,.9])
        assert False
    except:
        assert True
