import numpy as np

def glass(x, y):
    """Return Effect size through Glass D
    """
    x = np.array(x)
    y = np.array(y)

    # mean difference: 
    md = np.abs(x.mean() - y.mean()) 

    # here, we want same as in R that is unbiased variance
    # so we use ddof = 1
    g1 = md / x.std(ddof=1)
    g2 = md / y.std(ddof=1)

    return g1, g2
