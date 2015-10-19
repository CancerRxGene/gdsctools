import numpy as np

def cohens(x, y):
    """Effect size through Cohen's D
    
    :param x:
    :param y:
    :return: the effect size of the interaction between the values 
        in the set given by the union of x and y and the factor
        inducing the dychotomy.

    """
    x = np.array(x)
    y = np.array(y)

    Nx = len(x) - 1.  # note the dot to cast to float
    Ny = len(y) - 1.
    # mean difference: 
    md = np.abs(x.mean() - y.mean()) 

    # here, we want same as in R that is unbiased variance
    # so we use ddof = 1
    xv = x.var(ddof=1)
    yv = y.var(ddof=1)
    csd = Nx * xv + Ny * yv
    csd /= Nx + Ny  # make sure this is float
    csd = np.sqrt(csd)

    return md / csd
