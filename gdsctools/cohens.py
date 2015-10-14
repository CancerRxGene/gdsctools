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

    Nx = len(x) - 1
    Ny = len(y) - 1
    # mean difference: 
    md = np.abs(x.mean() - y.mean()) 

    # here, we want same as in R that is unbiased variance
    # so we use ddof = 1
    csd = Nx * x.var(ddof=1) + Ny * y.var(ddof=1)
    csd /= float(Nx + Ny)  # make sure this is float
    csd = np.sqrt(csd)

    return md / csd
