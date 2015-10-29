import numpy as np


__all__ = ['cohens']


def cohens(x, y):
    r"""Effect size metric through Cohen's *d* metric

    :param x: first vector
    :param y: second vector
    :return: absolute effect size value

    The Cohen's effect size *d* is defined as the difference
    between two means divided by a standard deviation of the data.

    .. math::

        d = \frac{\bar{x}_1 - \bar{x}_2}{s}

    For two independent samples, the *pooled standard deviation* is used
    instead, which is defined as:

    .. math::

        s = \sqrt{  \frac{(n_1-1)s_1^2 + (n_2-1)s_2^2}{n_1+n_2-2} }


    A Cohen's *d* is frequently used in estimating sample sizes for
    statistical testing: a lower *d* value indicates the necessity of
    larger sample sizes, and vice versa.

    .. note:: we return the absolute value

    :references: https://en.wikipedia.org/wiki/Effect_size

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
