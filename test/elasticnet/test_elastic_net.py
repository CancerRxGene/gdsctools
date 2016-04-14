from sklearn import linear_model
import numpy as np
from easydev import assert_list_almost_equal as alae


X = np.array([[0,0],[1,1],[2,2],[3,4]])
Y = [0, 1, 2, 4]
myalphas = [0.08633579, 0.07866596, 0.07167749, 0.06530986, 0.05950791, 0.05422139, 0.04940451, 0.04501555, 0.0410165, 0.0373727, 0.03405262, 0.03102747, 0.02827108, 0.02575955, 0.02347115, 0.02138603, 0.01948616, 0.01775506, 0.01617775, 0.01474056, 0.01343105, 0.01223788, 0.0111507, 0.0101601, 0.009257504, 0.008435093, 0.007685742, 0.007002962, 0.006380838, 0.005813982, 0.005297484, 0.00482687, 0.004398064, 0.004007352, 0.00365135, 0.003326974, 0.003031415, 0.002762113, 0.002516734, 0.002293154, 0.002089437, 0.001903817, 0.001734687, 0.001580582, 0.001440168, 0.001312227, 0.001195653, 0.001089434, 0.0009926518, 0.0009044673, 0.000824117, 0.0007509047, 0.0006841964, 0.0006234143, 0.0005680319, 0.0005175695, 0.00047159, 0.0004296953, 0.0003915223, 0.0003567406, 0.0003250487, 0.0002961723, 0.0002698612, 0.0002458874, 0.0002240435, 0.0002041401, 0.0001860048, 0.0001694807, 0.0001544245, 0.0001407058, 0.0001282059, 0.0001168165, 0.0001064388, 9.698307e-05, 8.836736e-05, 8.051705e-05, 7.336414e-05, 6.684667e-05, 6.090819e-05, 5.549728e-05]


"""Elastic Net::


    1 / (2 * n_samples) * ||y - Xw||^2_2 
        + alpha * l1_ratio * ||w||_1
        + 0.5 * alpha * (1 - l1_ratio) * ||w||^2_2

and glmnet::      

    1 / (2 * nobs) * RSS
    lambda * [ (1-alpha)/2||beta||_2^2 + alpha||beta||_1. ]

donc python  ----- R 

       l1_ratio      alpha
       alpha         lambda

RIDGE    l1_ratio=0   alpha=0
LASSO    l1_ratio=1   alpha=1

"""


def test_ridge():
    """

    Ridge from sklearn is the same as the R package penalised
    but differs from glmnet, whose regularisation on beta_0 seems
    to be different as compared to penalized R package and sklearn package.


    glmnet and ElasticNet gives different results and different from 
    Ridge and penalized
    """
    R_code_penalized = """
    library(penalized)

    x = matrix(c(0,1,2,3,0,1,2,4), ncol=2)

    coefficients(
        penalized(c(0,1,2,4), x, lambda1=0,lambda2=1, standardize=F))
    0.04615385  0.40000000  0.63076923

    """ 
    ridge = linear_model.Ridge(alpha=1)

    alae(ridge.fit(X, Y).coef_ , np.array([0.4, 0.63076923]), places=6)
    alae([ridge.fit(X, Y).intercept_], [0.04615384615384599], places=10)

    # In ElasticNet, the doc says that l1_ratio < 0.01 is not reliable except
    # if one set its own alphas. It appears indeed that if we set l1_ratio to 0, 
    # we do not get the same results as in Ridge or penalized function.

    # Note however, that using alpha =0.01 gives better results and of course 
    # alpha=0 (pure OLS) gives the same results.

    # en = linear_model.ElasticNet(l1_ratio=0., alpha=1); en.fit(X,Y);
    # print(en.intercept_, en.coef_) 
    # (0.33103448275862091, array([ 0.35862069,  0.50344828]))
    # THis depends on the value of alpha. 
    # If alpha is small, we get closer results !

    # glmnet:
    # fit = glmnet(x, c(0,1,2,4), alpha=0.0, lambda=1, standardize=F);
    # fit$a0, fit$beta
    # 0.2193916, [0.3822875, 0.5469584]
    # one can play with lambda e.g. 0.370576 gives:
    # 0.0461538, [0.4007532, 0.6301236)


def test_lasso():
    """

    Lasso, ElasticNet (using l1_ratio=1) and glmnet (alpha=1) gives the same 
    results.

    scikit-learn  Lasso documentation::

        (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1


    In glmnet::

      lambda * [ (1-alpha)/2||beta||_2^2 + alpha||beta||_1. ]
    
    so if we set alpha = 1, the first term (ridge) disappears

    Using the penalised function, we get different results...that is
    intercept is 0.2 and coeff are [0 (0.8857143)
    """

    R_code_glmnet = """
 
    fit = glmnet(x, c(0,1,2,4), alpha=1, lambda=1, standardize=F)
    c(fit$a0, fit$beta)
    $s0
    [1] 0.8
    V1 = ( ., V2 0.5428571)
    """

    lasso = linear_model.Lasso(alpha=1)

    alae(lasso.fit(X, Y).coef_ , np.array([0., 0.54285714]), places=6)
    alae([lasso.fit(X, Y).intercept_], [0.8], places=10)

    # The code using ElasticNet class is 
    en = linear_model.ElasticNet(l1_ratio=1, alpha=1)
    en.fit(X, Y)

    alae([0.8], [en.intercept_], places=6)
    alae(en.coef_, np.array([ 0.        ,  0.54285714]), places=6)



def _test_elastic_net():
    en = linear_model.ElasticNet(l1_ratio=0.5, alpha=1)
    en.fit(X, Y); print(en.intercept_, en.coef_)
    (0.59088726449594753, array([ 0.13641303,  0.54542468]))

    """

    fit = glmnet(x, c(0,1,2,4), alpha=0.5, lambda=1, standardize=F)
    c(fit$a0,fit$beta)
    0.5534667,  0.07282901, 0.62130846

    coefficients(penalized(c(0,1,2,4), x, lambda1=0.5,  lambda2=0.5, standardize=F), "all")
    # nonzero coefficients: 3
    (Intercept)                         
    0.1304348   0.2173913   0.7391304 


    """

    """
    Effect of alpha now st to 0.1
    SKLEarn:
    0.034782344700849954 ,     array([ 0.18960879,  0.81760256])

    glmnet

    coefficients(glmnet(x, y, alpha=0.5, lambda=0.1, standardize=F))
                   s0
    (Intercept) 0.0388365
    V1          0.1237032
    V2          0.8717764

    and set to 0.01 :
    sklearn: -0.002135 0.0498 0.9584
    glmnet: -0.00041117, 0.035, 0.9700


    """

    pass


#def test_enet_path():
#    alphas, coefs, _ = enet_path(X, Y, l1_ratio=0.5, alphas=myalphas)











