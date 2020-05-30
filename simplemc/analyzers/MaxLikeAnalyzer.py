#
# This is an analyzer that takes a Likelihood function
# and then tries to maximize it and get the errors from the
# second derivative matrix. It kinda works, but not very well.
#

import scipy as sp
from scipy.optimize import minimize
import scipy.linalg as la

try:
    import numdifftools as nd
except:
    sys.exit('install numdifftools')


class MaxLikeAnalyzer:
    def __init__(self, like, withErrors=False):
        self.like   = like
        self.params = like.freeParameters()
        self.vpars  = [p.value for p in self.params]
        self.sigma  = sp.array([p.error for p in self.params])
        bounds = [p.bounds for p in self.params]
        print("Minimizing...", self.vpars, "with bounds", bounds)

        self.res = minimize(self.negloglike, self.vpars, bounds=bounds, method='L-BFGS-B')
        print(self.res, 'with noErrors =', withErrors)

        if (withErrors == 'Yes'):
            hess    = nd.Hessian(self.negloglike)(self.res.x)
            eigvl, eigvc = la.eig(hess)
            print (hess, eigvl,)
            self.cov = la.inv(hess)
            print('Covariance matrix \n', self.cov)
            # set errors:
            for i, pars in enumerate(self.params):
                pars.setError(sp.sqrt(self.cov[i, i]))
        # update with the final result
        self.result(self.negloglike(self.res.x))


    def negloglike(self, x):
        for i, pars in enumerate(self.params):
            pars.setValue(x[i])
        self.like.updateParams(self.params)
        loglike = self.like.loglike_wprior()

        if (sp.isnan(loglike)):
            return self.lastval+10
        else:
            self.lastval = -loglike
            return -loglike


    def result(self, loglike):
        print ("------")
        print("Done.")
        print("Optimal loglike : ", loglike)
