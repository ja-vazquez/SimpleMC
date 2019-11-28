#
# This is an analyzer that takes a Likelihood function
# and then tries to maximize it and get the errors from the
# second derivative matrix. It kinda works, but not very well.
#

import scipy as sp
from scipy.optimize import minimize
import scipy.linalg as la
#try:
import numdifftools as nd
#except:
#    pass


class MaxLikeAnalyzer:
    def __init__(self, like, noErrors=False):
        self.like   = like
        self.params = like.freeParameters()
        self.vpars  = [p.value for p in self.params]
        self.sigma  = sp.array([p.error for p in self.params])
        print("Minimizing...", self.vpars)
        bounds = [p.bounds for p in self.params]
        print(bounds)
        res = minimize(self.negloglike, self.vpars, bounds=bounds, method='L-BFGS-B')
        print(res)

        if (not noErrors):
            hess    = nd.Hessian(self.negloglike, step=self.sigma)(res.x)
            print ('Working?', '**'*20)
            eigvl, eigvc = la.eig(hess)
            print (hess, eigvl,)
            self.cov = la.inv(hess)
            print ('Thought so', '**'*20)
            print(self.cov)
            # set errors:
            for i, pars in enumerate(self.params):
                pars.setError(sp.sqrt(self.cov[i, i]))
        # update with the final result
        self.negloglike(res.x)
        print("Done.")



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
