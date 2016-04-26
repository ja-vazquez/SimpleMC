#
# This is an analyzer that takes a Likelihood function
# and then tries to maximize it and get the errors from the
# second derivative matrix. It kinda works, but not very well.
#

from scipy import *
from scipy.optimize import minimize
import scipy.linalg as la
try:
    import numdifftools as nd
except:
    pass

class MaxLikeAnalyzer:
    def __init__ (self,like,noErrors=False):
        self.like=like
        self.params=like.freeParameters()
        self.vpars=[p.value for p in self.params]
        print "Minimizing...", self.vpars
        bounds=[p.bounds for p in self.params]
        print bounds
        res=minimize(self.negloglike,self.vpars,bounds=bounds,method='L-BFGS-B')#, bounds=bounds)
        print res
        if (not noErrors):
            hessfun=nd.Hessian(self.negloglike)
            hess=hessfun(res.x)
            self.cov=la.inv(hess)
            print self.cov
            ## set errors:
            for i in range(len(self.params)):
                self.params[i].setError(sqrt(self.cov[i,i]))
        ## update with the final result
        self.negloglike(res.x) 
        print "Done."

    def negloglike(self,x):
        for i in range(len(self.params)):
            self.params[i].setValue(x[i])
        self.like.updateParams(self.params)
        loglike=self.like.loglike_wprior()
        #print x,loglike

        if (isnan(loglike)):
            return self.lastval+10
        else:
            self.lastval=-loglike
            return -loglike
        
    
