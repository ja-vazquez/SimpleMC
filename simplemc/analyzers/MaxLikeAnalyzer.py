#
# This is an analyzer that takes a Likelihood function
# and then tries to maximize it and get the errors from the
# second derivative matrix. It kinda works, but not very well.
#
from matplotlib.patches import Ellipse
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import scipy.linalg as la
import scipy as sp
import numpy as np
import sys

try:
    import numdifftools as nd
except:
    sys.exit('install numdifftools')


class MaxLikeAnalyzer:
    def __init__(self, like, model, withErrors=False, param1=False, param2=False):
        self.like   = like
        self.model  = model
        self.params = like.freeParameters()
        self.vpars  = [p.value for p in self.params]
        self.sigma  = sp.array([p.error for p in self.params])
        bounds = [p.bounds for p in self.params]
        print("Minimizing...", self.vpars, "with bounds", bounds)

        self.res = minimize(self.negloglike, self.vpars, bounds=bounds, method='L-BFGS-B')
        print(self.res, 'with noErrors =', withErrors)

        if (withErrors):
            hess    = nd.Hessian(self.negloglike)(self.res.x)
            eigvl, eigvc = la.eig(hess)
            print ('Hessian', hess, eigvl,)
            self.cov = la.inv(hess)
            print('Covariance matrix \n', self.cov)
            # set errors:
            for i, pars in enumerate(self.params):
                pars.setError(sp.sqrt(self.cov[i, i]))
        # update with the final result
        self.result(self.negloglike(self.res.x))


        if param1 and param2:
            for i, par in enumerate(self.like.freeParameters()):
                if param1 == par.name: par1 = i
                elif param2 == par.name: par2 = i


            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111) #, aspect='equal')

            val,vec = la.eig(self.cov[[par1,par2], :][:, [par1,par2]])
            vec=vec.T
            mn = self.res.x

            vec[0]*=sp.sqrt(11.83*sp.real(val[0]))
            vec[1]*=sp.sqrt(11.83*sp.real(val[1]))

            plt.plot(mn[par1],mn[par2],'bo', label=self.model)
            plt.plot([mn[par1]-vec[0][0],mn[par1]+vec[0][0]],
                [mn[par2]-vec[0][1],mn[par2]+vec[0][1]],'r-')
            plt.plot([mn[par1]-vec[1][0],mn[par1]+vec[1][0]],
                [mn[par2]-vec[1][1],mn[par2]+vec[1][1]],'r-')

            def eigsorted(cov):
                vals, vecs = sp.linalg.eigh(cov)
                order = vals.argsort()[::-1]
                return vals[order], vecs[:,order]

            vals, vecs = eigsorted(self.cov[[par1,par2], :][:, [par1,par2]])
            theta = sp.degrees(np.arctan2(*vecs[:,0][::-1]))

            sigmas = [2.3, 5.99, 11.83]
            for sigs in sigmas:
                w, h = 2  * sp.sqrt(vals) * sp.sqrt(sigs)
                ell = Ellipse(xy=(mn[par1], mn[par2]),  width = w, height = h, angle=theta, color='green',  lw=4)
                ell.set_facecolor('none')
                ax.add_artist(ell)

            plt.legend(loc='best')
            plt.title('+'.join(self.like.compositeNames()), fontsize=10)
            plt.show()



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
