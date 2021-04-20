

from simplemc.plots.Plot_elipses import plot_elipses
from simplemc.cosmo.Derivedparam import AllDerived
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import scipy.linalg as la
import scipy as sp
import sys


try:
    import numdifftools as nd
except:
    sys.exit('install numdifftools')



class MaxLikeAnalyzer:
    def __init__(self, like, model, compute_errors=False, compute_derived= False,
                 show_contours=False, plot_param1=None, plot_param2=None):
        """
        This is an analyzer that takes a Likelihood function
        and then tries to maximize it and get the errors from the
        second derivative matrix. It kinda works, but not very well.
        Parameters
        ----------
        like
        model
        compute_errors
        plot_param1
        plot_param1

        Returns
        -------

        """
        self.like = like
        self.model = model
        self.params = like.freeParameters()
        self.vpars = [p.value for p in self.params]
        self.sigma = sp.array([p.error for p in self.params])
        bounds = [p.bounds for p in self.params]
        print("Minimizing...", self.vpars, "with bounds", bounds)

        self.res = minimize(self.negloglike, self.vpars, bounds=bounds, method='L-BFGS-B')
        print(self.res, 'with Errors =', compute_errors)


        if compute_derived:
            for par, val in zip(self.params, self.res.x):
                par.setValue(val)
            self.like.updateParams(self.params)

            print('--'*5, 'Derived parameters', '--'*5)
            for der_par in AllDerived().listDerived(self.like):
                print(der_par.name, ': ', der_par.value)
            print('--'*20)
            # for errors, consider the df = J cov_x J^t


        if (compute_errors):
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




        if show_contours and compute_errors:
            param_names = [par.name for par in self.params]
            if (plot_param1 in param_names) and (plot_param2 in param_names):
                idx_param1 = param_names.index(plot_param1)
                idx_param2 = param_names.index(plot_param2)
            else:
                sys.exit('\n Not a base parameter, derived-errors still on construction')

            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111)
            plot_elipses(self.res.x, self.cov, idx_param1, idx_param2, ax=ax)
            plt.show()



    def negloglike(self, x):
        for i, pars in enumerate(self.params):
            pars.setValue(x[i])
        self.like.updateParams(self.params)
        loglike = self.like.loglike_wprior()

        if sp.isnan(loglike):
            return self.lastval+10
        else:
            self.lastval = -loglike
            return -loglike


    def result(self, loglike):
        print ("------")
        print("Done.")
        print("Optimal loglike : ", loglike)

