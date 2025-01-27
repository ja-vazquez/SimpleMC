
#TODO: get back self.postprocess() on Driver

import pyswarms as ps
from pyswarms.single.global_best import GlobalBestPSO
#from pyswarms.utils.plotters import (plot_cost_history,plot_contour,plot_surface)

import scipy as sp

class PSO_optimizer():
    #explain some input params

    def __init__(self, like, model, more_params=True):

        self.like = like
        self.model = model

        self.params = like.freeParameters()
        self.vpars = [p.value for p in self.params]
        self.sigma = sp.array([p.error for p in self.params])
        self.bounds = [p.bounds for p in self.params]
        self.pso_bounds = list(zip(*self.bounds))
        print("Minimizing...", self.vpars, "with pso bounds", self.pso_bounds)
        self.cov = None

        self.dimensions = len(self.params)
        self.nparticles = 30
        self.iterations = 80

    def main(self):
        # limits in the format requested by the code:

        min_bound = self.pso_bounds[0]
        max_bound = self.pso_bounds[1]
        bounds = (min_bound, max_bound)

        options = {'c1': 0.5, 'c2': 0.5, 'w': 0.9}

        optimizer = ps.single.GlobalBestPSO(n_particles=self.nparticles,
                                            dimensions=self.dimensions,
                                            options=options, bounds=bounds)

        cost, pos = optimizer.optimize(self.negloglike, iters=self.iterations)

        return {'population': self.nparticles, 'no_generations': self.iterations, 'param_fit':pos,
            'best_params': pos, 'cov': self.cov, 'maxlike': cost}


    def negloglike(self, xvals):
        nloglike = []
        for x in xvals:
            for i, j in enumerate(x):
                self.params[i].setValue(j)
            self.like.updateParams(self.params)
            loglike = self.like.loglike_wprior()
            nloglike.append(-loglike)
        return nloglike



if  __name__ == '__main__' :
    g= PSO_optimizer()
    g.main()
