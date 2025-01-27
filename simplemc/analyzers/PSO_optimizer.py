
#TODO: get back self.postprocess() on Driver

import pyswarms as ps
from pyswarms.single.global_best import GlobalBestPSO
#from pyswarms.utils.plotters import (plot_cost_history,plot_contour,plot_surface)

import scipy as sp

class PSO_optimizer():
    #explain some input params

    def __init__(self, like, model, outputname='pso_output',
                 nparticles=30, iterations= 80,
                 opt_c1=0.5, opt_c2=0.5, opt_w= 0.9,
                 plot_fitness=False, more_params=True):

        self.like = like
        self.model = model
        self.outputname = outputname

        self.params = like.freeParameters()
        self.vpars = [p.value for p in self.params]
        self.sigma = sp.array([p.error for p in self.params])
        self.bounds = [p.bounds for p in self.params]
        self.pso_bounds = list(zip(*self.bounds))
        print("Minimizing...", self.vpars, "with pso bounds", self.pso_bounds)
        self.cov = None

        self.plot_fitness = plot_fitness

        #PSO Algorithm hyperparams
        self.dimensions = len(self.params)
        self.nparticles = nparticles
        self.iterations = iterations

        self.opt_c1 = opt_c1
        self.opt_c2 = opt_c2
        self.opt_w = opt_w

    def main(self):
        # limits in the format requested by the code:

        min_bound = self.pso_bounds[0]
        max_bound = self.pso_bounds[1]
        bounds = (min_bound, max_bound)

        options = {'c1': self.opt_c1, 'c2': self.opt_c2, 'w': self.opt_w}

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
