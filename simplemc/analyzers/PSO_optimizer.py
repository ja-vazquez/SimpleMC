


import sys

try:
    import pyswarms as ps
    from pyswarms.single.global_best import GlobalBestPSO
    #from pyswarms.utils.plotters import (plot_cost_history,plot_contour,plot_surface)
except:
    sys.exit("*error: Install pyswarms library to use pso algorithms.")

from simplemc.plots.Plot_elipses import plot_elipses
import matplotlib.pyplot as plt
import scipy.linalg as la
import numpy as np
import scipy as sp

class PSO_optimizer():
    """
    explain some input params
    """

    def __init__(self,
                 like,
                 model,
                 outputname='pso_output',
                 technique=None,
                 nparticles=30,
                 iterations=80,
                 opt_c1=0.5,
                 opt_c2=0.5,
                 opt_w=0.9,
                 opt_k=3,
                 opt_p=1,
                 nproc=1,
                 use_oh_strategy=False,
                 early_stop=False,
                 ftol=-np.inf,
                 ftol_iter=1,
                 plot_fitness=False,
                 compute_errors=False,
                 show_contours=False,
                 plot_param1=None,
                 plot_param2=None,
                 save_like=False,
                 more_params=True):

        """Initialize the PSO Optimizer

        Attributes
        ----------
        n_particles : int
        number of particles in the swarm.
        """

        print('-- output --', outputname)

        self.like = like
        self.model = model
        self.outputname = outputname
        self.params = like.freeParameters()
        self.vpars = [p.value for p in self.params]
        self.npars = [p.name for p in self.params]
        self.sigma = sp.array([p.error for p in self.params])
        self.bounds = [p.bounds for p in self.params]
        self.pso_bounds = list(zip(*self.bounds))
        self.cov = None

        print("Minimizing...", self.npars, "starting", self.vpars)
        print("with bounds", self.bounds)
        print("#particles=", nparticles, "#iterations=", iterations, "\n")

        #PSO Algorithm hyperparams
        self.dimensions = len(self.params)
        self.nparticles = nparticles
        self.iterations = iterations
        self.technique = technique
        if use_oh_strategy:
            self.oh_strategy = {"w": 'exp_decay', 'c1': 'lin_variation'}

        self.nprocesses = nproc
        self.early_stop = early_stop
        self.ftol = ftol if early_stop else -np.inf
        self.ftol_iter = ftol_iter

        self.opt_c1 = opt_c1
        self.opt_c2 = opt_c2
        self.opt_w = opt_w
        self.opt_k = opt_k
        self.opt_p = opt_p

        self.plot_fitness = plot_fitness
        self.compute_errors = compute_errors
        self.show_contours = show_contours
        self.plot_param1 = plot_param1
        self.plot_param2 = plot_param2

        self.save_like = save_like


    def main(self, fout=None):
        """Optimize with the swarm

        Performs the optimization of the likelihood
        function :code:`f` for a number of iterations :code:`iter.`

        Parameters
        ----------
        objective_func : callable
            objective function to be evaluated
        iters : int
            number of iterations
        n_processes : int
            number of processes to use for parallel particle evaluation (default: None = no parallelization)
        kwargs : dict
            arguments for the objective function

        Returns
        -------
        dictionary

        """

        # limits in the format requested by the code:
        min_bound = self.pso_bounds[0]
        max_bound = self.pso_bounds[1]
        bounds = (min_bound, max_bound)

        if self.technique == 'Global':
            options = {'c1':self.opt_c1, 'c2':self.opt_c2, 'w':self.opt_w}
            optimizer = ps.single.GlobalBestPSO(n_particles=self.nparticles,
                                            dimensions=self.dimensions, ftol=self.ftol,
                                            options=options, bounds=bounds, center= self.vpars,
                                            oh_strategy=self.oh_strategy)

        elif self.technique == 'Local':
            options = {'c1':self.opt_c1, 'c2':self.opt_c2, 'w':self.opt_w,
                       'k':self.opt_k, 'p':self.opt_p}
            optimizer = ps.single.LocalBestPSO(n_particles=self.nparticles,
                                                dimensions=self.dimensions, ftol=self.ftol,
                                                options=options, bounds=bounds, center=self.vpars,
                                                oh_strategy=self.oh_strategy)

        else:
            sys.exit('no correct strategy')


        if self.early_stop:
            optimizer.ftol_iter = self.ftol_iter

        cost, pos = optimizer.optimize(self.negloglike, iters=self.iterations,
                                       n_processes=self.nprocesses)

        # Printing -only for testing purposes as it doubles the time
        self.formstr = '%g ' + '%g ' * len(self.vpars) + '\n'

        for i, j in enumerate(optimizer.pos_history[1::]):
            for k in j:
                if self.save_like:
                    s = np.concatenate((self.negloglike2(k), k), axis=None)
                else:
                    s=np.concatenate((i,k), axis=None)
                fout.write(self.formstr%tuple(s))


        if self.plot_fitness:
            self.plotting(optimizer)

        if self.compute_errors:
            try:
                #import numdifftools as nd
                import statsmodels.tools.numdiff as nd
            except:
                sys.exit("*'error': Install numdifftools to compute errors.")

            #hess = nd.Hessian(self.negloglike2, step=self.sigma*0.1)(pos)
            hess = nd.approx_hess3(pos, self.negloglike2, epsilon=self.sigma*0.1)
            #eigvl, eigvc = la.eig(hess)
            #print('Hessian', hess)
            #print('Eigen vals', eigvl)
            #print('Eigen vecs', eigvc)

            self.cov = la.inv(hess)
            print('Covariance matrix \n', self.cov)
            # set errors:
            for i, pars in enumerate(self.params):
                pars.setError(sp.sqrt(self.cov[i, i]))

            with open('{}.cov'.format(self.outputname), 'w') as f:
                np.savetxt(f, self.cov, fmt='%.4e', delimiter=',')


        if self.compute_errors and self.show_contours:
            param_names = [par.name for par in self.params]
            param_Ltx_names = [par.Ltxname for par in self.params]
            if (self.plot_param1 in param_names) and (self.plot_param2 in param_names):
                idx_param1 = param_names.index(self.plot_param1)
                idx_param2 = param_names.index(self.plot_param2)
                param_Ltx1 = param_Ltx_names[idx_param1]
                param_Ltx2 = param_Ltx_names[idx_param2]
            else:
                sys.exit('\n Not a base parameter, derived-errors still on construction')

            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111)
            plot_elipses(pos, self.cov, idx_param1, idx_param2, param_Ltx1, param_Ltx2, ax=ax)

            # Adding history
            import matplotlib as mpl

            niters = 5
            lh = len(optimizer.pos_history)

            c = np.arange(1, niters + 1)
            d = c*int(lh/niters)
            norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
            cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.cool) #jet, cool
            cmap.set_array([])

            #fig, ax = plt.subplots(dpi=100)
            for i, yi in enumerate(optimizer.pos_history[1::int(lh/niters)]):
                ax.plot(list(zip(*yi))[idx_param1], list(zip(*yi))[idx_param2], 'o', c=cmap.to_rgba(i))
            cbar = fig.colorbar(cmap, ticks=c, ax=ax)
            cbar.ax.set_yticklabels(d)
            cbar.set_label('Iterations', rotation=270)

            plt.savefig('pso_plot.pdf')
            plt.show()


        return {'population': self.nparticles, 'no_generations': self.iterations, 'param_fit':pos,
            'best_params': pos, 'cov': self.cov, 'maxlike': cost}


    def plotting(self, optimizer):
        plt.figure(figsize=(6, 6))

        plt.plot(np.arange(len(optimizer.cost_history[1::])), optimizer.cost_history[1::])
        plt.title("Cost History")
        plt.xlabel("Iterations", fontsize=20)
        plt.ylabel("Cost", fontsize=20)
        plt.grid()
        plt.legend(loc="upper right")

        # plt.savefig('GA_fitness.pdf')
        plt.show()


    def negloglike(self, xvals):
        nloglike = []
        for x in xvals:
            for i, pars in enumerate(self.params):
                pars.setValue(x[i])
            self.like.updateParams(self.params)
            loglike = self.like.loglike_wprior()
            nloglike.append(-loglike)
        return nloglike


    def negloglike2(self, x):
        #original likelihood, standard prior
        for i, pars in enumerate(self.params):
            pars.setValue(x[i])
        self.like.updateParams(self.params)
        loglike = self.like.loglike_wprior()
        return -loglike


if  __name__ == '__main__' :
    g= PSO_optimizer()
    g.main()
