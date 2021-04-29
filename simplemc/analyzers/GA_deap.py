

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from simplemc.plots.Plot_elipses import plot_elipses

try:
    # Importamos libreria de algoritmos evolutivos
    from deap import base, creator, tools, algorithms
except:
    import warnings
    warnings.warn("Pleas install DEAP library if you want to use ga_deap genetic algorithms.")

# We import an independent module to implement elitism in the GA.
from simplemc.analyzers import elitism

import scipy.linalg as la

import random
import sys

try:
    import numdifftools as nd
except:
    sys.exit('install numdifftools')

class GA_deap:
    def __init__(self, like, model, plot_fitness=False, compute_errors=False, \
                 show_contours=False, plot_param1=None, plot_param2=None):

        self.like = like
        self.model = model
        self.params = like.freeParameters()
        self.vpars = [p.value for p in self.params]
        self.sigma = sp.array([p.error for p in self.params])
        self.bounds = [p.bounds for p in self.params]
        print("Minimizing...", self.vpars, "with bounds", self.bounds)

        self.plot_fitness = plot_fitness
        self.compute_errors = compute_errors
        self.show_contours = show_contours
        self.plot_param1 = plot_param1
        self.plot_param2 = plot_param2

        # Genetic Algorithm constants:
        self.POPULATION_SIZE = 20    # 10-20
        self.P_CROSSOVER = 0.7       # probability for crossover
        self.P_MUTATION = 0.3        # (try also 0.5) probability for mutating an individual
        self.MAX_GENERATIONS = 20    # 100- 300
        self.HALL_OF_FAME_SIZE = 1
        self.CROWDING_FACTOR = 20.0  # crowding factor for crossover and mutation

        self.RANDOM_SEED = 42        # set the random seed

        self.DIMENSIONS = len(self.params)        # number of dimensions
        self.BOUND_LOW, self.BOUND_UP = 0.0, 1.0  # boundaries for all dimensions




    def main(self):
        toolbox= self.GA()

        # create initial population (generation 0):
        population = toolbox.populationCreator(n=self.POPULATION_SIZE)

        # prepare the statistics object:
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("min", np.min)
        stats.register("avg", np.mean)

        # define the hall-of-fame object:
        hof = tools.HallOfFame(self.HALL_OF_FAME_SIZE)

        # perform the Genetic Algorithm flow with elitism:
        population, logbook = elitism.eaSimpleWithElitism(population, toolbox, cxpb=self.P_CROSSOVER,\
                                                          mutpb=self.P_MUTATION, ngen=self.MAX_GENERATIONS,\
                                                          stats=stats, halloffame=hof, verbose=True)

        # print info for best solution found:
        best = hof.items[0]
        print("-- Best Fitness = ", best.fitness.values[0])
        print("- Best solutions are:")
        best_params = [self.change_prior(i, x) for i, x in enumerate(best)]
        for i, x in enumerate(best_params):
            print("-- Best %s = "%self.params[i].name , x)


        #for i in range(self.HALL_OF_FAME_SIZE):
        #    print(i, ": ", hof.items[i].fitness.values[0], " -> ", self.old_prior(i, hof.items[i]) )

        if self.plot_fitness:
            self.plotting(population, logbook, hof)


        if self.compute_errors:
            hess = nd.Hessian(self.negloglike2)(best_params)
            eigvl, eigvc = la.eig(hess)
            print ('Hessian', hess, eigvl,)
            self.cov = la.inv(hess)
            print('Covariance matrix \n', self.cov)
            # set errors:
            #for i, pars in enumerate(self.params):
            #    pars.setError(sp.sqrt(self.cov[i, i]))
        # update with the final result
        #self.result(self.negloglike(self.res.x))

        if self.show_contours and self.compute_errors:
            param_names = [par.name for par in self.params]
            if (self.plot_param1 in param_names) and (self.plot_param2 in param_names):
                idx_param1 = param_names.index(self.plot_param1)
                idx_param2 = param_names.index(self.plot_param2)
            else:
                sys.exit('\n Not a base parameter, derived-errors still on construction')

            fig = plt.figure(figsize=(6,6))
            ax = fig.add_subplot(111)
            plot_elipses(best_params, self.cov, idx_param1, idx_param2, ax=ax)
            plt.show()
        return population, logbook, hof



    # helper function for creating random real numbers uniformly distributed within a given range [low, up]
    # it assumes that the range is the same for every dimension
    def randomFloat(self, low, up):
        return [random.uniform(l, u) for l, u in zip([low]*self.DIMENSIONS, [up]*self.DIMENSIONS)]




    def GA(self):
        random.seed(self.RANDOM_SEED)
        toolbox = base.Toolbox()

        # define a single objective, minimizing fitness strategy:
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))

        # create the Individual class based on list:
        creator.create("Individual", list, fitness=creator.FitnessMin)

        # create an operator that randomly returns a float in the desired range and dimension:
        toolbox.register("attrFloat", self.randomFloat, self.BOUND_LOW, self.BOUND_UP)

        # create the individual operator to fill up an Individual instance:
        toolbox.register("individualCreator", tools.initIterate, creator.Individual, toolbox.attrFloat)

        # create the population operator to generate a list of individuals:
        toolbox.register("populationCreator", tools.initRepeat, list, toolbox.individualCreator)

        ## ----- here the connection
        toolbox.register("evaluate", self.negloglike)
        ## -----

        # genetic operators:
        toolbox.register("select", tools.selTournament, tournsize=2)
        toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=self.BOUND_LOW, \
                         up=self.BOUND_UP, eta=self.CROWDING_FACTOR)
        toolbox.register("mutate", tools.mutPolynomialBounded, low=self.BOUND_LOW, \
                         up=self.BOUND_UP, eta=self.CROWDING_FACTOR, indpb=1.0/self.DIMENSIONS)

        return toolbox



    def plotting(self, pop, log, hof):

        # extract statistics
        gen, avg, min_, max_ = log.select("gen", "avg", "min", "max")

        plt.figure(figsize=(10, 7))

        plt.plot(gen, min_, label="minimum")

        plt.title("Fitness Evolution")
        plt.xlabel("Generation")
        plt.ylabel("Fitness")
        plt.legend(loc="upper right")
        #plt.savefig('GA_wwCDM_150.pdf')
        plt.show()



    def negloglike(self, x):
        for i, pars in enumerate(self.params):
            new_par = self.change_prior(i, x[i])
            pars.setValue(new_par)
        self.like.updateParams(self.params)
        loglike = self.like.loglike_wprior()

        if sp.isnan(loglike):
            print('-1-'*10,loglike,'--'*10)
            return self.lastval+10
        else:
            self.lastval = -loglike
        return -loglike,


    def negloglike2(self, x):
        for i, pars in enumerate(self.params):
            pars.setValue(x[i])
        self.like.updateParams(self.params)
        loglike = self.like.loglike_wprior()

        if sp.isnan(loglike):
            return self.lastval+10
        else:
            self.lastval = -loglike
        return -loglike



    def change_prior(self, i, x):
        new_par = self.bounds[i][0] + (self.bounds[i][1] - self.bounds[i][0])*x
        return new_par

    def old_prior(self, i, x):
        old_par = (x - self.bounds[i][0])/(self.bounds[i][1] - self.bounds[i][0])
        return old_par


    def result(self, loglike):
        print ("------")
        print("Done.")
        print("Optimal loglike : ", loglike)


if  __name__ == '__main__' :
    g= GA_deap()
    g.GA()