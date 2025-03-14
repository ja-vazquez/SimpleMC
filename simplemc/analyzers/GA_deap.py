
#TODO make processes a variable
#TODO check/change select, mate, mutate


from simplemc.plots.Plot_elipses import plot_elipses
import matplotlib.pyplot as plt
import scipy.linalg as la
import multiprocessing
import scipy as sp
import numpy as np
import random
import sys

try:
    # Importamos libreria de algoritmos evolutivos
    from deap import base, creator, tools, algorithms
except:
    sys.exit("*error: Install DEAP library to use genetic algorithms (ga_deap).")

# We import an independent module for elitism in the GA.
from simplemc.analyzers import elitism


class GA_deap:
    """
        Genetic algorithms from deap library

        :param like: likelihood function
        :param model: theoretical model model
        :param plot_fitness: True if you want to generate a plot of the fitness.
        :param compute_errors: True if you want to compute errors.
        :param show_contours: True if you want to show the contours in a plot.
        :param plot_param1: a parameter to plot in x-axis.
        :param plot_param2: a parameter to plot in y-axis.

    """
    def __init__(self, like, model, outputname='deap_output',
                 population=20, crossover=0.7, mutation=0.3, max_generation=20,
                 hof_size=1, crowding_factor=1, plot_fitness=False,
                 compute_errors=False, show_contours=False,
                 plot_param1=None, plot_param2=None, sharing=False):

        self.like = like
        self.model = model
        self.outputname = outputname
        self.params = like.freeParameters()
        self.vpars = [p.value for p in self.params]
        self.npars = [p.name for p in self.params]
        self.sigma = sp.array([p.error for p in self.params])
        self.bounds = [p.bounds for p in self.params]
        self.cov = None

        print("Minimizing...", self.npars, "starting", self.vpars,  "with bounds", self.bounds)

        self.plot_fitness = plot_fitness
        self.compute_errors = compute_errors
        self.show_contours = show_contours
        self.plot_param1 = plot_param1
        self.plot_param2 = plot_param2

        # Genetic Algorithm constants:
        self.POPULATION_SIZE = population    # 10-20
        self.P_CROSSOVER = crossover       # probability for crossover
        self.P_MUTATION = mutation        # (try also 0.5) probability for mutating an individual
        self.MAX_GENERATIONS = max_generation    # 100- 300
        self.HALL_OF_FAME_SIZE = hof_size
        self.CROWDING_FACTOR = crowding_factor  # crowding factor for crossover and mutation

        self.RANDOM_SEED = 42        # set the random seed

        self.DIMENSIONS = len(self.params)        # number of dimensions
        self.BOUND_LOW, self.BOUND_UP = 0.0, 1.0  # boundaries for all dimensions

        self.sharing = sharing
        # sharing constants:
        if self.sharing:
            DISTANCE_THRESHOLD = 0.1
            SHARING_EXTENT = 5.0



    def main(self, pool):
        toolbox = self.GA()
        
        # using multiprocess
        toolbox.register("map", pool.map)

        # create initial population (generation 0):
        population = toolbox.populationCreator(n=self.POPULATION_SIZE)

        # prepare the statistics object:
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        
        if self.sharing:
            stats.register("max", np.max)
        else:
            stats.register("min", np.min)
        stats.register("avg", np.mean)

        # define the hall-of-fame object:
        hof = tools.HallOfFame(self.HALL_OF_FAME_SIZE)

        # perform the Genetic Algorithm flow with elitism:
        population, logbook, gens = elitism.eaSimpleWithElitism(population, toolbox, cxpb=self.P_CROSSOVER,
                                                              mutpb=self.P_MUTATION, ngen=self.MAX_GENERATIONS,
                                                              stats=stats, halloffame=hof, verbose=True,
                                                              outputname=self.outputname, bounds=self.bounds)

        # print info for best solution found:
        best = hof.items[0]
        print("-- Best Fitness = ", best.fitness.values[0])
        print("- Best solutions are:")
        best_params = [self.change_prior(i, x) for i, x in enumerate(best)]

        for i, x in enumerate(best_params):
            print("-- Best %s = "%self.params[i].name , x)
            # res.append("{}: {:.5f}".format(self.params[i].name, x))
        # res.append("Best Fitness: {:.5f}".format(best.fitness.values[0]))

        #for i in range(self.HALL_OF_FAME_SIZE):
        #    print(i, ": ", hof.items[i].fitness.values[0], " -> ", self.old_prior(i, hof.items[i]) )

        with open('{}.maxlike'.format(self.outputname), 'w') as f:
            np.savetxt(f, best_params, fmt='%.4e', delimiter=',')

        if self.plot_fitness:
            self.plotting(population, logbook, hof)

        if self.compute_errors:
            try:
                import numdifftools as nd
            except:
                sys.exit("*'error': Install numdifftools to compute errors.")

            hess = nd.Hessian(self.negloglike2, step=self.sigma*0.05)(best_params)
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

        # update with the final result
        #self.result(self.negloglike(self.res.x))

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
            plot_elipses(best_params, self.cov, idx_param1, idx_param2, param_Ltx1, param_Ltx2, ax=ax)
            plt.savefig('ga_plot.pdf')
            plt.show()

        return {'population': len(population), 'no_generations': gens, 'param_fit': best_params,
                'best_fitness': best.fitness.values[0], 'cov': self.cov, 'maxlike': best.fitness.values[0]}




    # helper function for creating random real numbers uniformly distributed within a given range [low, up]
    # it assumes that the range is the same for every dimension
    def randomFloat(self, low, up):
        return [random.uniform(l, u) for l, u in zip([low]*self.DIMENSIONS, [up]*self.DIMENSIONS)]


    # wraps the tools.selTournament() with fitness sharing
    # same signature as tools.selTournament()
    def selTournamentWithSharing(self, individuals, k, tournsize, fit_attr="fitness"):

        # get orig fitnesses:
        origFitnesses = [ind.fitness.values[0] for ind in individuals]

        # apply sharing to each individual:
        for i in range(len(individuals)):
            sharingSum = 1

            # iterate over all other individuals
            for j in range(len(individuals)):
                if i != j:
                    # calculate eucledean distance between individuals:
                    distance = math.sqrt(
                        ((individuals[i][0] - individuals[j][0]) ** 2) + ((individuals[i][1] - individuals[j][1]) ** 2))

                    if distance < DISTANCE_THRESHOLD:
                        sharingSum += (1 - distance / (SHARING_EXTENT * DISTANCE_THRESHOLD))

            # reduce fitness accordingly:
            individuals[i].fitness.values = origFitnesses[i] / sharingSum,

        # apply original tools.selTournament() using modified fitness:
        selected = tools.selTournament(individuals, k, tournsize, fit_attr)

        # retrieve original fitness:
        for i, ind in enumerate(individuals):
            ind.fitness.values = origFitnesses[i],

        return selected




    def GA(self):
        random.seed(self.RANDOM_SEED)
        toolbox = base.Toolbox()

        if  self.sharing:
            creator.create("FitnessMax", base.Fitness, weights=(1.0,))
            creator.create("Individual", list, fitness=creator.FitnessMax)
        else:
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
        if self.sharing:
            toolbox.register("select", selTournamentWithSharing, tournsize=2)
        else:
            toolbox.register("select", tools.selTournament, tournsize=2)
        
        toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=self.BOUND_LOW, \
                         up=self.BOUND_UP, eta=self.CROWDING_FACTOR)
        toolbox.register("mutate", tools.mutPolynomialBounded, low=self.BOUND_LOW, \
                         up=self.BOUND_UP, eta=self.CROWDING_FACTOR, indpb=1.0/self.DIMENSIONS)
                                                                    #indpb probability of each attribute to be mutated.
        return toolbox



    def plotting(self, pop, log, hof):
        # extract statistics
        gen, avg, min_, max_ = log.select("gen", "avg", "min", "max")
        #print(pop)
        plt.figure(figsize=(6, 6))
        plt.plot(gen, min_, label="minimum")

        plt.title("Fitness Evolution")
        plt.xlabel("Generation", fontsize=20)
        plt.ylabel("Fitness", fontsize=20)
        plt.legend(loc="upper right")
        #plt.savefig('GA_fitness.pdf')
        plt.show()



    def negloglike(self, x):
        for i, pars in enumerate(self.params):
            new_par = self.change_prior(i, x[i])
            pars.setValue(new_par)
        self.like.updateParams(self.params)
        loglike = self.like.loglike_wprior()

        if self.sharing:
            loglike = -loglike

        if sp.isnan(loglike):
            print('-1-'*10,loglike,'--'*10)
            return self.lastval+10
        else:
            self.lastval = -loglike
        return -loglike,


    def negloglike2(self, x):
        #original likelihood, standard prior
        for i, pars in enumerate(self.params):
            pars.setValue(x[i])
        self.like.updateParams(self.params)
        loglike = self.like.loglike_wprior()

        if self.sharing:
            loglike = -loglike

        if sp.isnan(loglike):
            return self.lastval+10
        else:
            self.lastval = -loglike
        return -loglike



    def change_prior(self, i, x):
        # Intput priors are on the range [0, 1]
        new_par = self.bounds[i][0] + (self.bounds[i][1] - self.bounds[i][0])*x
        return new_par

    def old_prior(self, i, x):
        old_par = (x - self.bounds[i][0])/(self.bounds[i][1] - self.bounds[i][0])
        return old_par


    def result(self, loglike):
        print ("------")
        print("Done.")
        print("Optimal loglike : ", loglike)
        return {'maxlike': loglike, 'param_fit': self.res.x, 'cov': self.cov}

if  __name__ == '__main__' :
    g= GA_deap()
    g.GA()
