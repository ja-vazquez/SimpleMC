
import numpy as np
import scipy.constants as cte
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Importamos libreria de algoritmos evolutivos
from deap import base, creator, tools, algorithms

# importamos modulo independiente para acoplar a DEAP, adjunto en la carpeta fuente.
import elitism

# We import an independent module to implement elitism in the GA.
import random


class GA_deap:
    def __init__(self, like=None, model=None):

        """
        self.like   = like
        self.model  = model
        self.params = like.freeParameters()
        self.vpars  = [p.value for p in self.params]
        self.sigma  = sp.array([p.error for p in self.params])
        bounds = [p.bounds for p in self.params]
        print("Minimizing...", self.vpars, "with bounds", bounds)
        """
        pass

    def GA(self):
        # problem constants:
        DIMENSIONS = 3                  # number of dimensions
        BOUND_LOW, BOUND_UP = 0.0, 1.0  # boundaries for all dimensions


        # Genetic Algorithm constants:
        POPULATION_SIZE = 20    # 10-20
        P_CROSSOVER = 0.7       # probability for crossover
        P_MUTATION = 0.3        # (try also 0.5) probability for mutating an individual
        MAX_GENERATIONS = 15    # 100- 300
        HALL_OF_FAME_SIZE = 1
        CROWDING_FACTOR = 20.0  # crowding factor for crossover and mutation

        # set the random seed:
        RANDOM_SEED = 42
        random.seed(RANDOM_SEED)

        toolbox = base.Toolbox()

        # define a single objective, minimizing fitness strategy:
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))

        # create the Individual class based on list:
        creator.create("Individual", list, fitness=creator.FitnessMin)

        # helper function for creating random real numbers uniformly distributed within a given range [low, up]
        # it assumes that the range is the same for every dimension
        def randomFloat(low, up):
            return [random.uniform(l, u) for l, u in zip([low] * DIMENSIONS, [up] * DIMENSIONS)]

        # create an operator that randomly returns a float in the desired range and dimension:
        toolbox.register("attrFloat", randomFloat, BOUND_LOW, BOUND_UP)

        # create the individual operator to fill up an Individual instance:
        toolbox.register("individualCreator", tools.initIterate, creator.Individual, toolbox.attrFloat)


        print('working')


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


if  __name__ == '__main__' :
    g= GA_deap()
    g.GA()