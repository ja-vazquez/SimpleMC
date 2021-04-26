from .Individual import Individual
import numpy as np
import copy
import sys
from simplemc import logger
try:
    import pandas as pd
except:
    import warnings
    warnings.warn("Please install Pandas library if you want to use SimpleGenetic.")

class Population:
    def __init__(self, n_individuals, n_variables, lower_bounds=None,
                 upper_bounds=None):
        """
        It creates several objects of the Individual class and generates a Population.

        Parameters
        ----------
        n_individuals : int
            Number of initial individuals

        n_variables : int
            Number of variables or free parameters

        lower_bounds : list
            List with the lower bounds for each free parameter.

        upper_bounds : list
            List with the upper bounds for each free parameter.
        """

        self.n_individuals = n_individuals
        self.n_variables = n_variables
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds
        self.individuals = []
        self.optimized = False
        self.iter_optimization = None
        self.best_individual = None
        self.best_fitness = None
        self.best_function_value = None
        self.best_value_variables = None
        self.historical_individuals = []
        self.historical_best_value_variables = []
        self.historical_best_fitness = []
        self.historical_best_function_value = []
        self.difference_abs = []
        self.results_df = None
        self.fitness_optimal = None
        self.value_variables_optimal = None
        self.function_value_optimal = None

        if self.lower_bounds is not None \
                and not isinstance(self.lower_bounds, np.ndarray):
            self.lower_bounds = np.array(self.lower_bounds)

        if self.upper_bounds is not None and not isinstance(self.upper_bounds, np.ndarray):
            self.upper_bounds = np.array(self.upper_bounds)

        for i in np.arange(n_individuals):
            individual_i = Individual(
                n_variables=self.n_variables,
                lower_bounds=self.lower_bounds,
                upper_bounds=self.upper_bounds)
            self.individuals.append(individual_i)



    def eval_population(self, target_function, optimization):
        """

        Parameters
        ----------
        target_function : method
            Function to optimize. Usually (in simplemc context) the likelihood.

        optimization : str
            {"maximize", "minimize}
            Default: maximize

        Returns
        -------

        """
        for i in np.arange(self.n_individuals):
            self.individuals[i].calculate_fitness(
                target_function=target_function,
                optimization=optimization
            )

        self.best_individual = copy.deepcopy(self.individuals[0])
        for i in np.arange(self.n_individuals):
            if self.individuals[i].fitness > self.best_individual.fitness:
                self.best_individual = copy.deepcopy(self.individuals[i])

        self.best_fitness = self.best_individual.fitness
        self.best_value_variables = self.best_individual.value_variables
        self.best_function_value = self.best_individual.function_value



    def cross_individuals(self, parental_1, parental_2):
        """
        Method to cross two differents individuals
        Parameters
        ----------
        parental_1
        parental_2

        Returns
        -------

        """

        if parental_1 not in np.arange(self.n_individuals):
            raise Exception(
                "index of parental_1 should be between 0 and "
                "the number of individuals of the population."
            )
        if parental_2 not in np.arange(self.n_individuals):
            raise Exception(
                "index of parental_2 should be between 0 and "
                "the number of individuals of the population."
            )

        parental_1 = self.individuals[parental_1]
        parental_2 = self.individuals[parental_2]

        offspring = copy.deepcopy(parental_1)
        offspring.value_variables = np.repeat(None, offspring.n_variables)
        offspring.fitness = None

        inheritance_parent_1 = np.random.choice(
            a=[True, False],
            size=offspring.n_variables,
            p=[0.5, 0.5],
            replace=True
        )
        inheritance_parent_2 = np.logical_not(inheritance_parent_1)

        offspring.value_variables[inheritance_parent_1] \
            = parental_1.value_variables[inheritance_parent_1]

        offspring.value_variables[inheritance_parent_2] \
            = parental_2.value_variables[inheritance_parent_2]

        offspring = copy.deepcopy(offspring)

        return (offspring)



    def select_individual(self, n, return_indexs=True,
                          method_selection="tournament"):
        """
        This method selects an individual given a selection method.

        Parameters
        ----------
        n : int
        return_indexs : bool
        method_selection : str

        Returns
        -------

        """

        if method_selection not in ["roulette", "rank", "tournament"]:
            raise Exception(
                "Selection method should be roulette, rank o tournament"
            )

        array_fitness = np.repeat(None, self.n_individuals)
        for i in np.arange(self.n_individuals):
            array_fitness[i] = copy.copy(self.individuals[i].fitness)

        if method_selection == "roulette":
            selection_probability = array_fitness / np.sum(array_fitness)
            ind_selected = np.random.choice(
                a=np.arange(self.n_individuals),
                size=n,
                p=list(selection_probability),
                replace=True
            )
        elif method_selection == "rank":
            order = np.flip(np.argsort(a=array_fitness) + 1)
            ranks = np.argsort(order) + 1
            selection_probability = 1 / ranks
            selection_probability = selection_probability / np.sum(selection_probability)
            ind_selected = np.random.choice(
                a=np.arange(self.n_individuals),
                size=n,
                p=list(selection_probability),
                replace=True
            )
        elif method_selection == "tournament":
            ind_selected = np.repeat(None, n)
            for i in np.arange(n):
                # Random selection of a pair of Individuals
                candidates_a = np.random.choice(
                    a=np.arange(self.n_individuals),
                    size=2,
                    replace=False
                )
                candidates_b = np.random.choice(
                    a=np.arange(self.n_individuals),
                    size=2,
                    replace=False
                )
                # Selection of best fitness of each pair.
                if array_fitness[candidates_a[0]] > array_fitness[candidates_a[1]]:
                    winner_a = candidates_a[0]
                else:
                    winner_a = candidates_a[1]

                if array_fitness[candidates_b[0]] > array_fitness[candidates_b[1]]:
                    winner_b = candidates_b[0]
                else:
                    winner_b = candidates_b[1]

                # For each couple, comparison of winners
                if array_fitness[winner_a] > array_fitness[winner_b]:
                    ind_final = winner_a
                else:
                    ind_final = winner_b

                ind_selected[i] = ind_final

        if (return_indexs):
            return (ind_selected)
        else:
            if n == 1:
                return (copy.deepcopy(self.individuals[int(ind_selected)]))
            if n > 1:
                return (
                    [copy.deepcopy(self.individuals[i]) for i in ind_selected]
                )



    def create_new_generation(self, method_selection="tournament",
                              elitism=0.1, prob_mut=0.01,
                              distribution="uniform",
                              media_distribution=1, sd_distribution=1,
                              min_distribution=-1, max_distribution=1):

        news_individuals = []

        if elitism > 0:
            n_elitism = int(np.ceil(self.n_individuals * elitism))
            array_fitness = np.repeat(None, self.n_individuals)
            for i in np.arange(self.n_individuals):
                array_fitness[i] = copy.copy(self.individuals[i].fitness)
            rank = np.flip(np.argsort(array_fitness))
            elite = [copy.deepcopy(self.individuals[i]) for i in rank[:n_elitism]]
            news_individuals = news_individuals + elite
        else:
            n_elitism = 0

        for i in np.arange(self.n_individuals - n_elitism):
            # select parents
            index_parents = self.select_individual(
                n=2,
                return_indexs=True,
                method_selection=method_selection)

            offspring = self.cross_individuals(
                parental_1=index_parents[0],
                parental_2=index_parents[1])

            offspring.mutate(
                prob_mut=prob_mut,
                distribution=distribution,
                min_distribution=min_distribution,
                max_distribution=max_distribution)

            news_individuals = news_individuals + [offspring]

        self.individuals = copy.deepcopy(news_individuals)
        self.best_individual = None
        self.best_fitness = None
        self.best_value_variables = None
        self.best_function_value = None

        # New generation created! Known number of new individuals.



    def optimize(self, target_function, optimization, n_generations=50,
                 method_selection="tournament", elitism=0.1, prob_mut=0.01,
                 distribution="uniform", media_distribution=1,
                 sd_distribution=1, min_distribution=-1, max_distribution=1,
                 stopping_early=False, rounds_stopping=3,
                 tolerance_stopping=0.1, outputname="geneticOutput"):
        """
        This is the crucial method of the genetic modules.
        Within it, a target function is optimized.

        Parameters
        ----------
        target_function : method
            Usually the likelihood function.

        optimization : str
            Options {"maximize", "minimize"}. Default is "maximize"

        n_generations : int
            Number of generations to evolve the original population

        method_selection : str
            Type of selection method {"tournament", "rank", "roulette"}.
            Default: "tournament"

        elitism : float
            Value of elitism.
            Default: 0.01.

        prob_mut : float
            Probability of mutation.
            Default: 0.4

        distribution : str
            {"uniform", "gaussian", "random"}
            Default: uniform

        media_distribution : float
            Media value for gaussian distributions

        sd_distribution : float
            Standard deviation for gaussian distributions
            Default: 1.0

        min_distribution : float
            Minimum value for uniform distributions
            Default: -1.0

        max_distribution : float
            Maximum value for uniform distributions
            Default: 1.0

        stopping_early : bool
            It needs a not None value for "rounds_stopping" and "tolerance_stopping".
            Default: True

        rounds_stopping : int
            Rounds to consider to stopping early with the tolerance_stopping value.

        tolerance_stopping : float
            Value to stopping early criteria. This value is the difference between the
            best fit for the latest rounds_stopping generations.

        outputname : str
            Name for the output text file.
            Default: "geneticOutput"

        """

        if stopping_early and (rounds_stopping is None or tolerance_stopping is None):
            raise Exception("To activate Stopping early it is necessary to indicate a",
                            "value of rounds_stopping and tolerance_stopping.")

        f = open("{}.txt".format(outputname), "+w")
        f.write("# Generation, best_fitness\n")
        for i in np.arange(n_generations):

            self.eval_population(
                target_function=target_function,
                optimization=optimization)

            self.historical_individuals.append(copy.deepcopy(self.individuals))
            self.historical_best_fitness.append(copy.deepcopy(self.best_fitness))
            self.historical_best_value_variables.append(
                copy.deepcopy(self.best_value_variables)
            )
            self.historical_best_function_value.append(
                copy.deepcopy(self.best_function_value)
            )

            if i == 0:
                self.difference_abs.append(None)
            else:
                difference = abs(self.historical_best_fitness[i] \
                                 - self.historical_best_fitness[i - 1])
                self.difference_abs.append(difference)
            f.write("{} {} \n".format(i, self.best_fitness))
            sys.stdout.write("\r{}it | best fitness: {:.4f} | best individual: {}".format(i + 1, self.best_fitness,
                                                                                          self.best_value_variables))
            sys.stdout.flush()
            if stopping_early and i > rounds_stopping:
                latest_n = np.array(self.difference_abs[-(rounds_stopping):])
                if all(latest_n < tolerance_stopping):
                    print("Algorithm stopped in the {} generation"
                          "for lack of absolute minimum change of {}"
                          "while {} consecutive generations".format(i + 1, tolerance_stopping,
                                                                    rounds_stopping))

                    break

            self.create_new_generation(
                method_selection=method_selection,
                elitism=elitism,
                prob_mut=prob_mut,
                distribution=distribution)

        f.close()
        self.optimized = True
        self.iter_optimization = i

        index_value_optimal = np.argmax(np.array(self.historical_best_fitness))
        self.fitness_optimal = self.historical_best_fitness[index_value_optimal]
        self.function_value_optimal = self \
            .historical_best_function_value[index_value_optimal]
        self.value_variables_optimal = self \
            .historical_best_value_variables[index_value_optimal]

        self.results_df = pd.DataFrame(
            {
                "best_fitness": self.historical_best_fitness,
                "best_function_value": self.historical_best_fitness,
                "best_value_variables": self.historical_best_value_variables,
                "difference_abs": self.difference_abs
            }
        )
        self.results_df["generation"] = self.results_df.index

        logger.info("\nSummary\n------\nNumber of generations: {}\n"
                    "Optimal value of variables: {}\n"
                    "Target function value: {}\n".format(self.iter_optimization + 1,
                                                         self.value_variables_optimal,
                                                         self.function_value_optimal))

        return ["Number of generations: {}\n"
                "Optimal value of variables: {}\n"
                "Target function value: {}\n".format(self.iter_optimization + 1,
                                                     self.value_variables_optimal,
                                                     self.function_value_optimal)]
