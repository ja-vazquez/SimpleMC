import numpy as np
import random
import warnings


class Individual:
    def __init__(self, n_variables, lower_lims=None, upper_lims=None):
        """
        This class generates the individuals within the bounds of the parameter space.

        Parameters
        -----------
        n_variables : int
            Number of variables or free parameters

        lower_lims : list
            List with the lower bounds
        """
        self.n_variables = n_variables
        self.lower_lims = lower_lims
        self.upper_lims = upper_lims
        self.value_variables = np.repeat(None, n_variables)
        self.fitness = None
        self.function_value = None

        if self.lower_lims is not None and not isinstance(self.lower_lims, np.ndarray):
            self.lower_lims = np.array(self.lower_lims)

        if self.upper_lims is not None \
                and not isinstance(self.upper_lims, np.ndarray):
            self.upper_lims = np.array(self.upper_lims)

        if self.lower_lims is not None and len(self.lower_lims) != self.n_variables:
            raise Exception(
                "lower_lims must have a value for each variable." +
                "If you do not want to limit any variable, use None." +
                "Ex: lower_lims = [10, None, 5]"
            )
        elif self.upper_lims is not None \
                and len(self.upper_lims) != self.n_variables:
            raise Exception(
                "upper_lims must have a value for each variable." +
                "If you do not want to limit any variable, use None." +
                "Ex: lower_lims = [10, None, 5]"
            )
        elif (self.lower_lims is None) or (self.upper_lims is None):
            warnings.warn(
                "It is highly recommended to indicate" +
                " the limits within which the solution of each variable." +
                "By default: [-10^3, 10^3]."
            )
        elif any(np.concatenate((self.lower_lims, self.upper_lims)) == None):
            warnings.warn(
                "By default the limits are: [-10^3, 10^3]."
            )

        if self.lower_lims is None:
            self.lower_lims = np.repeat(-10 ** 3, self.n_variables)

        if self.upper_lims is None:
            self.upper_lims = np.repeat(+10 ** 3, self.n_variables)

        if self.lower_lims is not None:
            self.lower_lims[self.lower_lims == None] = -10 ** 3

        if self.upper_lims is not None:
            self.upper_lims[self.upper_lims == None] = +10 ** 3

        for i in np.arange(self.n_variables):
            self.value_variables[i] = random.uniform(
                self.lower_lims[i],
                self.upper_lims[i]
            )

    # def __repr__(self):
    #     """
    #     Info for print individual object.
    #     """
    #     text = ("Individual \n --------- \n  Variables values:"
    #              "{} \n Target function value: {} \n Fitness:"
    #              "{} \n lower limits for each variable: {}"
    #              "upper bounds for each variable:  {} \n").format(
    #              self.value_variables, self.function_value,
    #              self.fitness, self.lower_lims, self.upper_lims)
    #
    #     return(text)

    def calculate_fitness(self, target_function, optimization):
        if not optimization in ["maximize", "minimize"]:
            raise Exception(
                "Arg should be: 'maximize' or 'minimize'"
            )

        self.function_value = target_function(*self.value_variables)
        if optimization == "maximize":
            self.fitness = self.function_value
        elif optimization == "minimize":
            self.fitness = -self.function_value

    def mutate(self, prob_mut=0.01, distribution="uniform", media_distribution=1,
              sd_distribution=1, min_distribution=-1, max_distribution=1):
        if not distribution in ["normal", "uniform", "random"]:
            raise Exception(
                "Arg should be: 'normal', 'uniform' or 'random'"
            )

        pos_mutated = np.random.uniform(
            low=0,
            high=1,
            size=self.n_variables
        )
        pos_mutated = pos_mutated < prob_mut

        if distribution in ["normal", "uniform"]:
            if distribution == "normal":
                factor_mut = np.random.normal(
                    loc=media_distribution,
                    scale=sd_distribution,
                    size=np.sum(pos_mutated)
                )
            if distribution == "uniform":
                factor_mut = np.random.uniform(
                    low=min_distribution,
                    high=max_distribution,
                    size=np.sum(pos_mutated)
                )
            self.value_variables[pos_mutated] = \
                self.value_variables[pos_mutated] + factor_mut

            for i in np.flatnonzero(pos_mutated):
                if self.value_variables[i] < self.lower_lims[i]:
                    self.value_variables[i] = self.lower_lims[i]
                if self.value_variables[i] > self.upper_lims[i]:
                    self.value_variables[i] = self.upper_lims[i]

        if distribution == "random":
            for i in np.flatnonzero(pos_mutated):
                self.value_variables[i] = random.uniform(
                    self.lower_lims[i],
                    self.upper_lims[i]
                )

        self.fitness = None
        self.function_value = None

