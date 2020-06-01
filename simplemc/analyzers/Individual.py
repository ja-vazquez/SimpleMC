import numpy as np
import random
import warnings


class Individual:
    def __init__(self, n_variables, limits_inf=None, limits_sup=None,
                 verbose=False):
        self.n_variables = n_variables
        self.limits_inf = limits_inf
        self.limits_sup = limits_sup
        self.value_variables = np.repeat(None, n_variables)
        self.fitness = None
        self.function_value = None

        if self.limits_inf is not None and not isinstance(self.limits_inf, np.ndarray):
            self.limits_inf = np.array(self.limits_inf)

        if self.limits_sup is not None \
                and not isinstance(self.limits_sup, np.ndarray):
            self.limits_sup = np.array(self.limits_sup)

        if self.limits_inf is not None and len(self.limits_inf) != self.n_variables:
            raise Exception(
                "limits_inf must have a value for each variable." +
                "If you do not want to limit any variable, use None." +
                "Ex: limits_inf = [10, None, 5]"
            )
        elif self.limits_sup is not None \
                and len(self.limits_sup) != self.n_variables:
            raise Exception(
                "limits_sup must have a value for each variable." +
                "If you do not want to limit any variable, use None." +
                "Ex: limits_inf = [10, None, 5]"
            )
        elif (self.limits_inf is None) or (self.limits_sup is None):
            warnings.warn(
                "It is highly recommended to indicate" +
                " the limits within which the solution of each variable." +
                "By default: [-10^3, 10^3]."
            )
        elif any(np.concatenate((self.limits_inf, self.limits_sup)) == None):
            warnings.warn(
                "By default the limits are: [-10^3, 10^3]."
            )

        if self.limits_inf is None:
            self.limits_inf = np.repeat(-10 ** 3, self.n_variables)

        if self.limits_sup is None:
            self.limits_sup = np.repeat(+10 ** 3, self.n_variables)

        if self.limits_inf is not None:
            self.limits_inf[self.limits_inf == None] = -10 ** 3

        if self.limits_sup is not None:
            self.limits_sup[self.limits_sup == None] = +10 ** 3

        for i in np.arange(self.n_variables):
            self.value_variables[i] = random.uniform(
                self.limits_inf[i],
                self.limits_sup[i]
            )
        if verbose:
            print("New individal created")
            print("----------------------")
            print("Variables values: {}".format(self.value_variables))
            print("Target function value: {}".format(self.function_value))
            print("Fitness: {}".format(self.fitness))
            print("Lower bounds: {}".format(self.limits_inf))
            print("Upper limits for each variable: {}".format(self.limits_sup))

    def __repr__(self):
        """
        Info for print individual object.
        """
        text = ("Individual \n --------- \n  Variables values:"
                 "{} \n Target function value: {} \n Fitness:"
                 "{} \n lower limits for each variable: {}"
                 "upper bounds for each variable:  {} \n").format(
                 self.value_variables, self.function_value,
                 self.fitness, self.limits_inf, self.limits_sup)

        return(text)

    def calculate_fitness(self, target_function, optimization, verbose=False):
        if not optimization in ["maximize", "minimize"]:
            raise Exception(
                "Arg should be: 'maximize' or 'minimize'"
            )

        self.function_value = target_function(*self.value_variables)
        if optimization == "maximize":
            self.fitness = self.function_value
        elif optimization == "minimize":
            self.fitness = -self.function_value

        if verbose:
            print("The individual has been evaluated")
            print("-----------------------------")
            print("Value target function: {}".format(self.function_value))
            print("Fitness: {}".format(self.fitness))
            print("")

    def mutate(self, prob_mut=0.01, distribution="uniform", media_distribution=1,
              sd_distribution=1, min_distribution=-1, max_distribution=1,
              verbose=False):
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
                if self.value_variables[i] < self.limits_inf[i]:
                    self.value_variables[i] = self.limits_inf[i]
                if self.value_variables[i] > self.limits_sup[i]:
                    self.value_variables[i] = self.limits_sup[i]

        if distribution == "random":
            for i in np.flatnonzero(pos_mutated):
                self.value_variables[i] = random.uniform(
                    self.limits_inf[i],
                    self.limits_sup[i]
                )

        self.fitness = None
        self.function_value = None

        if verbose:
            print("Individual has been mutated")
            print("Total mutuations: {}".format(np.sum(pos_mutated)))
            print("Values of variables: {}".format(self.value_variables))
