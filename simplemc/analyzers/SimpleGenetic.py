from .Population import Population


class SimpleGenetic:
    def __init__(self, target_function, n_variables, bounds, n_individuals=50,
                optimization="maximize",
                n_generations=250, method_selection="tournament", elitism=0.01,
                prob_mut=0.1, distribution="uniform", media_distribution=1,
                sd_distribution=1, min_distribution=-1, max_distribution=1,
                stopping_early=True, rounds_stopping=500, tolerance_stopping=0.1,
                outputname="geneticOutput"):

        self.target_function = target_function
        # These bounds are a list where every input is the limit of a param
        bounds = bounds

        self.lower_bounds = []
        self.upper_bounds = []

        for bound in bounds:
            self.lower_bounds.append(bound[0])
            self.upper_bounds.append(bound[1])

        self.n_individuals = n_individuals
        self.n_variables = n_variables
        
        self.distribution = distribution
        self.elitism = elitism
        self.max_distribution = max_distribution
        self.media_distribution = media_distribution
        self.method_selection = method_selection
        self.min_distribution = min_distribution
        self.n_generations = n_generations
        self.optimization = optimization
        self.stopping_early = stopping_early
        self.prob_mut = prob_mut
        self.rounds_stopping = rounds_stopping
        self.sd_distribution = sd_distribution
        self.tolerance_stopping  = tolerance_stopping
        self.outputname = outputname

        self.optimize()

    def optimize(self):
        population = Population(n_individuals=self.n_individuals,
                n_variables=self.n_variables,
                lower_bounds=self.lower_bounds,
                upper_bounds=self.upper_bounds)

        o = population.optimize(target_function=self.target_function,
                            optimization=self.optimization,
                            n_generations=self.n_generations,
                            method_selection=self.method_selection,
                            elitism=self.elitism,
                            prob_mut=self.prob_mut,
                            distribution=self.distribution,
                            media_distribution=self.media_distribution,
                            sd_distribution=self.sd_distribution,
                            min_distribution=self.min_distribution,
                            max_distribution=self.max_distribution,
                            stopping_early=self.stopping_early,
                            rounds_stopping=self.rounds_stopping,
                            outputname=self.outputname)
        return o