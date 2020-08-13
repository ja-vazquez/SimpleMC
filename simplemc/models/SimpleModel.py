from simplemc import logger


#TODO There's SimpleModel, SimpleCosmology and Generic, unify them


class SimpleModel:
    def __init__(self, parameters, function):
        """
        Parameters:
            parameters: List of Parameter class objects
            function: This function recivies a vector in the parameter space
                        and returns its evaluation in the function.
                        Example:
                            def straight_line(parameters, x):
                                m, b = parameters
                                return m*x +b
            cosmological: boolean value. If True, then the LCDMCosmology
                          class will be the base of the model.
        """
        self.parameters = parameters
        self.function = function
        SimpleModel.updateParams(self, [])

    # my free params (parameters/priors see ParamDefs.py)
    def freeParameters(self):
        return self.parameters

    def printFreeParameters(self):
        logger.info("Free parameters and its bounds:")
        self.printParameters(self.freeParameters())

    def printParameters(self, params):
        for p in params:
            logger.info('{} = {} +/- {}'.format(p.name, p.value, p.error))

    def updateParams(self, pars):
        return True

    def genericModel(self, x):
        values = []
        for param in self.parameters:
            values.append(param.value)
        return self.function(values, x)

    def prior_loglike(self):
        return 0


