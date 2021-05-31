from simplemc import logger
from simplemc.models.LCDMCosmology import LCDMCosmology



class SimpleModel:
    """
    This is a generic model

    """
    def __init__(self, parameters, function):
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


