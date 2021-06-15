from simplemc import logger
from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
import math as N
import numpy as np
import sys

class SimpleModel:
    """
        This is a generic model

        Parameters
        ----------
        parameters : list
            List of Parameter objects
        function : function
            model or function. It should be in terms of the parameters list.
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
        l = []
        for p in params:
            logger.info('{} = {} +/- {}'.format(p.name, p.value, p.error))
            l.append("{}: {} = +/- {}".format(p.name, p.value, p.error))
        return l

    def updateParams(self, pars):
        return True

    def genericModel(self, x):
        values = []
        for param in self.parameters:
            values.append(param.value)
        return self.function(values, x)

    def prior_loglike(self):
        return 0


class SimpleCosmoModel(LCDMCosmology):
    def __init__(self,  extra_params=None, RHSquared=None):
        """
        This is a simple cosmological model based on slightly deviations of LCDMCosmology
        class, RHSquared must to have an analytical form.

        Parameters
        ----------
        extra_params : list
            List of Parameter objects with the extra parameters.
            LCDM already have:
                        h, NuContrib, Ocb, Omrad and Om
        RHSquared : str
            model or function. It should be in terms of scale factor a
            Example: 'Ocb/a**3+Omrad/a**4+NuContrib+(1.0-Om-newparameter)'
                      where newparameter = Parameter('newparameter', value, step, (b1,b2), '$Latex_name')
        """
        self.extra_params = extra_params
        self.RHSquared = RHSquared
        LCDMCosmology.__init__(self)

    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        for parameter in self.extra_params:
            l.append(parameter)
        return l

    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        return True

    # this is relative hsquared as a function of a
    # i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        h = self.h
        NuContrib = self.NuDensity.rho(a) / h ** 2
        Ocb = self.Ocb
        Omrad = self.Omrad
        Om = self.Om
        for param in self.extra_params:
            exec("%s = %f" % (param.name, param.value))
        if self.RHSquared:
            return eval(self.RHSquared)
        else:
            sys.exit('Please set a string with the RHSquared (H(z)^2/H(z=0)^2) '
                     'expression in terms of a (scale factor)')