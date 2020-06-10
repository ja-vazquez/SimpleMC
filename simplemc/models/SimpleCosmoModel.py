from .LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
from simplemc.cosmo.paramDefs import Ok_par, w_par, wa_par
import math as N


#TODO There's SimpleModel, SimpleCosmology and Generic, unify them


class SimpleCosmoModel(LCDMCosmology):
    def __init__(self):
        # Create a list with your parameters,
        # pull it from paramDefs or create them
        self.miparam = Parameter("nombre", 0.5, 0.005, (0.01, 1), "\LaTeX")
        # with the Parameter class
        # self.nuevoparm = Parameter("nombre", 0.5)
        self.parameters = [Ok_par, w_par, wa_par]
        self.Ok = Ok_par.value
        self.w0 = w_par.value
        self.wa = wa_par.value

        LCDMCosmology.__init__(self)
    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        for parameter in self.parameters:
            l.append(parameter)
        return l

    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        return True

    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    # 2) Write your own RHSquared function
    # from LCDM you have self.h, self.Omrad, self.Om, self.Ocb
    def RHSquared_a(self, a):
        NuContrib = self.NuDensity.rho(a) / self.h ** 2
        rhow = a ** (-3 * (1.0 + self.w0 + self.wa)) * N.exp(-3 * self.wa * (1 - a))
        return (self.Ocb / a ** 3 + self.Ok / a ** 2 + self.Omrad / a ** 4 + NuContrib + (
                    1.0 - self.Om - self.Ok) * rhow)
