# This is a CDM cosmology with constant eos w for DE

from LCDMCosmology import LCDMCosmology
from ParamDefs import w_par


class wCDMCosmology(LCDMCosmology):
    def __init__(self, varyw=True):
        # two parameters: Om and h

        self.varyw = varyw
        self.w = w_par.value
        LCDMCosmology.__init__(self)


    # my free parameters. We add w on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varyw): l.append(w_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "w":
                self.w = p.value
        return True


    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        NuContrib = self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om)*a**(-3*(1.0+self.w)))
