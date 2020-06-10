

from simplemc.models.LCDMCosmology import LCDMCosmology
from  simplemc.cosmo.paramDefs import wq_par, Oq_par
import math as N

class DGPCDMCosmology(LCDMCosmology):
    def __init__(self):
        """
        Simple Class for Brane Cosmology,
        here flat-DGP model
        Returns
        -------

        """
        self.varywq = True
        self.varyOq = True

        self.wq = wq_par.value
        self.Oq = Oq_par.value
        LCDMCosmology.__init__(self)


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varywq):  l.append(wq_par)
        if (self.varyOq): l.append(Oq_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "wq":
                self.wq = p.value
            elif p.name == "Oq":
                self.Oq = p.value
        return True


    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        NuContrib = self.NuDensity.rho(a)/self.h**2
        if self.Ocb + self.Oq <= 1:
            return self.Ocb/a**3
        else:
            Orc = (0.5*(self.Ocb + self.Oq -1 ))**2
            return  (N.sqrt(self.Omrad/a**4 +NuContrib + self.Ocb/a**3 + Orc + self.Oq/a**(3*(self.wq+1))) - N.sqrt(Orc) )**2