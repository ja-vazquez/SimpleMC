## This is a CDM cosmology with w

from LCDMCosmology import *

class wLCDMCosmology(LCDMCosmology):
    def __init__(self):
        ## two parameters: Om and h
        self.w=-1
        LCDMCosmology.__init__(self)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        return [w_par]+LCDMCosmology.freeParameters(self)


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="w":
                self.w=p.value
        return True

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om)*a**(-3*(1.0+self.w)))


