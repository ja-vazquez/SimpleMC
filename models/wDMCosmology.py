## This is a CDM cosmology with w

from LCDMCosmology import *

class wDMCosmology(LCDMCosmology):
    def __init__(self):
        ## two parameters: Om and h
        self.wDM=0
        LCDMCosmology.__init__(self)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        return [wDM_par]+LCDMCosmology.freeParameters(self)


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="wDM":
                self.wDM=p.value
        return True

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
	self.Odm = self.Ocb -self.Obh2/(self.h**2)
	self.Ob = self.Ocb -self.Odm 
        return (self.Ob/a**3+ self.Odm*a**(-3*(1.0+self.wDM))+self.Omrad/a**4+NuContrib+(1.0-self.Om))


