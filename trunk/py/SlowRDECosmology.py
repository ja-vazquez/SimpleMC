#
# Jose Vazquez knows what this is... :)
# I do not.
#

import math as N
from LCDMCosmology import *

class SlowRDECosmology(LCDMCosmology):
    def __init__(self, varyw=True, varyOk=True):
        ## two parameters: Om and h

        self.varyw=varyw
        self.varyOk=varyOk

        self.Ok=Ok_par.value   
        self.dw0=dw_par.value
      
        LCDMCosmology.__init__(self)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyw): l.append(dw_par)
        if (self.varyOk): l.append(Ok_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="dw":
                self.dw0=p.value
            elif p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False
        return True

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        onepz=1.0/a
        NuContrib=self.NuDensity.rho(a)/self.h**2
        Ode = 1.0-self.Om-self.Ok
        rhow= (onepz**3.0/(self.Om*onepz**3.0+Ode))**(self.dw0/Ode)
        return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+NuContrib+Ode*rhow)

