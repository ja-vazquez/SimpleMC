## This is LCDM cosmology with optional
## curvature which you can set up with 
## setVaryOk()

from numpy import linspace
from scipy.interpolate import InterpolatedUnivariateSpline
from LCDMCosmology import *

class SplineLCDMCosmology(LCDMCosmology):
    def __init__(self):
        ## two parameters: Om and h
        self.w=-0.5
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

    def Spline(self, a):
        z=1.0/a-1.0
        x=[0,1,1.5,2.34,3]
        y=[1.0,1.01,0.9,2.0, 1.01]
        s= InterpolatedUnivariateSpline(x, y) 
        #xs = linspace(0,3,100) 
        ys = s(z)
        return ys

    def Rho_de(self,a):
        z=1.0/a-1.0
        if (z > 4):
          return 1.0
        else:
          #return (1.0-self.Om)*a**(-3*(1.0+self.w))
          return self.Spline(a)

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om)*self.Rho_de(a))


