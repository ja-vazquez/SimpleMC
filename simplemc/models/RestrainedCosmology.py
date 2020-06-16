

import math as N
from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import *

#Scalar field emulator via anisotropically deformed vacuum energy: Application to dark energy
#https://arxiv.org/abs/2004.14863

class RestrainedCosmology(LCDMCosmology):
    def __init__(self, varyOk=False, varyweff = True, varywcpl = False):
        """
        generalization of the usual vacuum energy, called `deformed vacuum energy',
        which yields anisotropic pressure whilst preserving zero inertial mass density
        Parameters
        ----------
        varyOk
        varyweff
        varywcpl

        Returns
        -------

        """

        self.varyweff = varyweff
        self.varywcpl = varywcpl
        self.varyOk   = varyOk

        self.Ok   = Ok_par.value
        self.weff = weff_par.value
        self.wcpl = wcpl_par.value
        LCDMCosmology.__init__(self)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyweff): l.append(weff_par)
        if (self.varywcpl): l.append(wcpl_par)
        if (self.varyOk): l.append(Ok_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="weff":
                self.weff=p.value
            if p.name=="wcpl":
                self.wcpl=p.value
            elif p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False
        return True

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        #NuContrib=self.NuDensity.rho(a)/self.h**2

        Omegasig = (1- self.Ocb- self.Omrad)*(1 + self.weff)*0.5
        Omegadv  = (1- self.Ocb- self.Omrad)*(1 - self.weff)*0.5
        Omegaeff = Omegasig + Omegadv

        rhow     = Omegaeff*a**(-3*(1.0+self.weff + self.wcpl))*N.exp(-3*self.wcpl*(1-a))

        return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4 + rhow)

