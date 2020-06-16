

from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import *

#Observational constraints on the free parameters of an
# interacting Bose-Einstein gas as a dark-energy model
#https://arxiv.org/abs/1802.07232

class IBEGCosmology(LCDMCosmology):
    def __init__(self, varyOk=False, varyOi0=True, varyxx=True):
        """
        Dark energy is modelled by a Bose-Einstein gas of particles
        with an attractive interaction.
        Parameters
        ----------
        varyOk
        varyOi0
        varyxx

        Returns
        -------

        """

        self.varyOk  = varyOk
        self.varyOi0 = varyOi0
        self.varyxx  = varyxx

        self.Ok  = Ok_par.value
        self.Oi0 = Oi0_par.value
        self.xx  = xx_par.value
        LCDMCosmology.__init__(self)


    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyOk): l.append(Ok_par)
        if (self.varyOi0): l.append(Oi0_par)
        if (self.varyxx): l.append(xx_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="Oi0":
                self.Oi0=p.value
            elif p.name=="xx":
                self.xx=p.value
            elif p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False
        return True


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        Oc0 = (1.0-self.Om-self.Ok - self.Oi0/(1 - 2*self.xx))*(2- 5.*self.xx)/2.
        rhow = 2.*Oc0*a**(5*self.xx-5.)/(2- 5*self.xx) + self.Oi0*a**(6*self.xx-6.)/(1- 2*self.xx)
        return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+NuContrib+rhow)

