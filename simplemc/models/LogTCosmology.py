


from simplemc.models.LCDMCosmology import LCDMCosmology
from  simplemc.cosmo.paramDefs import Ok_par, alpha_par
import numpy as np

#Screening Lambda in a new modified gravity model
#https://arxiv.org/abs/1903.11519

class LogTCosmology(LCDMCosmology):
    def __init__(self, varyalpha=True, varyOk=False):
        """
        We study a new model of Energy-Momentum Squared Gravity
        (EMSG), called Energy-Momentum Log Gravity (EMLG)
        Parameters
        ----------
        varyalpha
        varyOk

        Returns
        -------

        """

        self.varyalpha= varyalpha
        self.varyOk   = varyOk

        self.Ok    = Ok_par.value
        self.alpha = alpha_par.value
        LCDMCosmology.__init__(self)


    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyalpha): l.append(alpha_par)
        if (self.varyOk): l.append(Ok_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="alpha":
                self.alpha=p.value
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
        z   = 1./a - 1.
        alp = self.alpha
        term1   = (1.0 + alp)**2*(1+z)**3 - 2.*alp
        sqterm1 = np.sqrt(-4*alp**2 + term1**2)

        return self.Omrad/a**4 + 1 - self.Om*(1 - 0.5*(term1 + sqterm1) + \
               alp*np.log(0.5*(term1 + sqterm1)))

#        return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+(1.0-self.Om-self.Ok)*rhow)

