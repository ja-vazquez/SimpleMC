


import numpy as np
from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import Ok_par, ggama_par, glambda_par

#Graduated dark energy: Observational hints of a spontaneous sign
# switch in the cosmological constant
#https://arxiv.org/abs/1912.08751

class GraduatedCosmology(LCDMCosmology):
    def __init__(self, varyggama=True, varyglambda=False, varyOk=False):
        """
        Introduces the graduated dark energy (gDE) characterised by a minimal dynamical
        deviation from the null inertial mass density of the Lambda
        Parameters
        ----------
        varyggama
        varyglambda
        varyOk

        Returns
        -------

        """

        self.varyggama   = varyggama
        self.varyglambda = varyglambda
        self.varyOk      = varyOk

        self.Ok      = Ok_par.value
        self.ggama   = ggama_par.value
        self.glambda = glambda_par.value
        LCDMCosmology.__init__(self)


    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyggama): l.append(ggama_par)
        if (self.varyglambda): l.append(glambda_par)
        if (self.varyOk): l.append(Ok_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="ggama":
                self.ggama=p.value
            elif p.name=="glambda":
                self.glambda=p.value
            elif p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False
        return True


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        #For glambda=0, it becomes inertial mass density
        z = 1./a - 1.
        #NuContrib= self.NuDensity.rho(a)/self.h**2
        term = (1. - 3.*self.ggama*(self.glambda-1)*np.log(1+z))
        if self.glambda == 1:
            rhow = 1.
        elif self.glambda ==0:
            rhow = (1. + 3.*self.ggama*np.log(1+z))
        else:
            rhow= np.sign(term)*np.abs(term)**(1./(1-self.glambda))
        return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+(1.0-self.Om-self.Ok)*rhow)

