

from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import bd_par, Osig_par

#TODO still need to be tested

#Anisotropic massive Brans-Dicke gravity extension of the standard LCDM model
# https://arxiv.org/abs/1903.06679

class AnisotropicCosmology(LCDMCosmology):
    def __init__(self, varybd=True, varyOsig=True, bd_model=False):
        """
         Anisotropic massive Brans-Dicke (BD) gravity extension of the standard
         LCDM model, wherein the extension is characterized by two additional degrees
         of freedom; the BD parameter, w, and the present day density parameter
         corresponding to the shear scalar, Omega_sigma,0
        Parameters
        ----------
        varybd
        varyOsig

        Returns
        -------

        """

        self.varybd   = varybd
        self.varyOsig = varyOsig

        self.bd    = bd_par.value
        self.Osig  = Osig_par.value

        self.bd_model = bd_model
        LCDMCosmology.__init__(self)



    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varybd):   l.append(bd_par)
        if (self.varyOsig): l.append(Osig_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="bd":
                self.bd=p.value
            elif p.name=="Osig":
                self.Osig=p.value
            elif p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False
        return True


    def Ode_fz(self, z):
        sigma  = 10**self.Osig
        bdp    = 10**self.bd

        Ode    = (1.0 - self.Om - sigma)
        gamma  = 6.*(1 + bdp)**2/((3.*bdp + 4.)*(2.*bdp + 3))

        Ode_fz = Ode + sigma*(gamma - 1.)*((1+z)**(6 + 2./(1 + bdp)) - 1.) \
                + self.Ocb*((1 + z)**3*(gamma*(1 + z)**(1./(1. + bdp))-1.) + 1 - gamma)
        return Ode_fz



    def Osigma(self, z):
        bdp    = 10**self.bd
        sigma  = 10**self.Osig

        Osigma = sigma*(1 + z)**(6+2./(1. + bdp))
        return Osigma



    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        z     = 1./a - 1.0
        NuContrib = self.NuDensity.rho(a)/self.h**2


        if self.bd_model:
            H = self.Omrad/a**4 +NuContrib + self.Ocb/a**3  + self.Ode_fz(z) + self.Osigma(z)
        else:
            sigma = 10**self.Osig #4.0E-21
            H = self.Omrad/a**4 +NuContrib + self.Ocb/a**3 + sigma/a**6 + (1-self.Omrad-self.Ocb-sigma)

        return H


