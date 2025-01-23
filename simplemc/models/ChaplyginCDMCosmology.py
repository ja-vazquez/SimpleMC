

from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
from simplemc.cosmo.paramDefs import Ok_par





class ChaplyginCDMCosmology(LCDMCosmology):
    """
        This is CDM cosmology with Chaplygin gas.
        This class inherits LCDMCosmology class as the rest of the cosmological
        models already included in SimpleMC.

        :param varyas: variable As parameter
        :param varyalpha: variable beta parameter


    """
    def __init__(self, varyas=True, varyalpha=True, varybeta=True, varyOk=False, fixOm=True):
        # Chaplygin cosmology
        as_par = Parameter("as", 0.5, 0.1, (0.0, 1.0), "A_{s}")
        alpha_par = Parameter("alpha", 0.3, 0.1, (0.0, 1.0), "\\alpha")
        beta_par = Parameter("beta", 0.3, 0.1, (0.0, 1.0), "\\beta")

        self.varyOk = varyOk
        self.varyas  = varyas
        self.varyalpha = varyalpha
        self.varybeta = varybeta

        self.Ok = Ok_par.value
        self.as_chap = as_par.value
        self.alpha_chap = alpha_par.value
        self.beta_chap = beta_par.value
        LCDMCosmology.__init__(self, fixOm=fixOm)


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varyOk): l.append(Ok_par)
        if (self.varyas):  l.append(as_par)
        if (self.varyalpha): l.append(alpha_par)
        if (self.varybeta): l.append(beta_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "as":
                self.as_chap = p.value
            if p.name == "alpha":
                self.alpha_chap = p.value
            if p.name == 'beta':
                self.beta_chap = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False
        return True


    def auxiliar(self, z):
        return self.as_chap + (1-self.as_chap)*(1+z)**( -3*(1+self.beta_chap)*(1+self.alpha_chap))

    def rho_chap(self, z):
        rhow=  self.auxiliar(z)**(1/(self.alpha_chap))
        return rhow

    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        z= 1./a-1
        NuContrib = 0 #self.NuDensity.rho(a)/self.h**2
        standard = self.Obh2/self.h**2/a**3 + self.Ok/a**2+self.Omrad/a**4+NuContrib
        return standard + (1.0-self.Obh2/self.h**2-self.Ok)*self.rho_chap(z)


    def eos(self, z):
        self.beta_chap - self.alpha_chap/self.auxiliar(z)