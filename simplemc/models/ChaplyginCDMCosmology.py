

from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import As_par, Calpha_par, Ok_par





class ChaplyginCDMCosmology(LCDMCosmology):
    """
        This is CDM cosmology with Chaplygin gas.
        This class inherits LCDMCosmology class as the rest of the cosmological
        models already included in SimpleMC.

        :param varyas: variable As parameter
        :param varyalpha: variable beta parameter


    """
    def __init__(self, varyas=True, varyalpha=True, varyOk=False, fixOm=True):
        self.varyas  = varyas
        self.varyalpha = varyalpha
        self.varyOk = varyOk

        self.Ok = Ok_par.value
        self.As = As_par.value
        self.alpha = Calpha_par.value
        LCDMCosmology.__init__(self, fixOm=fixOm)


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varyas):  l.append(As_par)
        if (self.varyalpha): l.append(Calpha_par)
        if (self.varyOk): l.append(Ok_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "as":
                self.As = p.value
            elif p.name == "alpha":
                self.alpha = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False
        return True


    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        NuContrib = 0 #self.NuDensity.rho(a)/self.h**2
        rhow = (self.As + (1-self.As)*a**(-3*(1+self.alpha)))
        return (self.Obh2/self.h**2/a**3 + self.Ok/a**2+self.Omrad/a**4+NuContrib+(1.0-self.Obh2/self.h**2-self.Ok)*rhow)
