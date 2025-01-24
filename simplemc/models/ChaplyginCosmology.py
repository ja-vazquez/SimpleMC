

from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
from simplemc.cosmo.paramDefs import Ok_par
from scipy.interpolate import interp1d
from scipy import integrate
import numpy as np





class ChaplyginCosmology(LCDMCosmology):
    """
        This is CDM cosmology with Chaplygin gas.
        This class inherits LCDMCosmology class as the rest of the cosmological
        models already included in SimpleMC.

        :param varyas: variable As parameter
        :param varyalpha: variable beta parameter


    """
    def __init__(self, varyas=True, varyalpha=True, varybeta=True,
                 usesigmoid=False, varyOk=False, fixOm=True):
        # Chaplygin cosmology
        self.as_par = Parameter("as", 1.0, 0.1, (0.8, 1.2), "A_{s}")
        self.alpha_par = Parameter("alpha", 5, 0.1, (2.0, 20.0), "\\alpha")
        self.beta_par = Parameter("beta", 1., 0.1, (-0.2, 2.5), "\\beta")

        self.varyOk = varyOk
        self.varyas  = varyas
        self.varyalpha = varyalpha
        self.varybeta = varybeta

        self.Ok = Ok_par.value
        self.as_chap = self.as_par.value
        self.alpha_chap = self.alpha_par.value
        self.beta_chap = self.beta_par.value

        self.use_sigmoid = usesigmoid
        self.zvals = np.linspace(0, 3, 100)

        LCDMCosmology.__init__(self, fixOm=fixOm)
        if self.use_sigmoid:
            self.updateParams([])


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varyOk): l.append(Ok_par)
        if (self.varyas):  l.append(self.as_par)
        if (self.varyalpha): l.append(self.alpha_par)
        if (self.varybeta): l.append(self.beta_par)
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
        if self.use_sigmoid:
            self.initialize()
        return True


    def sigmoid(self, z):
        return 1./(1 + np.exp( -self.alpha_chap*(z-self.beta_chap) ))

    def integrand(self, z):
        return ( 1+self.as_chap*(self.sigmoid(z)-1) )/(1+z)

    def initialize(self):
        #in this case, we consider the pressure is given by
        # P=As*(\sigm(z)-1)\rho
        integ= lambda z: integrate.quad(self.integrand, 0, z)[0]
        self.rhow = interp1d(self.zvals, list(map(integ, self.zvals)))

    def auxiliar(self, z):
        return self.as_chap + (1-self.as_chap)*(1+z)**( 3*(1+self.beta_chap)*(1+self.alpha_chap))

    def rho_chap(self, z):
        rhow=  self.auxiliar(z)**(1/(self.alpha_chap))
        return rhow

    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        z= 1./a-1
        NuContrib = 0 #self.NuDensity.rho(a)/self.h**2
        standard = self.Obh2/self.h**2/a**3 + self.Ok/a**2+self.Omrad/a**4+NuContrib
        if self.use_sigmoid:
            rhow= self.rhow(z)
        else:
            rhow= self.rho_chap(z)
        return standard + (1.0-self.Obh2/self.h**2-self.Ok)*rhow


    def eos(self, z):
        if self.use_sigmoid:
            w = self.as_chap*( self.sigmoid(z)-1 )
        else:
            w = self.beta_chap - self.alpha_chap/self.auxiliar(z)
        return w