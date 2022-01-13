
import math as N
import numpy as np
from LCDMCosmology import *
from scipy.integrate import quad
from scipy.interpolate import interp1d
from ParamDefs import *
class eos_tanh(LCDMCosmology):
    def __init__(self):

        self.parvals = zbin_eos_par
        self.zs = [i.value for i in self.parvals]
        names   = [i.name  for i in self.parvals]
        self.index = dict((j,i) for i,j in enumerate(names))

        #self.zvals = np.logspace(np.log10(0.01),np.log10(2.261), len(self.parvals)+1)

        LCDMCosmology.__init__(self)

    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        #l.append(zbin_rho_par)
        l+= self.parvals
        return l

    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            for i in range(len(self.parvals)):
                if p.name == ("zbin_eos"+str(i)):
                    self.zs[i] = p.value
        return True

    def luisfunction(self, z):
        def luisfunction2(z):
            w_i = self.zs
            z_i = np.linspace(0.0,3.0,len(self.parvals)+1)
            def bines(w_2,w_1,z_2,z_1,eta):
                return (w_2-w_1)*(1.0+np.tanh((z_2-z_1)/eta))/2.0
            w = self.zs[0]
            for jj in range(len(self.zs)-1):
                w+=bines(self.zs[jj+1],self.zs[jj],z,z_i[jj+1],eta=0.15)
            rhow=w
            return rhow
        rhow=luisfunction2(z)
        resultado = quad(lambda b: 3.0*(1.0+rhow)/(1.0+b), 0.0, z )
        return luisfunction2(z), resultado[0]


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        z= 1./a - 1
        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om)*(np.exp(self.luisfunction(z)[1]))  )


