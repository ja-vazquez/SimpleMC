
import math as N
import numpy as np
from simplemc.models.LCDMCosmology import LCDMCosmology
from scipy.integrate import quad
from scipy.interpolate import interp1d
from simplemc.cosmo.Parameter import Parameter

import matplotlib.pyplot as plt

class TanhCosmology(LCDMCosmology):
    def __init__(self):


        self.Nbins_eos = 3
        mean_eos = -1
        self.params = [Parameter("zbin_eos%d"%i, mean_eos, 0.2, (-3.5, 0), "zbin_eos%d"%i) for i in range(self.Nbins_eos)]
        self.pvals = [i.value for i in self.params]
        self.z_i = np.linspace(0.0, 3.0, self.Nbins_eos+1)

        self.zinter = np.linspace(0.0, 3.0, 50)

        LCDMCosmology.__init__(self, mnu=0)
        self.updateParams([])


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        #l.append(zbin_rho_par)
        l+= self.params
        return l

    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            for i in range(self.Nbins_eos):
                if p.name == ("zbin_eos"+str(i)):
                    self.pvals[i] = p.value

        self.initialize()
        return True

     
    def bines(self, w_2, w_1, z_2, z_1, eta):
        return (w_2-w_1)*(1.0+np.tanh((z_2-z_1)/eta))/2.0


    def de_eos(self, z):
        w = self.pvals[0]
        for jj in range(self.Nbins_eos - 1):
            w+=self.bines(self.pvals[jj+1], self.pvals[jj], z, self.z_i[jj+1], eta=0.15)
        return w


    def de_rhow(self, z):
        eos = self.de_eos(z)
        resultado = quad(lambda b: 3.0*(1.0+ eos)/(1.0+b), 0.0, z )
        return resultado[0]
   

    def initialize(self):
        #w_inter = [self.de_eos(z) for z in self.zinter]
        rhow = [self.de_rhow(z) for z in self.zinter]
        self.rhow_inter = interp1d(self.zinter, rhow)

        #plt.plot(self.zinter, rhow)
        #plt.show()
        return True

    


    def luisfunction(self, z):
        def luisfunction2(z):
            #w_i = self.pvals
            #z_i = np.linspace(0.0,3.0, self.Nbins_eos+1)
            #def bines(w_2,w_1,z_2,z_1,eta):
            #    return (w_2-w_1)*(1.0+np.tanh((z_2-z_1)/eta))/2.0
            w = self.pvals[0]
            for jj in range(self.Nbins_eos - 1):
                w+=self.bines(self.pvals[jj+1],self.pvals[jj], z, self.z_i[jj+1],eta=0.15)
            rhow=w
            return rhow
        rhow=luisfunction2(z)
        resultado = quad(lambda b: 3.0*(1.0+rhow)/(1.0+b), 0.0, z )
        return luisfunction2(z), resultado[0]


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        z= 1./a - 1
        if z>= 3.0:
            rhow = (1.0-self.Om)
        else:
            rhow = (1.0-self.Om)*np.exp(self.rhow_inter(z)) 
        #return (self.Ocb/a**3+self.Omrad/a**4 +(1.0-self.Om)*(np.exp(self.luisfunction(z)[1]))  )
        return self.Ocb/a**3 + self.Omrad/a**4 + rhow

