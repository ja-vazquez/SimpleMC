#Fourier cosmology model 
## designed by D Tamayo 
# 12.Oct. 2018

import math as N
from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import *
from scipy import interpolate, integrate
import numpy as np
import scipy as sci
from scipy import special
import matplotlib.pyplot as plt

#Fourier-series expansion of the dark-energy equation of state
#https://arxiv.org/abs/1901.08679

class wFourierCosmology(LCDMCosmology):
    def __init__(self, varya0=False, varya1=True, varyb1=True, varya2=True, \
                       varyb2=True, varyOk=False):
        """
        A general description for the dark-energy equation-of-state w(z) in
        the form of Fourier series.
        Parameters
        ----------
        varya0
        varya1
        varyb1
        varya2
        varyb2
        varyOk

        Returns
        -------

        """

        self.varya0 = varya0
        self.varya1 = varya1
        self.varyb1 = varyb1
        self.varya2 = varya2
        self.varyb2 = varyb2
        self.varyOk = varyOk

        self.Ok = Ok_par.value
        self.a0 = a0_par.value
        self.a1 = a1_par.value
        self.b1 = b1_par.value
        self.a2 = a2_par.value
        self.b2 = b2_par.value

        LCDMCosmology.__init__(self)

        self.zini = 3.2
        self.zmed = 3.0
        self.zfin = 0.0

        self.aini = 1.0/(1.0 + self.zini)
        self.amed = 1.0/(1.0 + self.zmed)
        self.afin = 1.0/(1.0 + self.zfin)

        T         = self.afin - self.amed
        self.theta = 2.0*N.pi/T

        self.npoints = 30
        self.ascale  = np.linspace(self.aini, self.afin, self.npoints)

        self.updateParams([])


    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varya0): l.append(a0_par)
        if (self.varya1): l.append(a1_par)
        if (self.varyb1): l.append(b1_par)
        if (self.varya2): l.append(a2_par)
        if (self.varyb2): l.append(b2_par)
        if (self.varyOk): l.append(Ok_par)
        return l

    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="a0":
                self.a0=p.value
            elif p.name=="a1":
                self.a1=p.value
            elif p.name=="b1":
                self.b1=p.value
            elif p.name=="a2":
                self.a2=p.value
            elif p.name=="b2":
                self.b2=p.value
            elif p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False

        #initialize calculations
        self.calc_gs()
        #plt.plot(self.ascale,  [self.RHSquared_a(i) for i in self.ascale])
        #plt.show()
        return True



    def interp(self, x, y):
        return interpolate.interp1d(x, y, kind='cubic')


    def calc_gs(self):
        #compute g to the calculate rhos'
        a     = self.ascale
        lena  = np.ones(len(a))

        theta = self.theta
        amed  = self.amed
        aini  = self.aini
        amthe = amed*theta
        gamma = (1.0 + 0.5*self.a0 + self.b1 + self.b2)/(amed-aini)


        Si_1, Ci_1         = sci.special.sici(theta)
        Si_amed, Ci_amed   = sci.special.sici(amthe)
        Si_a, Ci_a         = sci.special.sici(a*theta)
        Si_2a, Ci_2a       = sci.special.sici(a*2*theta)
        Si_2, Ci_2         = sci.special.sici(2*theta)
        Si_2amed, Ci_2amed = sci.special.sici(2*amthe)

        f1 = self.b1*N.cos(amthe)     - self.a1*N.sin(amthe)
        f2 = self.a1*N.cos(amthe)     + self.b1*N.sin(amthe)
        f3 = self.b2*N.cos(2.0*amthe) - self.a2*N.sin(2.0*amthe)
        f4 = self.a2*N.cos(2.0*amthe) + self.b2*N.sin(2.0*amthe)

        ga  = (Ci_a*f1 + Si_a*f2 + Ci_2a*f3 + Si_2a*f4)*lena
        g1 = (Ci_1*f1   + Si_1*f2   + Ci_2*f3   + Si_2*f4)*lena
        g2 = (amed*gamma + f1*(Ci_1-Ci_amed) + f2*(Si_1-Si_amed) + f3*(Ci_2-Ci_2amed) + f4*(Si_2-Si_2amed))*lena
        g3 = (aini - amed)*gamma - f1*(Ci_1-Ci_amed) - f2*(Si_1-Si_amed) - f3*(Ci_2-Ci_2amed) - f4*(Si_2-Si_2amed)


        rho_C = aini**(3.0*aini*gamma)*amed**(-3.0*(1.0 + 0.5*self.a0 +aini*gamma))*np.exp(-3.0*g3)
        rho_L = a**(3.0*aini*gamma)*amed**(-3.0*(1.0 + 0.5*self.a0 +aini*gamma))*np.exp(-3.0*(a*gamma-g2))
        rho_F = a**(-3.0*(1.0 + 0.5*self.a0))*np.exp(-3.0*(ga-g1))

        self.rhoC       =  rho_C
        self.inter_rhoL =  self.interp(a, rho_L)
        self.inter_rhoF =  self.interp(a, rho_F)

        return True



    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):

        backg = self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4

        if a > self.amed:
            return  backg +(1.0-self.Om-self.Ok)*self.inter_rhoF(a)
        elif self.amed >= a >= self.aini:
            return  backg +(1.0-self.Om-self.Ok)*self.inter_rhoL(a)
        else:
            return  backg +(1.0-self.Om-self.Ok)*self.rhoC
 
