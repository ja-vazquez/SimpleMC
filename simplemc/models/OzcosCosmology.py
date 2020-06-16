

from simplemc.cosmo.paramDefs import *
from simplemc.models.LCDMCosmology import LCDMCosmology
import math as ma

#TODO still in progress

class OzcosCosmology(LCDMCosmology):
    def __init__(self, varynk=True, varykk=True, varyw=False, varywa=False, varywb=False, varyOk=False, fixOm=False):


        self.varyw=varyw
        self.varywa=varywa
        self.varyOk=varyOk
        self.varywb=varywb
        self.varynk=varynk
        self.varykk=varykk

        self.nk=nk_par.value
        self.kk=kk_par.value

        self.Ok=Ok_par.value
        self.w0=w_par.value
        self.wa=wa_par.value
        self.wb=wb_par.value

        LCDMCosmology.__init__(self,fixOm=fixOm)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyw): l.append(w_par)
        if (self.varywa): l.append(wa_par)
        if (self.varywb): l.append(wb_par)
        if (self.varyOk): l.append(Ok_par)

        if (self.varynk): l.append(nk_par)
        if (self.varykk): l.append(kk_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="w":
                self.w0=p.value
            elif p.name=="wa":
                self.wa=p.value
            elif p.name=="wb":
                self.wb=p.value
            elif p.name=="nk":
                self.nk=p.value
            elif p.name=="kk":
                self.kk=p.value
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

       z = 1./a - 1.
       if self.varynk:
            n = self.nk
            k = self.kk

            fac1 = (ma.sin(k/2.0))**4
            fac2 = 3.0*ma.cos(k/2.0)/ma.sin(k/2.0) + 2.0*n

            Om0 = (1.0/81)*fac1*fac2**4
            Ox0 = -(8.0*n/81)*fac1*fac2**3
            #Ol0 = (2.0/27)*fac1*fac2**2*(3.0 + 4.0*n**2)
            Oy0 = -(8.0*n/81)*fac1*fac2*(9.0 + 4.0*n**2)
            Oz0 = (1.0/81)*fac1*(9.0 + 4*n**2)**2
            Ol0 = 1 - self.Omrad - Om0 - Ox0 - Oy0 - Oz0

            new_term = Om0*(1+z)**3  + Ox0*(1+z)**(1.5) + Oy0*(1+z)**(-1.5) + Oz0*(1+z)**(-3) + Ol0
       else:      
            Ol0 = 1 - self.Omrad - self.Ocb - (1-ma.sqrt(self.Ocb))**2
            new_term = self.Ocb*(1+z)**3 + (1-ma.sqrt(self.Ocb))**2*(1+z)**(3*(1+self.w0)) + Ol0

       return (self.Omrad/a**4 +NuContrib  +new_term )



    def w_de(self, z):
        a = 1./(1+z)
        NuContrib=self.NuDensity.rho(a)/self.h**2
        n = self.nk
        k = self.kk
   
        fac1 = (ma.sin(k/2.0))**4
        fac2 = 3.0*ma.cos(k/2.0)/ma.sin(k/2.0) + 2.0*n
        self.Om0 = (1.0/81)*fac1*fac2**4
        self.Ox0 = -(8.0*n/81)*fac1*fac2**3
        self.Oy0 = -(8.0*n/81)*fac1*fac2*(9.0 + 4.0*n**2)
        self.Oz0 = (1.0/81)*fac1*(9.0 + 4*n**2)**2
        self.Ol0 = 1 - self.Omrad - self.Om0 - self.Ox0 - self.Oy0 - self.Oz0

        top = 0.5*self.Ox0*(1+z)**1.5 - 0.5*self.Oy0*(1+z)**(-1.5) - self.Oz0*(1+z)**(-3)
        bot = self.Ox0*(1+z)**1.5 + self.Ol0 + self.Oy0*(1+z)**(-1.5) + self.Oz0*(1+z)**(-3)

        return -1 + top/bot
