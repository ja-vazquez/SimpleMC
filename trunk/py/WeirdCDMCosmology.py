## This is what is called Oscillating model or whatever
## It is really a scam to decrease lya BAO chi2.


import math as N

from LCDMCosmology import *

class WeirdCDMCosmology(LCDMCosmology):
    def __init__(self, varymu=True, varyAmp=True, varysig=True, varyCos=True):
        ## two parameters: Om and h

        ## we start with false here...
        varyOk=False
        ## this is my "original cosmology" -- outside gaussian not much will change.

        self.varyOk=varyOk
        self.varymu=varymu
        self.varyAmp=varyAmp
        self.varysig=varysig
        self.varyCos=varyCos

        self.Ok=Ok_par.value   
        self.mu=mu_par.value
        self.Amp=Amp_par.value
        self.sig=sig_par.value

        ## auxiliary self cosmology
        self.oC=LCDMCosmology()

        LCDMCosmology.__init__(self)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=[]
        if (self.varyCos): l+=LCDMCosmology.freeParameters(self)
        if (self.varyOk): l.append(Ok_par)
        if (self.varymu): l.append(mu_par)
        if (self.varyAmp): l.append(Amp_par)
        if (self.varysig): l.append(sig_par)
        return l


    def updateParams(self,pars):
        LCDMCosmology.updateParams(self,pars)
        self.oC.updateParams(pars)

        for p in pars:
            if p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False
            elif p.name=="mu":
                self.mu=p.value
            elif p.name=="Amp":
                self.Amp=p.value
            elif p.name=="sig":
                self.sig=p.value
        return True


    def Weird(self,a):
        #return W(z) and dW(z)/dz
        z=1.0/a-1
        lmu=log(self.mu+1)
        lz=log(z+1)
        W=self.Amp*N.exp(-(lz-lmu)**2/(2*self.sig**2))
        WP=W*(-(lz-lmu)/self.sig**2)*a ## a=1/(1+z) = dlog(1+z)/dz
        return W,WP


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        W,WP=self.Weird(a)
        H2=LCDMCosmology.RHSquared_a(self,a)
        HI=1/sqrt(H2)
        HI=HI*(1+W) + self.oC.Da_z(1/a-1)*WP
        return 1/HI**2

    def Da_z(self,z):
       W,WP=self.Weird(1/(1.0+z))
       return self.oC.Da_z(z)*(1+W)
