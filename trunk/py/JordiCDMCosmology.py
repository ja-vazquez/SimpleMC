#
# This is a cosmological model that Jordi Miralda came up with.
# The world is full of mystery.
#


import math as N
from LCDMCosmology import *

class JordiCDMCosmology(LCDMCosmology):
    def __init__(self):
        ## two parameters: Om and h
        self.Ok=Ok_par.value   
        self.q=q_par.value
        self.za=za_par.value
        self.zb=zb_par.value
        self.wd=wd_par.value 
        self.Od=Od_par.value
        LCDMCosmology.__init__(self)
        

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
         return [Ok_par, q_par,za_par, zb_par,wd_par,Od_par]+LCDMCosmology.freeParameters(self)

    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="q":
                 self.q=p.value
            if p.name=="za":
                 self.za=p.value
            if p.name=="zb":
                 self.zb=p.value
            if p.name=="wd":
                 self.wd=p.value
            if p.name=="Od":
                 self.Od=p.value 
            elif p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
        return True


    def Az(self,a):
        z=1.0/a-1.0
        if((self.zb<z)&(z<self.za)):
          return 1-(self.q*(self.za-z))/(self.za-self.zb)
        elif(z>self.za):
          return 1.0
        elif(z<self.zb):
          return 1-self.q
        else:
          stop("BAD MODEL")
    
    def Qdrz(self,a):
        z=1.0/a-1.0
        if(z>self.zb):
          zmax=z
        else:
          zmax=self.zb        
        if(z<self.za):
          return self.q*(1+z)*N.log((1.0+self.za)/(1.0+zmax))/(self.za-self.zb) 
        else:
          return 0.0

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        Omnuh2=self.NuDensity.omnuh2today
        Omnu=Omnuh2/self.h**2
        NuContrib=self.NuDensity.rho(a)/self.h**2
        Olambda= (1.0-self.Az(1.0)*self.Om-self.Ok-self.Omrad-self.Omnu-self.Qdrz(1.0)-self.Od)
        Omegad = self.Od*a**(-3*(1.0+self.wd))
        return (self.Omrad/a**4+self.Az(a)*self.Om/a**3+self.Ok/a**2+NuContrib+Olambda+self.Qdrz(a)+Omegad)

    def prior_loglike(self):
        if (self.za<=self.zb):
            return -100
        else:
            return LCDMCosmology.prior_loglike(self)

