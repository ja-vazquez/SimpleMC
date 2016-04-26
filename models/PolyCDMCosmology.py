## This is LCDM cosmology with optional
## curvature which you can set up with 
## setVaryOk()

from LCDMCosmology import *

class PolyCDMCosmology(LCDMCosmology):
    def __init__(self, varyOk=True):
        ## Ok, LCDM has Omega_m, we also have Omega_1 and Omega_2
        ## and Lambda is then what remains
        ##
        self.Ok =Ok_par.value
        self.Om1=Om1_par.value
        self.Om2=Om2_par.value
        self.varyOk=varyOk
        LCDMCosmology.__init__(self,mnu=0)

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)+[Om1_par,Om2_par]
        if self.varyOk:
            l.append(Ok_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="Om1":
                self.Om1=p.value
            elif p.name=="Om2":
                self.Om2=p.value
            elif p.name=="Ok":
                self.Ok=p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok)>1.0):
                   return False
        return True

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        ## ok, there is really no point in adding neutrinos to this
        return (self.Om/a**3+self.Om1/a**2+self.Ok/a**2+self.Om2/a+(1-self.Om-self.Om1-self.Om2-self.Ok))

    def prior_loglike(self):
        return (-self.Ok**2/(2*0.1**2) ## A 0.1 prior in Ok as discussed at the telecon
                +LCDMCosmology.prior_loglike(self))
