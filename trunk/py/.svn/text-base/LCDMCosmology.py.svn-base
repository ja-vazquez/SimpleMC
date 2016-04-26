## This is LCDM cosmology 
## I didn't indend it that way, but it is used as a base class
## for most other cosmologies, mostly because it treats Neutrinos and Radiation
## hassle.

from BaseCosmology import *
import CosmoApprox as CA

class LCDMCosmology(BaseCosmology,RadiationAndNeutrinos):

    ## possible options: "Anderson", "Cuesta", "CuestaNeff", "EH"
    rd_approx="Cuesta"

    def __init__(self,Obh2=Obh2_par.value, Om=Om_par.value, h=h_par.value, mnu=mnu_par.value, 
                 Nnu=Nnu_par.value, degenerate_nu=False, disable_radiation=False, fixOm=False):
        ## two parameters: Om and h
        self.Om=Om
        self.Obh2=Obh2
        self.fixOm=fixOm
        BaseCosmology.__init__(self,h)
        RadiationAndNeutrinos.__init__(self, mnu, Nnu, degenerate_nu, disable=disable_radiation)
        if (self.rd_approx=="Anderson"):
            self.rd_func_=CA.rd_anderson_approx
        elif (self.rd_approx=="Cuesta"):
            self.rd_func_=CA.rd_cuesta_approx
        elif (self.rd_approx=="CuestaNeff"):
            self.rd_func_=CA.rd_cuesta_Nnu_approx
        elif (self.rd_approx=="EH"):
            self.rd_func_=self.rd_EH_approx
        else:
            Error("Bad rd Approx specified")

        # Force rd update
        LCDMCosmology.updateParams(self,[])
        RadiationAndNeutrinos.updateParams(self,[])
        ## by default we keep Obh2 prior
        self.noObh2prior=False

    def setNoObh2prior(self,val=True):
        print "Disabling obh2 prior."
        self.noObh2prior=val

    # my free parameters, note that h is defined in Base
    # to change parameters/priors see ParamDefs.py    
    def freeParameters(self):
        Om_par.setValue(self.Om)
        l=[]
        if not self.fixOm:
            l.append(Om_par)
        if (not self.varyPrefactor):
            Obh2_par.setValue(self.Obh2)
            l+=[Obh2_par]
            
        l+=BaseCosmology.freeParameters(self)
        l+=RadiationAndNeutrinos.freeParameters(self)
        return l

    def updateParams(self,pars):
        BaseCosmology.updateParams(self,pars)
        RadiationAndNeutrinos.updateParams(self,pars)

        for p in pars:
            if p.name=="Om":
                self.Om=p.value
            elif p.name=="Obh2":
                self.Obh2=p.value
        self.Ocb=self.Om-self.Omnu-self.Omrad
        if (not self.varyPrefactor):
            Nnu=self.Nnu()
            # if we have disabled radiation, let use just use std value
            # to preven rd calculator from crashing. Not ever realistically
            # used
            if (Nnu==0):
                Nnu=3.04
            self.setrd(self.rd_func_(self.Obh2, self.Ocb*self.h**2,self.Omnuh2, Nnu))
        return True

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    #@autojit
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om))

    ## Obh2 prior
    def prior_loglike(self):
        ## hack
        #och2=self.Ocb*self.h**2-self.Obh2
        #if (self.och2<0.1103000) or (self.och2>0.1289000):
        #    if (rd_approx=="tabulated_Nnu"):
        #        return -1e50

        if (self.varyPrefactor or self.noObh2prior):
            return 0 
        else:
            ## Cooke et al, http://arxiv.org/abs/1308.3240
            ## 2.202 +/- 0.046
            return -(self.Obh2-0.02202)**2/(2*0.00046**2)

    ## this returns the Wang+Wang variables in a vec
    def WangWangVec(self):
        Omh2=self.Ocb*self.h**2+self.Omnuh2
        omt=Omh2/self.h**2
        zstar=self.CA.z_lastscattering(Omh2,self.Obh2)
        Dastar=self.Da_z(zstar)*self.c_/(self.h*100)
        R=sqrt(omt)*self.h*100*Dastar/self.c_
        la=pi*Dastar/self.CA.soundhorizon_star(Omh2, self.Obh2)
        #print la, R, self.Obh2
        #stop()
        return array([la,R,self.Obh2])

    ## this returns the "SimpleCMB" variables in a vec
    def CMBSimpleVec(self):
        Ocbh2=self.Ocb*self.h**2
        zstar=1090
        Dastar=self.Da_z(zstar)*self.c_/(self.h*100)
        return array([self.Obh2, Ocbh2, Dastar/self.rd])
            

        
                
