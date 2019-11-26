# This is LCDM cosmology
# I didn't intend it that way, but it is used as a base class
# for most other cosmologies, mostly because it treats Neutrinos and Radiation
# hassle.

import sys
import scipy as sp
import CosmoApprox as CA
from BaseCosmology import BaseCosmology
from RadiationAndNeutrinos import RadiationAndNeutrinos
from ParamDefs import Obh2_par, Om_par, h_par, mnu_par, Nnu_par



class LCDMCosmology(BaseCosmology, RadiationAndNeutrinos):
    # possible options: "Anderson", "Cuesta", "CuestaNeff", "EH"
    rd_approx = "Cuesta"

    def __init__(self, Obh2=Obh2_par.value, Om=Om_par.value, h=h_par.value, mnu=mnu_par.value,
                 Nnu=Nnu_par.value, degenerate_nu=False, disable_radiation=False, fixOm=False):
        # two parameters: Om and h
        self.Om    = Om
        self.Obh2  = Obh2
        self.fixOm = fixOm

        BaseCosmology.__init__(self, h)
        RadiationAndNeutrinos.__init__(
            self, mnu, Nnu, degenerate_nu, disable=disable_radiation)
        if (type(self.rd_approx) == type(147.)):
            self.rd_func_ = lambda t1, t2, t3, t4: self.rd_approx
        elif (self.rd_approx == "Anderson"):
            self.rd_func_ = CA.rd_anderson_approx
        elif (self.rd_approx == "Cuesta"):
            self.rd_func_ = CA.rd_cuesta_approx
        elif (self.rd_approx == "CuestaNeff"):
            self.rd_func_ = CA.rd_cuesta_Nnu_approx
        elif (self.rd_approx == "EH"):
            self.rd_func_ = self.rd_EH_approx
        else:
            print("Bad rd Approx specified")
            sys.exit(1)

        # Force rd update
        LCDMCosmology.updateParams(self, [])
        RadiationAndNeutrinos.updateParams(self, [])
        # by default we keep Obh2 prior
        self.noObh2prior = False


    def setNoObh2prior(self, val=True):
        print("Disabling obh2 prior.")
        self.noObh2prior = val


    # my free parameters, note that h is defined in Base
    # to change parameters/priors see ParamDefs.py
    def freeParameters(self):
        Om_par.setValue(self.Om)
        l = []
        if not self.fixOm:
            l.append(Om_par)
        if (not self.varyPrefactor):
            Obh2_par.setValue(self.Obh2)
            l += [Obh2_par]

        l += BaseCosmology.freeParameters(self)
        l += RadiationAndNeutrinos.freeParameters(self)
        return l


    def updateParams(self, pars):
        BaseCosmology.updateParams(self, pars)
        RadiationAndNeutrinos.updateParams(self, pars)

        for p in pars:
            if p.name == "Om":
                self.Om = p.value
            elif p.name == "Obh2":
                self.Obh2 = p.value
        self.Ocb = self.Om-self.Omnu-self.Omrad
        if (not self.varyPrefactor):
            Nnu = self.Nnu()
            # if we have disabled radiation, let use just use std value
            # to prevent rd calculator from crashing. Not ever realistically
            # used
            if (Nnu == 0):
                Nnu = 3.04
            self.setrd(self.rd_func_(
                self.Obh2, self.Ocb*self.h**2, self.Omnuh2, Nnu))
        return True


    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    # @autojit
    def RHSquared_a(self, a):
        NuContrib = self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om))


    # Obh2 prior
    def prior_loglike(self):
        if (self.varyPrefactor or self.noObh2prior):
            return 0
        ## put back in if needed
        # Cooke et al, http://arxiv.org/abs/1308.3240
        # 2.202 +/- 0.046
        #return -(self.Obh2-0.02202)**2/(2*0.00046**2)
        return 0


    # this returns the Wang+Wang variables in a vec
    def WangWangVec(self):
        Omh2   = self.Ocb*self.h**2+self.Omnuh2
        omt    = Omh2/self.h**2
        zstar  = self.CA.z_lastscattering(Omh2, self.Obh2)
        Dastar = self.Da_z(zstar)*self.c_/(self.h*100)
        R      = sp.sqrt(omt)*self.h*100*Dastar/self.c_
        la     = sp.pi*Dastar/self.CA.soundhorizon_star(Omh2, self.Obh2)
        # print la, R, self.Obh2
        return sp.array([la, R, self.Obh2])


    # this returns the "SimpleCMB" variables in a vec
    def CMBSimpleVec(self):
        Ocbh2  = self.Ocb*self.h**2
        zstar  = 1090
        Dastar = self.Da_z(zstar)*self.c_/(self.h*100)
        return sp.array([self.Obh2, Ocbh2, Dastar/self.rd])
