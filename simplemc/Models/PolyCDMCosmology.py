# This is LCDM cosmology with optional
# free parameters on the Hubble function

from LCDMCosmology import LCDMCosmology
from ParamDefs import Ok_par, Om1_par, Om2_par


class PolyCDMCosmology(LCDMCosmology):
    def __init__(self, polyvary=['Om1','Om2','Ok'], Ok_prior=0.1):
        # Ok, LCDM has Omega_m, we also have Omega_1 and Omega_2
        # and Lambda is then what remains
        ##
        self.Ok  = Ok_par.value
        self.Om1 = Om1_par.value
        self.Om2 = Om2_par.value
        self.varyOm1  = 'Om1' in polyvary
        self.varyOm2  = 'Om2' in polyvary
        self.varyOk   = 'Ok' in polyvary
        self.Ok_prior = Ok_prior
        LCDMCosmology.__init__(self, mnu=0)


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if self.varyOm1: l.append(Om1_par)
        if self.varyOm2: l.append(Om2_par)
        Ok_par.setError(0.7)
        if self.varyOk:  l.append(Ok_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "Om1":
                self.Om1 = p.value
            elif p.name == "Om2":
                self.Om2 = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False
        return True


    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        # ok, there is really no point in adding neutrinos to this
        # as 3/23 -- note Om1 and Om2 used to be swapped
        return (self.Om/a**3+self.Om2/a**2+self.Ok/a**2+self.Om1/a+(1-self.Om-self.Om1-self.Om2-self.Ok))


    def prior_loglike(self):
        return (-self.Ok**2/(2*self.Ok_prior**2)  # A 0.1 prior in Ok as discussed at the telecon
                + LCDMCosmology.prior_loglike(self))
