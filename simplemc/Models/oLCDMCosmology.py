# This is LCDM cosmology with optional
# curvature which you can set up with
# setVaryOk()

from LCDMCosmology import LCDMCosmology
from ParamDefs import Ok_par

class oLCDMCosmology(LCDMCosmology):
    # zeroDE forces Ol to zero.
    def __init__(self, zeroDE=False, kwargs_LCDM={}):
        # two parameters: Om and h
        self.Ok     = Ok_par.value
        self.zeroDE = zeroDE
        LCDMCosmology.__init__(self, **kwargs_LCDM)


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        Ok_par.setValue(self.Ok)
        l = []
        if not self.zeroDE:
            l.append(Ok_par)
        return l+LCDMCosmology.freeParameters(self)


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "Ok":
                self.Ok = p.value
                # se comments to BaseCosmology class
                # in short, WE hold the Ok parameter
                # (maybe you want to parameterize it differently),
                # but the BaseCosmology class needs to know about Ok in order
                # to get Da and its sin/sinhs right.
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False
        if self.zeroDE:
            # set curvature so that Ode=0
            self.Ok = 1-self.Om
            self.setCurvature(self.Ok)
        return True


    # this is relative Hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        NuContrib = self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+NuContrib+(1.0-self.Om-self.Ok))
