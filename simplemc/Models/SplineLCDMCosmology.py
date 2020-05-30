##
# This si a cosmology where w(z) is defined by splines.
##


from scipy.interpolate import InterpolatedUnivariateSpline
from LCDMCosmology import LCDMCosmology
from ParamDefs import Sp1_par, Sp2_par, Sp3_par, Sp4_par
import math as N


class SplineLCDMCosmology(LCDMCosmology):
    def __init__(self, varySp1=True, varySp2=True, varySp3=True, varySp4=True):

        self.varySp1 = varySp1
        self.varySp2 = varySp2
        self.varySp3 = varySp3
        self.varySp4 = varySp4

        self.Sp1 = Sp1_par.value
        self.Sp2 = Sp2_par.value
        self.Sp3 = Sp3_par.value
        self.Sp4 = Sp4_par.value

        # Nodes are equally-spaced in log10(z)
        self.zmin   = 0.1                 # first-node position
        self.zmax   = 2.5                 # last-node position
        self.Nnodes = 6                   # number of nodes used

        self.lzmin = N.log10(self.zmin)
        self.lzmax = N.log10(self.zmax)

        LCDMCosmology.__init__(self)



    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varySp1): l.append(Sp1_par)
        if (self.varySp2): l.append(Sp2_par)
        if (self.varySp3): l.append(Sp3_par)
        if (self.varySp4): l.append(Sp4_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "Sp1":
                self.Sp1 = p.value
            elif p.name == "Sp2":
                self.Sp2 = p.value
            elif p.name == "Sp3":
                self.Sp3 = p.value
            elif p.name == "Sp4":
                self.Sp4 = p.value
        return True


    def Spline(self, a):
        z = 1.0/a-1.0
        x = [10**(self.lzmin+(self.lzmax-self.lzmin)/(self.Nnodes-1)*i)
             for i in range(0, int(self.Nnodes))]

        y = [1.0, self.Sp1, self.Sp2, self.Sp3, self.Sp4, 1.0]
        s = InterpolatedUnivariateSpline(x, y)
        ys = s(z)
        return ys


    def Rho_de(self, a):
        z = 1.0/a-1.0
        if (z > self.zmax or z < self.zmin):
            return 1.0
        else:
            return self.Spline(a)


    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        NuContrib = self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om)*self.Rho_de(a))
