#
# This module calculates the predictions for the evolution
# of neutrino energy densities.
#


import scipy as sp
from scipy import constants as ct
from scipy.special import zeta
from scipy.integrate import quad
from scipy.interpolate import interp1d
import CosmoApprox as CA
#from numba import autojit


class NuIntegral:
    def __init__(self):
        "this integrates the stupid integral"
        print("Initalizing nu density look up table...", end=' ')
        rat = 10**(sp.arange(-4, 5, 0.1))
        intg = []
        for r in rat:
            # min below is to supress the stupid overlow warning
            res = quad(lambda x: sp.sqrt(x*x+r**2) /
                       (sp.exp(min(x, 400))+1.0)*x**2, 0, 1000)
            intg.append(res[0]/(1+r))
        intg = sp.array(intg)
        intg *= 7/8./intg[0]  # the right normalization

        self.interpolator = interp1d(sp.log(rat), intg)
        # type this into maple
        # evalf(45*Zeta(3)/(2*Pi^4));  0.2776566337
        self.int_infty = 45*zeta(3)/(2*ct.pi**4)
        print("Done")  # self.int_infty,intg[-1]*(1+r)/r


    def SevenEights(self, mnuOT):
        if (mnuOT < 1e-4):
            return 7/8.
        elif (mnuOT > 1e4):
            return self.int_infty*mnuOT  # I don't think this ever matters
        else:
            return self.interpolator(sp.log(mnuOT))*(1+mnuOT)


class ZeroNuDensity:
    # fake class that returns zeros if want to disable this
    def __init__(self):
        return

    def rho(self, a):
        return 0.0


class NuDensity:
    I = NuIntegral()

    def __init__(self, TCMB, Nnu=3.046, mnu=0.06, degenerate=False):
        # self.I=NuIntegral()
        # one neutrino species
        self.mnu_ = mnu
        self.Nnu_ = Nnu

        # this factor accounts for Neff=3.046 vs Neff=3
        # We make Tnu hotter by this factor and hence don't need to include it below.
        # It's all very academic, but what do we really mean by deltaNeff=1? Is it
        # one neutrino worth of radiation at the nominal temperature, or heated on?
        # See CAMB notes Eq. 4-7.

        self.gfact    = (3.046/3.0)
        self.gfact_o4 = self.gfact**(0.25)
        # ideal neutrino temp
        self.Tnu0     = (4./11.)**(1./3.)*TCMB
        # actual neutrino temp
        self.Tnu      = self.Tnu0*self.gfact_o4
        # same for prefactors
        self.prefix0  = 4.48130979e-7 * TCMB**4 * ((4./11.)**(4./3.))
        self.prefix   = self.prefix0*self.gfact

        self.degenerate = degenerate
        self.set_mnuone_()


    def set_mnuone_(self):
        if self.degenerate:
            self.mnuone  = self.mnu_/self.Nnu_
        else:
            self.mnuone  = self.mnu_*1.0
        self.omnuh2today = self.rho(1)


    def setMnu(self, mnu):
        self.mnu_ = mnu
        self.set_mnuone_()


    def setNnu(self, Nnu):
        self.Nnu_ = Nnu
        self.set_mnuone_()


    # @autojit
    def rho(self, a):
        # This returns the density at a normalized so that
        # we get nuh2 at a=0
        # (1 eV) / (Boltzmann constant * 1 kelvin) = 11 604.5193

        if (self.mnuone == 0):
            return self.Nnu_*7/8.*self.prefix0/a**4

        mnuOT = self.mnuone/(self.Tnu/a)*(1./ct.value(u'Boltzmann constant in eV/K')) #11604.5193

        # Here for massive we use 1*prefix (accounting for 1.015 in Tnu)
        # For massles we use Neff*prefix0 (so we account
        if self.degenerate:
            return 3*self.I.SevenEights(mnuOT)*self.prefix/a**4 + (self.Nnu_-3.014)*7/8.*self.prefix0/a**4
        else:
            return ((self.I.SevenEights(mnuOT)*self.prefix+(self.Nnu_-1.015)*7/8.*self.prefix0))/a**4



if __name__ == "__main__":
    A = NuDensity(CA.Tcmb, 3.046, 0.06)
    print(A.omnuh2today, '=including massless neutrinos')
    B = NuDensity(CA.Tcmb, 2.030, 0.00)
    print(B.omnuh2today, end=' ')
    print(A.omnuh2today-B.omnuh2today, '=excluding massless neutrinos')

    A = NuDensity(CA.Tcmb, 1.015, 60.)
    print(A.omnuh2today/1000., '=assuming very cold')
    A = NuDensity(CA.Tcmb, 1.015, 0.06)

    print(A.omnuh2today, '=assuming real temperature')

    print(1/(A.I.int_infty/A.Tnu*11604.5193*A.prefix))
    print(1/(A.I.int_infty/A.Tnu*11604.5193*A.prefix*(3.046/3.)))
