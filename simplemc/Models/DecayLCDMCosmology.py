# This is a CDM cosmology with a decaying
# dark matter component.
##

import sys
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from LCDMCosmology import LCDMCosmology
from ParamDefs import xfrac_par, lambda_par


class DecayLCDMCosmology(LCDMCosmology):
    # note that if we don't varyOr, it will be set so that
    # density at early a is zero.
    def __init__(self, varylam=True, varyxfrac=True, xfrac=xfrac_par.value):

        self.varylam   = varylam
        self.varyxfrac = varyxfrac

        self.lam = lambda_par.value
        self.xfrac = xfrac

        LCDMCosmology.__init__(self)

        self.logar  = np.linspace(0.0, -7.1, 100)
        self.ilogar = self.logar[::-1]

        # force caching
        self.updateParams([])


    # my free parameters. We add lam, xfrac
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varylam):   l.append(lambda_par)
        if (self.varyxfrac): l.append(xfrac_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "lambda":
                self.lam = p.value
            if p.name == "xfrac":
                self.xfrac = p.value
        self.SolveEq()
        # and updated with relevant rd
        self.setrd(self.rd_func_(
            self.Obh2, self.Ocbh2_early, self.Omnuh2, self.Nnu()))
        assert(abs(self.RHSquared_a(1.0)-1) < 1e-4)
        return True


    def H2_rxrr_a(self, a, rx, rr):
        NuContrib = self.NuDensity.rho(a)/self.h**2
        return self.Ocb_std/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om-self.Or)+rx/a**3+rr/a**4


    def RHS_(self, y, lna):
        # we are solving rx, so that rhox=rx/a**3 and rhor=rr/a**4
        ##
        a  = np.exp(lna)
        H2 = self.H2_rxrr_a(a, y[0], y[1])
        H  = np.sqrt(abs(H2))
        factor = self.lam*y[0]/H
        return np.array([- factor, + factor * a])


    def FuncMin_(self, x):
        fractoday, self.Or = x
        self.Odm_dec  = self.Odm*fractoday
        self.Odm_ndec = self.Odm*(1-fractoday)
        self.Ocb_std  = self.Ob+self.Odm_ndec
        yinit = np.array([self.Odm_dec, self.Or])
        sol   = odeint(self.RHS_, yinit, self.logar)
        rxe, rre = sol[-1, :]
        # we want Ore early be as small as possible
        eps = rre**2
        # we want early frac to be xfrac
        eps += ((rxe/(rxe+self.Odm_ndec))-self.xfrac)**2
        return eps


    def SolveEq(self):
        self.Odm = self.Ocb-self.Obh2/(self.h**2)
        self.Ob  = self.Ocb-self.Odm
        if (self.lam == 0):
            self.fractoday, self.Or = self.xfrac, 0.0
        if (self.xfrac == 1.0):
            res = minimize(lambda x: self.FuncMin_(
                [1.0, x[0]]), [0.001], tol=1e-5)
            self.Or = res.x[0]
            self.fractoday = 1.0
        else:
            res = minimize(self.FuncMin_, [self.xfrac, 0.001], tol=1e-5)
            self.fractoday, self.Or = res.x
            # print fmin.__doc__
            # print "lam=",self.lam
            # print "res=",res
            # stop
        # stupid interp1d doesn't take inverses
        self.Ocb_std = self.Ob+self.Odm*(1-self.fractoday)
        self.Odm_dec = self.Odm*self.fractoday
        yinit = np.array([self.Odm_dec, self.Or])
        sol   = odeint(self.RHS_, yinit, self.logar)


        self.sol = sol
        self.rx  = interp1d(self.ilogar, sol[::-1, 0])
        self.rr  = interp1d(self.ilogar, sol[::-1, 1])
        # take early time solution
        self.Ocbh2_early = (self.Ocb_std+sol[-1, 0])*self.h**2


    def RHSquared_a(self, a):
        lna = np.log(a)
        return self.H2_rxrr_a(a, self.rx(lna), self.rr(lna))


    def WangWangVec(self):
        print("no WW with Decay")
        sys.exit(1)
        return None

    # this returns the "SimpleCMB" variables in a vec
    def CMBSimpleVec(self):
        zstar = 1090.
        Dastar = self.Da_z(zstar)*self.c_/(self.h*100)
        return np.array([self.Obh2, self.Ocbh2_early, Dastar/self.rd])

    def Om_z(self, a):
        lna = np.log(a)
        return self.Ocb_std/a**3 + self.rx(lna)/a**3
