# This is a CDM cosmology with \phi
from pylab import *
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from LCDMCosmology import *


class QuintCosmology(LCDMCosmology):
    def __init__(self, varylam=True, varyV0=True, varyA=True, varyB=True):
        # two parameters: Om and h

        self.varyV0 = varyV0
        self.varylam = varylam
        self.varyA = varyA
        self.varyB = varyB

        V0f = 1.0  # ln(10**15)
        self.V0 = V0f*V0_par.value
        self.lam = lam_par.value
        self.B = B_par.value
        self.A = A_par.value

        LCDMCosmology.__init__(self, mnu=0)

        self.lna = linspace(-15, 0, 200)
        self.Cte = sqrt(3.0)*self.h

        # force caching
        self.updateParams([])

    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varylam):
            l.append(lam_par)
        if (self.varyV0):
            l.append(V0_par)
        if (self.varyA):
            l.append(A_par)
        if (self.varyB):
            l.append(B_par)
        return l

    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "lam":
                self.lam = p.value
        for p in pars:
            if p.name == "V0":
                self.V0 = p.value

        self.Ini_phi()
        return True

    def Pot(self, x, i):
        funct1 = (x-self.B)**2 + self.A
        funct2 = 2.0*(x-self.B) - self.lam*funct1
        if i == 0:
            return funct1*self.V0*exp(-self.lam*x)
        if i == 1:
            return funct2*self.V0*exp(-self.lam*x)
        else:
            print('wrong choice')
            stop()

            # x ->\phi, y->\dotphi
    def RHS(self, x_vec, lna):
        x, y = x_vec
        return [sqrt(3.0)*y/self.hub(lna, x_vec), -3*y - self.Pot(x, 1)*self.h/(self.Cte*self.hub(lna, x_vec))]

    def solver(self, x0):
        y0 = [x0, 0]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result

    def Ini_phi(self):
        lowr, highr = -50, 50
        tol, tol1 = 100, 100
        Ttol = 1e-2
        count = 0
        # search initial conditions
        if True:
            while (abs(tol) > Ttol):
                mid = (lowr+highr)/2.0
                sol = self.solver(mid)
                sol1 = self.solver(highr)

                Omegal = (0.5*sol[-1, 1]**2+self.Pot(sol[-1, 0],
                                                     0)/self.Cte**2)/self.hub(0.0, sol[-1])**2
                Omegal1 = (0.5*sol1[-1, 1]**2+self.Pot(sol1[-1, 0],
                                                       0)/self.Cte**2)/self.hub(0.0, sol1[-1])**2
                tol = (1.0-self.Om) - Omegal
                tol1 = (1.0-self.Om) - Omegal1

                if(abs(tol) < Ttol):
                    # print 'Omega', Omegal, 'phi', mid
                    # print 'reach tolerance', abs(tol), count
                    break
                else:
                    if(tol*tol1 > 0):
                        highr = mid
                    else:
                        lowr = mid

                count += 1
                if (count > 50):
                    'No solution found!'
                    break
        #sol =self.solver(0.0)
        self.sol = sol
        self.hubble = interp1d(self.lna, (self.hub(self.lna, sol.T))**2)
        return self.sol

    def hub(self, lna, x_vec):
        a = exp(lna)
        x, y = x_vec
        NuContrib = self.NuDensity.rho(a)/self.h**2
        return sqrt(0.5*y**2 + self.Pot(x, 0)/self.Cte**2 + self.Ocb/a**3 + self.Omrad/a**4 + NuContrib)

    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        lna = log(a)
        return self.hubble(lna)
