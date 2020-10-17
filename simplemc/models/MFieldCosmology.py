## This is Multi field cosmology

import numpy as np
from simplemc.models.LCDMCosmology import LCDMCosmology
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from scipy.optimize import newton

import sys

class MFieldCosmology(LCDMCosmology):
    def __init__(self, a=2., b=2.1, c=2.):
        """Let's test with one field"""

        self.lna   = np.linspace(-6, 0, 500)
        self.z     = np.exp(-self.lna) - 1.
        self.zvals = np.linspace(0, 5, 200)

        #self.eps   = -1
        self.amp   = a
        self.shift = b
        self.phiini= c

        LCDMCosmology.__init__(self, mnu=0)
        self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False

        self.set_ini()
        return True



    def RHS(self, x_vec , lna):
        wphi, Ophi, phi, hub = x_vec

        lam = 3./phi
        eps = self.amp*np.tanh( ( phi - self.shift)/0.5)

        term   = np.sqrt(3*Ophi*np.abs((1+wphi)/eps))
        Pi     = -1.5*(wphi*Ophi + 1)

        Ophi_prime  = -3*Ophi*(1+ wphi + 2*Pi/3.)
        wphi_prime  = -(1-wphi)*(3*(1+wphi) - term*lam*eps)
        phi_prime   = term
        hub_prime   = hub*Pi

        return [wphi_prime, Ophi_prime, phi_prime, hub_prime]
        

    def solver(self, ini_Ophi):
        ini_wphi = -1 + 1.0e-4*np.sign(self.amp) #self.eps
        ini_phi  = self.phiini
        ini_hub  = 100*self.h*self.Om**0.5*np.exp(-1.5*self.lna[0])

        #print ([ini_wphi, 10**(-ini_Ophi), ini_phi, ini_hub])
        y0      = [ini_wphi, 10**(-ini_Ophi), ini_phi, ini_hub]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-5)
        return y_result


    def logatoz(self, func):
        #change functions from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])


    def rfunc(self, ini_Ophi0):
        #returns Ophi that's solution
        sol  = self.solver(ini_Ophi0).T
        return (1.0-self.Om) - sol[1][-1]


    def set_ini(self):
        try:
            Ophi0 = newton(self.rfunc, 7) #7.5
            x_vec = self.solver(Ophi0).T
            self.do = 1
            self.hub_SF   = interp1d(self.lna, x_vec[3])
            #useful only for plotting purposes
            #self.hub_SF_z = self.logatoz(x_vec[3])
            self.w_eos    = interp1d(self.lna, x_vec[0])
            self.Ophi    = interp1d(self.lna, x_vec[1])
            self.phi      = interp1d(self.lna, x_vec[2])
        except RuntimeError:
            print('--'*10)
            print('not working')
            print('--'*10)
            self.do = 2
            self.w_eos    = interp1d(self.lna, -1+np.zeros(len(self.lna)))
            self.phi    = interp1d(self.lna, np.zeros(len(self.lna)))


    def hubble(self, a):
        #NuContrib = self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3 + self.Ok/a**2 + self.Omrad/a**4 + (1.0-self.Om-self.Ok))



    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        lna = np.log(a)
        if self.do==2:
            return self.Ocb/a**3
        else:
            if (lna > self.lna[0]):
                hubble = (self.hub_SF(lna)/100./self.h)**2.
            else:
                hubble = self.hubble(a)
            return hubble



    def w_de(self, a):
        lna = np.log(a)
        return self.w_eos(lna)

    def phi_de(self, a):
        lna = np.log(a)
        return self.phi(lna)

    def Omegaphi(self, lna):
        return self.Ophi(lna)

    def Omegak(self, lna):
        return self.Oka(lna)


if __name__=='__main__':
    import matplotlib.pyplot as plt
    import numpy as np


    a, c, = -2, 1.9
    zl  = np.arange(0, 3.0, 0.05)
    lna = np.linspace(-6, 0, 500)


    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace': 0, 'height_ratios': [2, 1]}, figsize=(7,7))
    plt.suptitle('$V=V_0\phi^2$, $A=%d$, $\phi_{ini}=%.1f$'%(a, c) , y=0.95, fontsize=15)

    for i in np.linspace(2.0 ,2.2, 5):
        M = MFieldCosmology(a=a, b=i, c=c)
        y1=[M.w_de(1./(1+z))     for z in zl]
        ax1.plot(zl, y1, label='b=%.2f'%i)

    ax1.set_ylabel('$w(z)$', fontsize=20, labelpad = -5)
    ax1.grid()
    ax1.legend(loc='best')
    ax1.axhline(-1, color='k')

    for i in np.linspace(2.0 ,2.2, 5):
        M = MFieldCosmology(a=a, b=i, c=c)
        y2=[M.phi_de(1./(1+z))     for z in zl]
        ax2.plot(zl, y2)

    ax2.set_ylabel('$\phi(z)$', fontsize=20, labelpad = -5)
    ax2.set_xlabel("$z$", fontsize=20)
    ax2.grid()
    plt.tight_layout()
    plt.savefig('tanh_field_2.pdf')
    plt.show()