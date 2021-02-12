## This is Multi field cosmology
#Using the dynamical system approach

import numpy as np
from simplemc.models.LCDMCosmology import LCDMCosmology
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from scipy.optimize import newton

import sys

class MFieldCosmology(LCDMCosmology):
    def __init__(self):
        """Let's test with one field"""

        self.lna   = np.linspace(-6, 0, 500)
        self.z     = np.exp(-self.lna) - 1.
        self.zvals = np.linspace(0, 5, 200)

        self.phiini = 2.5
        self.psiini = 1.5

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


    def dphilnV(self, phi, psi):
        return 3./phi



    def dpsilnV(self, phi, psi):
        return 2./psi



    def RHS(self, x_vec , lna):
        xphi, xpsi, phi, psi, y = x_vec

        fac = np.sqrt(1.5)
        Pi  = -1.5*(1 + xphi**2 - xpsi**2 - y**2 )

        xphi_prime  = -xphi*(3 + Pi) - fac*y**2*self.dphilnV(phi, psi)
        xpsi_prime  = -xpsi*(3 + Pi) + fac*y**2*self.dpsilnV(phi, psi)
        phi_prime   = np.sqrt(6)*xphi
        psi_prime   = np.sqrt(6)*xpsi
        y_prime     = -Pi*y + fac*xphi*y*self.dphilnV(phi, psi) + fac*xpsi*y*self.dpsilnV(phi, psi)


        return [xphi_prime, xpsi_prime, phi_prime, psi_prime, y_prime]
        

    def solver(self, ini_y):
        ini_xphi = 0.01
        ini_xpsi = 0.5
        ini_phi  = self.phiini
        ini_psi  = self.psiini


        y0      = [ini_xphi, ini_xpsi, ini_phi, ini_psi,  10**(-ini_y)]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-5)
        return y_result


    def logatoz(self, func):
        #change functions from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])


    #shooting to find the initial condition
    def rfunc(self, ini_y):
        sol  = self.solver(ini_y).T
        Ophi = (sol[0][-1]**2 - sol[1][-1]**2 + sol[4][-1]**2)

        return (1.0-self.Om) - Ophi


    def set_ini(self):
        try:
            y0 = newton(self.rfunc, 4) #7.5
            x_vec = self.solver(y0).T
            self.do = 1
            self.hub_SF   = interp1d(self.lna, x_vec[3])
            #useful only for plotting purposes
            #self.hub_SF_z = self.logatoz(x_vec[3])

            pe  = x_vec[0]**2 - x_vec[1]**2 - x_vec[4]**2
            rho = x_vec[0]**2 - x_vec[1]**2 + x_vec[4]**2
            self.w_eos    = interp1d(self.lna, pe/rho)
            self.Ophi    = interp1d(self.lna, rho)
            #self.phi      = interp1d(self.lna, x_vec[2])
        except RuntimeError:
            print('--'*10)
            print('not working')
            print('--'*10)
            self.do = 2
            self.w_eos    = interp1d(self.lna, -1+np.zeros(len(self.lna)))
            #self.phi    = interp1d(self.lna, np.zeros(len(self.lna)))


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

    def Omegaphi(self, a):
        lna = np.log(a)
        return self.Ophi(lna)

    def Omegak(self, lna):
        return self.Oka(lna)


if __name__=='__main__':
    import matplotlib.pyplot as plt
    import numpy as np



    zl  = np.arange(0, 6.0, 0.05)
    lna = np.linspace(-6, 0, 500)


    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace': 0, 'height_ratios': [1, 1]}, figsize=(7,7))
    #plt.suptitle('$V=V_0\phi^2$, $A=%d$, $\phi_{ini}=%.1f$'%(a, c) , y=0.95, fontsize=15)

    #for i in np.linspace(2.0 ,2.2, 5):
    M = MFieldCosmology()
    y1=[M.w_de(1./(1+z))     for z in zl]
    y2=[M.Omegaphi(1./(1+z)) for z in zl]
    ax1.plot(zl, y1)
    ax2.plot(zl, y2)

    #ax1.set_ylabel('$w(z)$', fontsize=20, labelpad = -5)
    ax1.grid()
    #ax1.legend(loc='best')
    #ax1.axhline(-1, color='k')

    #for i in np.linspace(2.0 ,2.2, 5):
    #    M = MFieldCosmology(a=a, b=i, c=c)
    #    y2=[M.phi_de(1./(1+z))     for z in zl]
    #    ax2.plot(zl, y2)

    #ax2.set_ylabel('$\phi(z)$', fontsize=20, labelpad = -5)
    #ax2.set_xlabel("$z$", fontsize=20)
    ax2.grid()
    plt.tight_layout()
    #plt.savefig('tanh_field_2.pdf')
    plt.show()