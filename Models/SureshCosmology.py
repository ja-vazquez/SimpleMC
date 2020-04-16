## This is phiCDM cosmology

import numpy as np
from LCDMCosmology import *
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from ParamDefs import alpha_par, Ok_par
import matplotlib.pyplot as plt

class SureshCosmology(LCDMCosmology):
    def __init__(self, varyalpha = True, varyOk=False):
        ## two parameters: Om and h
        """Is better to start the chains at masses equal one, othewise
        may take much longer"""

        self.varyalpha = varyalpha
        self.varyOk    = varyOk

        self.alpha = alpha_par.value
        self.Ok    = Ok_par.value

        self.a_in  = 5e-5
        self.ndots = 2000
        self.zvals = np.linspace(0, 5, 300)


        LCDMCosmology.__init__(self, mnu=0)
        self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyalpha): l.append(alpha_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name  == "alpha":
                self.alpha= p.value

        self.start()

        """
        dataHz = np.loadtxt('data/Hz_all.dat')
        redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
        plt.errorbar(redshifts, obs, errors, xerr=None,
                     color='purple', marker='o', ls='None',
                     elinewidth =2, capsize=5, capthick = 1, label='$Datos$')
        plt.xlabel(r'$z$')
        plt.ylabel(r'$H(z) [km/s Mpc^{-1}]$')
        plt.plot(self.zvals, [100*self.h*np.sqrt(self.hub_phi(1./(1+z))) for z in self.zvals])
        plt.title('mquin %f'%(self.alpha))
        plt.show()
        """
        return True


    def init_conditions(self, a):
        #phi0 = (2*self.alpha*(self.alpha + 2)/3.)**0.5*a**(3./(2 + self.alpha))
        #phi1 = 2*(2*self.alpha/((self.alpha + 2)*3))**0.5*a**(-3*self.alpha/(2*(2 + self.alpha)) + 1)
        phi0 = (self.alpha*(self.alpha + 2)/2.)**0.5*a**(4./(2 + self.alpha))
        phi1 = (2*self.alpha/((self.alpha + 2)))**0.5*a**((2-self.alpha)/(2 + self.alpha) + 1)
        return phi0, phi1




    def model(self, Phi, a, V1cap):
        phi, phiprime = Phi
        t1  = 1 - (a*phiprime)**2/6.
        t2  = V1cap*self.alpha*phi**(-(self.alpha + 1))
        t3  = V1cap*phi**(-self.alpha)
        num = t1*((-2.0*self.Omrad/a**4 - 3*self.Ocb/a**3 -4*self.Ok/a**2 +\
              t2/3.*a*phiprime - 2*t3)*a*phiprime + 2*t1*t2)
        den = 2.0*(self.Omrad/a**4 + self.Ocb/a**3 + self.Ok/a**2 + t3 / 3)*a**2.0
        phidoubleprime = num/den - phiprime/a
        return np.array((phiprime, phidoubleprime))



    def solvephicdm(self, V1cap):
        phi0, phiprime0 = self.init_conditions(self.a_in)
        a = np.linspace(self.a_in, 1, self.ndots)
        V1capprev = 0
        V1capnew  = V1cap
        while abs(V1capnew - V1capprev) > 1e-3:
            V1capprev = V1capnew
            Phi_model = odeint(self.model, [phi0, phiprime0], a, args=(V1capnew,))
            phi1      = Phi_model[-1, 0]
            phiprime1 = Phi_model[-1, 1]
            V1capnew  = 3*phi1**self.alpha*(1-self.Om-self.Ok-phiprime1**2/6.)
        return Phi_model, V1capnew


    def getphiandphiprimes(self, zdata, Phi_model):
        if type(zdata) is float or type(zdata) is np.float64:
            alen  = 1
            adata = round(1. / (1 + zdata), 4)
            adata_pos = adata * self.ndots - 1
            return Phi_model[int(adata_pos), 0], Phi_model[int(adata_pos), 1]
        else:
            alen  = len(zdata)
            adata = 1. / (1 + zdata)
            for i in range(alen):
                adata[i] = round(adata[i], 4)
            adata_pos = adata * self.ndots - 1
            phi      = np.zeros(alen)
            phiprime = np.zeros(alen)
            for i in range(alen):
                phi[i]      = Phi_model[int(adata_pos[i]), 0]
                phiprime[i] = Phi_model[int(adata_pos[i]), 1]
            return phi, phiprime



    def start(self):
        V1cap = 1.0
        sol, V1capnew = self.solvephicdm(V1cap)
        phi, phiprime = self.getphiandphiprimes(self.zvals, sol)

        a  = 1./(1+ self.zvals)
        t1 = 1 - (a * phiprime)**2/6.
        t3 = V1capnew*phi**(-self.alpha)
        H2 = (self.Omrad/a**4 + self.Ocb/a**3 +self.Ok/a**2 + t3/3.)/t1
        self.hub_phi   = interp1d(a, H2)
        return True





    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        #NuContrib = self.NuDensity.rho(a)/self.h**2
        if a< 1./(1+ self.zvals[-1]):
            hubble = (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+(1.0-self.Om-self.Ok))
        else:
            hubble = self.hub_phi(a)
        return hubble
