## This is phiCDM cosmology

from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import alpha_par, Ok_par
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import odeint

import matplotlib.pyplot as plt
#TODO In construction, not updated on github but in my computer

class SureshCosmology(LCDMCosmology):
    def __init__(self, varyalpha = True, varyOk=False):
        ## two parameters: Om and h
        """Is better to start the chains at masses equal one, othewise
        may take much longer"""

        self.varyalpha = varyalpha
        self.varyOk    = varyOk

        self.alpha = alpha_par.value
        self.Ok    = Ok_par.value

        #Compute Vcap either analytically or find it
        self.fixVcap = False

        #Pick one for initial conditions
        #self.inic  = 'Matter'
        #self.a_in  = 5e-3
        self.inic  = 'Radiation'
        self.a_in  = 5e-5

        self.ndots = 1000
        self.zvals = np.linspace(0, 4, 50)
        self.ascal = np.linspace(self.a_in, 1, self.ndots)


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
        return True


    def init_conditions(self, a):
        if self.inic == 'Matter':
            phi0 = (2*self.alpha*(self.alpha + 2)/3.)**0.5*a**(3./(2 + self.alpha))
            phi1 = 2*(2*self.alpha/((self.alpha + 2)*3))**0.5*a**(-3*self.alpha/(2*(2 + self.alpha)))
        if self.inic == 'Radiation':
            phi0 = (0.5*self.alpha*(self.alpha + 2))**0.5*a**(4./(2 + self.alpha))
            phi1 = (2*self.alpha/((self.alpha + 2)))**0.5*a**((2-self.alpha)/(2 + self.alpha))
        else:
            print('??')
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
        V1capprev = 0
        V1capnew  = V1cap
        if self.fixVcap:
            Phi_model = odeint(self.model, [phi0, phiprime0], self.ascal, args=(V1cap,))
            return Phi_model, V1cap
        while abs(V1capnew - V1capprev) > 1e-3:
            V1capprev = V1capnew
            Phi_model = odeint(self.model, [phi0, phiprime0], self.ascal, args=(V1capnew,))
            phi1      = Phi_model[-1, 0]
            phiprime1 = Phi_model[-1, 1]
            V1capnew  = 3*phi1**self.alpha*(1-self.Om-self.Ok-phiprime1**2/6.)
        print ('Vcap -- ', V1capnew, ((self.alpha+6)/(self.alpha+2))*(0.5*self.alpha*(self.alpha+2))**(0.5*self.alpha))
        return Phi_model, V1capnew



    def start(self):
        if self.fixVcap:
            V1cap = ((self.alpha+6)/(self.alpha+2))*(0.5*self.alpha*(self.alpha+2))**(0.5*self.alpha)
            #kappa factor with units
            V1cap /= 1.42
        else:
            #Initial Guess
            V1cap = 1.8
        sol, V1capnew = self.solvephicdm(V1cap)
        phi, phiprime = sol.T

        a  = self.ascal
        t1 = 1 - (a * phiprime)**2/6.
        t3 = V1capnew*phi**(-self.alpha)
        H2 = (self.Omrad/a**4 + self.Ocb/a**3 +self.Ok/a**2 + t3/3.)/t1
        self.hub_phi   = interp1d(a, H2)

        wde = ((a*phiprime)**2/6. - V1capnew*phi**(-self.alpha)/3.)/\
              ((a*phiprime)**2/6. + V1capnew*phi**(-self.alpha)/3.)
        self.w_eos = interp1d(a, wde)
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



    def Hubble_a(self, a):
        hubble = 100*self.h*np.sqrt(self.hub_phi(a))
        return hubble



    def w_de(self, a):
        return self.w_eos(a)
