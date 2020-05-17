## This is phiCDM cosmology

import numpy as np
from LCDMCosmology import *
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from ParamDefs import phialp_par, philam_par, \
                      phibeta_par, epsilon_par, Ok_par
from scipy.optimize import newton

class PhiCosmology(LCDMCosmology):
    def __init__(self, poten='pow', varyalpha=False, varybeta=False,\
                 varyilam=False, varyepsilon=False, varyOk=False):
        """Is better to start the chains at masses equal one, othewise
        may take much longer"""


        self.poten    = poten  #Pow-law = pow, Exp = exp

        self.varyOk      = varyOk
        self.varyilam    = varyilam
        self.varybeta    = varybeta
        self.varyalpha   = varyalpha
        self.varyepsilon = varyepsilon

        self.Ok     = Ok_par.value
        self.alpha  = phialp_par.value
        self.beta   = phibeta_par.value
        self.ilam   = philam_par.value
        self.eps    = epsilon_par.value     #1=Quintes, -1=Panthom

        self.lna   = np.linspace(-5, 0, 400)
        self.z     = np.exp(-self.lna) - 1.
        self.zvals = np.linspace(0, 5, 200)

        self.ini_gamma = 1.0e-4

        LCDMCosmology.__init__(self, mnu=0)
        self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyOk)     : l.append(Ok_par)
        if (self.varyilam)   : l.append(philam_par)
        if (self.varybeta)   : l.append(phibeta_par)
        if (self.varyalpha)  : l.append(phialp_par)
        if (self.varyquipha) : l.append(epsilon_par)

        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name  == "philam":
                self.ilam= p.value
            elif p.name  == "phialp":
                self.alpha= p.value
            elif p.name  == "phibeta":
                self.beta= p.value
            elif p.name  == "epsilon":
                self.eps = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False
        self.set_ini()


        """ Testing
        import matplotlib.pyplot as plt
        dataHz = np.loadtxt('data/Hz_all.dat')
        redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
        plt.errorbar(redshifts, obs, errors, xerr=None,
                     color='purple', marker='o', ls='None',
                     elinewidth =2, capsize=5, capthick = 1, label='$Datos$')
        plt.xlabel(r'$z$')
        plt.ylabel(r'$H(z) [km/s Mpc^{-1}]$')
        plt.plot(self.zvals, [67.4*np.sqrt(self.RHSquared_a(1./(1+z))) for z in self.zvals])
        #plt.plot(self.zvals, [67.4*np.sqrt(self.hubble(1./(1+z))) for z in self.zvals])
        plt.title('mquin %f'%(self.alpha))
        plt.show()
        """
        return True


    def MG(self, lam):
        if self.poten == 'pow':
            return (self.alpha-1)/self.alpha*lam**2
        elif self.poten == 'exp':
            if self.alpha == 1:
                return 1
            else:
                fac = (self.alpha-1)/(self.beta*self.alpha)
                return 1 + fac*(lam/np.abs(self.beta*self.alpha))**(-self.alpha/(self.alpha-1))


    def RHS(self, x_vec, lna):
        gamma, Ophi, lam, Ok, hub = x_vec

        Mgamma= self.MG(lam)
        Pi    = -1.5*(-Ok/3. + 1 + (self.eps*gamma -1)*Ophi)

        Ophi_prime  = 3*Ophi*((1-self.eps*gamma)*(1-Ophi) - Ok/3.)
        gamma_prime = (2 - self.eps*gamma)*(-3*gamma + lam*np.sqrt(3*gamma*Ophi))
        Ok_prime    = -2*Ok - 2*Pi*Ok
        hub_prime   = -1.5*hub*(1 + (self.eps*gamma-1)*Ophi)

        if self.alpha ==0 or self.beta ==0: lam_prime =0
        else: lam_prime   = -np.sqrt(3)*lam**2*(Mgamma -1)*np.sqrt(gamma*Ophi)

        return [gamma_prime, Ophi_prime, lam_prime, Ok_prime, hub_prime]



    def solver(self, ini_Ophi):
        phi0 = self.ilam
        if self.varyilam:
            ini_lam  =  self.alpha/phi0
        else:
            if self.poten == 'pow':
                ini_lam  =  self.alpha/phi0
            elif self.poten == 'exp':
                if self.alpha == 1:
                    ini_lam = self.beta
                else:
                    ini_lam  = self.beta*self.alpha*phi0**(self.alpha-1)

        ini_lam  = np.abs(ini_lam)
        ini_hub  = 100*self.h*self.Ocb**0.5*np.exp(-1.5*self.lna[0])
        ini_Ok   = self.Ok*np.exp(-2*self.lna[0])/(self.Ocb**0.5*np.exp(-1.5*self.lna[0]))**2
        y0       = [self.ini_gamma, 10**(-ini_Ophi), ini_lam, ini_Ok, ini_hub]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result


    def logatoz(self, func):
        #change functions from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])


    def rfunc(self, ini_Ophi0):
        #returns lambda that's solution
        sol  = self.solver(ini_Ophi0).T
        return (1-self.Om) - sol[1][-1]


    def set_ini(self):
        try:
            Ophi0 = newton(self.rfunc, 6)
            x_vec = self.solver(Ophi0).T
            self.do = 1
            self.hub_SF   = interp1d(self.lna, x_vec[4])
            #self.hub_SF_z = self.logatoz(x_vec[3])
            self.w_eos    = interp1d(self.lna, x_vec[0])
        except RuntimeError:
            if np.abs(self.alpha) < 0.02: self.do = 0
            else:
                print('troubles', self.alpha, self.beta, self.ilam)
                self.do = 2



    def hubble(self, a):
        #NuContrib = self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3 + self.Ok/a**2 + self.Omrad/a**4 + (1.0-self.Om-self.Ok))


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        lna = np.log(a)
        if self.do==0:
            return self.hubble(a)
        elif self.do==1:
            if (lna > self.lna[0]):
                hubble = (self.hub_SF(lna)/100./self.h)**2.
            else:
                hubble = self.hubble(a)
        else:
            return 1.
        return hubble


    def w_de(self, a):
        lna = np.log(a)
        return self.eps*self.w_eos(lna)-1

