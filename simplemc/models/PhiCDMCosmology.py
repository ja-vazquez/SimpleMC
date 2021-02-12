## This is phiCDM cosmology

import numpy as np
from simplemc.models.LCDMCosmology import LCDMCosmology
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from scipy.optimize import newton
from simplemc.cosmo.paramDefs import phialp_par, philam_par, phimu_par, \
                      phibeta_par, Ok_par
import sys

class PhiCosmology(LCDMCosmology):
    def __init__(self, varyalpha=False, varybeta=False, varyilam=False,\
                       varymu=False, varyOk=False,
                       alpha=1, beta=1, mu=1, ilam=1, eps=1, curv=0):
        """Is better to start the chains at masses equal one, othewise
        may take much longer"""

         #eps => 1=Quintes, -1=Panthom

        self.varyOk      = varyOk
        self.varymu      = varymu
        self.varyilam    = varyilam
        self.varybeta    = varybeta
        self.varyalpha   = varyalpha

        self.Ok     = curv   #Ok_par.value
        self.alpha  = alpha
        self.beta   = beta   #phibeta_par.value
        self.mu     = mu     #phimu_par.value
        self.ilam   = ilam   #philam_par.value
        self.eps    = eps

        self.lna   = np.linspace(-6.5, 0, 500)
        self.z     = np.exp(-self.lna) - 1.
        self.zvals = np.linspace(0, 5, 200)

        LCDMCosmology.__init__(self, mnu=0)
        self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varyOk)      : l.append(Ok_par)
        if (self.varymu)      : l.append(phimu_par)
        if (self.varyilam)    : l.append(philam_par)
        if (self.varybeta)    : l.append(phibeta_par)
        if (self.varyalpha)   : l.append(phialp_par)

        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name    == "philam":
                self.ilam = p.value
            elif p.name  == "phialp":
                self.alpha= p.value
            elif p.name  == "phibeta":
                self.beta = p.value
            elif p.name  == "phimu":
                self.mu   = p.value
            elif p.name  == "Ok":
                self.Ok   = p.value
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


    def MGama(self, lam):
        if type(self.mu) == type('st'):
            tmp = self.alpha**2/lam**2
            if self.beta == 0:  #cosh
                return tmp
            elif self.beta == -1:
                return 0.5*(1+ tmp)
            elif self.beta == 1:
                return  0.5*(1- tmp)
            else: sys.exit('wrong potential')
        else:
            if self.beta==0:    #pow
                return 1 - 1./self.mu
            else:
                if self.mu == 0:
                    if self.alpha== 1:  #exp
                        return 1
                    elif self.alpha== 2: #exp_pow2
                        return 1+ 2*self.beta/lam**2
                    else:           #'exp_pow_a'
                        #fac = self.alpha*(self.alpha-1)*self.beta/lam**2
                        #return 1- fac*(-lam/(self.alpha*self.beta))**((self.alpha-2)/(self.alpha-1))
                        fac = (self.alpha-1)/(self.alpha*self.beta)
                        return 1+ fac*(-lam/(self.alpha*self.beta))**(self.alpha/(1-self.alpha))
                elif self.mu==2 and self.alpha==2:   #pow2_exp_pow2
                    fac= lam*np.sqrt(lam**2-16*self.beta)
                    return 1 + (4*self.beta/lam**2)*(-16*self.beta-lam**2+fac)/(-16*self.beta-2*lam**2+2*fac)
                else:       #pow_exp
                    if self.alpha ==1:
                        if self.mu ==0: return  1
                        else:
                            return 1 - (1+ self.beta/lam)**2/self.mu



    def RHS(self, x_vec , lna):
        wphi, Ophi, lam, Ok, hub = x_vec

        Mgamma = self.MGama(lam)
        term   = np.sqrt(3*np.abs((1+wphi)*self.eps*Ophi))
        Pi     = -1.5*(-Ok/3. + wphi*Ophi + 1)

        Ophi_prime  = -3*Ophi*(1+ wphi + 2*Pi/3.)
        Ok_prime    = -3*Ok*(1 - 1./3. + 2*Pi/3.)
        wphi_prime  = -(1-wphi)*(3*(1+wphi) - self.eps*lam*term)
        lam_prime   = -self.eps*lam**2*(Mgamma -1)*term
        hub_prime   = hub*Pi

        return [wphi_prime, Ophi_prime, lam_prime, Ok_prime, hub_prime]
        

    def solver(self, ini_Ophi):
        if type(self.mu) == type('st'):
            if self.beta == 0:                  #cosh
                ini_lam = -self.alpha*np.tanh(self.alpha*self.ilam)
            elif self.beta == -1:
                ini_lam = -self.alpha/np.tanh(0.5*self.alpha/self.ilam)
            elif self.beta == 1:
                ini_lam = self.alpha*np.tan(0.5*self.alpha*self.ilam)
            else: sys.exit('wrong potential')
        else:
            if self.beta==0:                        #pow
                ini_lam= -self.mu*self.ilam
            else:
                if self.mu == 0:
                    if self.alpha== 1:              #exp
                        ini_lam= -self.beta
                    elif self.alpha== 2:            #exp_pow2
                        ini_lam= -2*self.beta*self.ilam
                    else:                           #'exp_pow_a'
                        ini_lam= -self.alpha*self.beta*self.ilam**(self.alpha-1)
                elif self.mu==2 and self.alpha==2:   #pow2_exp_pow2
                        ini_lam= -self.ilam #2*(self.ilam + self.beta/self.ilam)
                else:                               #pow_exp
                    if self.alpha == 1:
                        ini_lam= -(self.mu*self.ilam + self.beta)


        #we'll use the sign of lambda to describe either quint or phant
        #ini_lam=-ini_lam
        self.eps = np.sign(ini_lam)
        ini_lam  =np.abs(ini_lam)
        ini_wphi = -1 + 1.0e-4*self.eps


        ini_hub = 100*self.h*self.Om**0.5*np.exp(-1.5*self.lna[0])
        ini_Ok  = self.Ok*np.exp(-2*self.lna[0])/(self.Om**0.5*np.exp(-1.5*self.lna[0]))**2
        y0      = [ini_wphi, 10**(-ini_Ophi), ini_lam, ini_Ok, ini_hub]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-5)
        return y_result


    def logatoz(self, func):
        #change functions from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])


    def rfunc(self, ini_Ophi0):
        #returns lambda that's solution
        sol  = self.solver(ini_Ophi0).T
        return (1.0-self.Om-self.Ok) - sol[1][-1]


    def set_ini(self):
        try:
            Ophi0 = newton(self.rfunc, 8)
            x_vec = self.solver(Ophi0).T
            self.do = 1
            self.hub_SF   = interp1d(self.lna, x_vec[4])
            #useful only for plotting purposes
            #self.hub_SF_z = self.logatoz(x_vec[3])
            self.w_eos    = interp1d(self.lna, x_vec[0])
            #self.Ophi    = interp1d(self.lna, x_vec[1])
            #self.Oka      = interp1d(self.lna, x_vec[3])
        except RuntimeError:
            self.do = 2
            #self.w_eos    = interp1d(self.lna, -1+np.zeros(len(self.lna)))



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

    def Omegaphi(self, lna):
        return self.Ophi(lna)

    def Omegak(self, lna):
        return self.Oka(lna)

