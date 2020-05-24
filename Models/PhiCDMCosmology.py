## This is phiCDM cosmology

import numpy as np
from LCDMCosmology import *
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from ParamDefs import phialp_par, philam_par, phimu_par, \
                      phibeta_par, epsilon_par, Ok_par, mu_par
from scipy.optimize import newton

class PhiCosmology(LCDMCosmology):
    def __init__(self, varyalpha=False, varybeta=False, varyilam=False,\
                       varyepsilon=False, varymu=False, varyOk=False,
                       alpha=1, beta=1, mu=1, eps=1, ilam=1):
        """Is better to start the chains at masses equal one, othewise
        may take much longer"""


        self.varyOk      = varyOk
        self.varymu      = varymu
        self.varyilam    = varyilam
        self.varybeta    = varybeta
        self.varyalpha   = varyalpha
        self.varyepsilon = varyepsilon

        #print('-'*10, 'very sensitive to initial conditions')
        #print('-'*10, 'hence we expect plenty of warning')

        self.Ok     = Ok_par.value
        self.eps    = eps   #1=Quintes, -1=Panthom
        self.alpha  = alpha
        self.beta   = beta
        self.mu     = mu
        self.ilam   = ilam


        self.lna   = np.linspace(-6, 0, 500)
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
        if (self.varyepsilon) : l.append(epsilon_par)

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
            elif p.name  == "epsilon":
                self.eps  = p.value
            elif p.name  == "Ok":
                self.Ok   = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False

        self.ini_gamma = 1.0e-4*self.eps
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
        if self.mu == type('st'):
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
                        fac = (self.alpha-1)/(self.beta*self.alpha)
                        return 1 + fac*(lam/np.abs(self.beta*self.alpha))**(-self.alpha/(self.alpha-1))
                elif self.mu==2 and self.alpha==2:   #pow2_exp_pow2
                    pass
                else:       #pow_exp
                    if self.alpha ==1:
                        if self.mu ==0: return  1
                        else:
                            return 1 - (1+ self.beta/lam)**2/self.mu



    def RHS(self, x_vec , lna):
        gamma, Ophi, lam, Ok, hub = x_vec

        Mgamma = self.MGama(lam)
        term   = np.sqrt(3*np.abs(gamma*Ophi))
        wphi   = self.eps*gamma -1
        Pi     = -1.5*(-Ok/3. + wphi*Ophi + 1)
        delta  = 1

        #gamma_prime = (2 - self.eps*gamma)*(self.eps*delta*lam*term -3*gamma)
        gamma_prime = (2 - self.eps*gamma)*(delta*lam*term -3*gamma)
        Ophi_prime  = -Ophi*(3*self.eps*gamma + 2*Pi)
        lam_prime   = -delta*lam**2*(Mgamma -1)*term
        Ok_prime    = -2*Ok - 2*Pi*Ok
        hub_prime   = hub*Pi

        return [gamma_prime, Ophi_prime, lam_prime, Ok_prime, hub_prime]
        

    def solver(self, ini_Ophi):
        if self.mu == type('st'):
            if self.beta == 0:                  #cosh
                ini_lam = self.alpha*np.tanh(self.alpha*self.ilam)
            elif self.beta == -1:
                ini_lam = self.alpha/np.tanh(0.5*self.alpha/self.ilam)
            elif self.beta == 1:
                ini_lam = -self.alpha*np.arctan(0.5*self.alpha*self.ilam)
            else: sys.exit('wrong potential')
        else:
            ini_lam=self.ilam
            """
            if self.beta==0:                        #pow
                ini_lam= self.mu*self.ilam
            else:
                if self.mu == 0:
                    if self.alpha== 1:              #exp
                        ini_lam= self.beta
                    elif self.alpha== 2:            #exp_pow2
                        ini_lam= 2*self.beta*self.ilam
                    else:                           #'exp_pow_a'
                        ini_lam= self.alpha*self.beta*self.ilam#**(self.alpha-1)
                elif self.mu==2 and self.alpha==2:   #pow2_exp_pow2
                        pass
                else:                               #pow_exp
                    if self.alpha == 1:
                        ini_lam= self.mu*self.ilam + self.beta

            ini_lam = np.abs(ini_lam)
            """


        ini_hub = 100*self.h*self.Ocb**0.5*np.exp(-1.5*self.lna[0])
        ini_Ok  = self.Ok*np.exp(-2*self.lna[0])/(self.Ocb**0.5*np.exp(-1.5*self.lna[0]))**2
        y0      = [self.ini_gamma, 10**(-ini_Ophi), ini_lam, ini_Ok, ini_hub]
        #if self.varyOk: y0.append(ini_Ok)
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
        return (1-self.Om) - sol[1][-1]


    def set_ini(self):
        try:
            Ophi0 = newton(self.rfunc, 8)
            x_vec = self.solver(Ophi0).T
            self.do = 1
            self.hub_SF   = interp1d(self.lna, x_vec[4])
            #self.hub_SF_z = self.logatoz(x_vec[3])
            self.w_eos    = interp1d(self.lna, x_vec[0])
        except RuntimeError:
            if np.abs(self.alpha) < 0.02: self.do = 0
            else:
                self.w_eos    = interp1d(self.lna, np.zeros(len(self.lna)))
                print('troubles', 'a=',self.alpha, 'b=',self.beta, 'l=', self.ilam, 'm=', self.mu)
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
            return hubble
        else:
            return self.Ocb/a**3


    def w_de(self, a):
        lna = np.log(a)
        return self.eps*self.w_eos(lna)-1

