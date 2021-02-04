


from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import *
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import numpy as np


class QuintomCosmology(LCDMCosmology):
    def __init__(self, vary_mquin=False, vary_mphan=False, vary_iniphi=False, vary_beta=False):

        """Is better to start the chains at masses equal one, othewise
        may take much longer"""

        self.vary_mquin  = vary_mquin
        self.vary_mphan  = vary_mphan
        self.vary_beta   = vary_beta
        self.vary_iniphi = vary_iniphi

        self.mquin    = 0 if (vary_mphan and (not vary_mquin)) else mquin_par.value
        self.mphan    = 0 if (vary_mquin and (not vary_mphan)) else mphan_par.value
        self.beta     = 0 if (not vary_beta) else beta_par.value
        self.iniphi   = iniphi_par.value

        self.lna   = np.linspace(-10, 0, 500)
        self.z     = np.exp(-self.lna) - 1.
        self.zvals = np.linspace(0, 3, 100)

        self.chatty = True

        self.exp_phi = 2.
        self.exp_psi = 2
        self.epsilon = 1

        LCDMCosmology.__init__(self, mnu=0)

        self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.vary_mquin)  : l.append(mquin_par)
        if (self.vary_mphan)  : l.append(mphan_par)
        if (self.vary_beta)   : l.append(beta_par)
        if (self.vary_iniphi) : l.append(iniphi_par)
        return l


    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name  == "mquin":
                self.mquin=p.value
            elif p.name== "mphan":
                self.mphan=p.value
            elif p.name== "iniphi":
                self.iniphi=p.value
            elif p.name=='beta':
                self.beta=p.value

        #Searches initial conditions and computes the Hubble function
        if self.chatty:
            print('-'*10)
            print ('mphi={}, mphan={}, beta={}, ini={}'.format(self.mquin, self.mphan, self.beta, self.iniphi))
        self.search_ini()

        return True


    def Vtotal(self, x, y, select):
        """Cuadratic potential and its derivatives wrt phi or psi"""
        if select == 0:
            Vtotal = 0.5*(self.mquin*x)**2  +  0.5*(self.mphan*y)**2 + self.beta*(x*y)**2
        elif select == 'phi':
            Vtotal = self.mquin**2*x + 2*self.beta*x*y**2
        elif select == 'psi':
            Vtotal = self.mphan**2*y + 2*self.beta*y*x**2
        return Vtotal


    def rhode(self, x_vec):
        quin, dotquin, phan, dotphan = x_vec
        Ode = dotquin**2 + dotphan**2 + self.Vtotal(quin, phan, 0)/(3*self.h**2)
        return Ode


    def eos(self, x_vec):
        quin, dotquin, phan, dotphan = x_vec
        w1 = dotquin**2 - dotphan**2 - self.Vtotal(quin, phan, 0)/(3*self.h**2)
        w2 = dotquin**2 - dotphan**2 + self.Vtotal(quin, phan, 0)/(3*self.h**2)
        return w1/w2



    def hubble(self, lna, x_vec=None, SF = True):
        a = np.exp(lna)
        if SF:
            Ode = self.rhode(x_vec)
        else:
            Ode = 1.0-self.Om

        return self.h*np.sqrt(self.Ocb/a**3 + self.Omrad/a**4 + Ode)


    def logatoz(self, func):
        # change function from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return np.interp(self.zvals, self.z[::-1], functmp[::-1])



    def RHS(self, x_vec, lna):
        factor = np.sqrt(6)*self.h
        quin, dotquin, phan, dotphan = x_vec
        hubble  = self.hubble(lna, x_vec)
        return [factor*dotquin/hubble, -3*self.h*dotquin + self.Vtotal(quin, phan, 'phi')/(factor*hubble),
                factor*dotphan/hubble, -3*self.h*dotphan - self.Vtotal(quin, phan, 'psi')/(factor*hubble)]



    def solver(self, quin0, dotquin0, phan_0, dotphan0):
        y0       = [quin0, dotquin0, phan_0, dotphan0]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result



    def calc_Ode(self, mid):
        "Select Quintess, Phantom or Quintom"
        #if (self.varymquin)   and (self.mphan == 0) and (self.beta ==0): quin0, phan0 = mid, 0
        #elif (self.varymphan) and (self.mquin == 0) and (self.beta ==0): quin0, phan0 = 0, mid
        #else : quin0, phan0 = mid, mid*self.iniphi
            #Still figuring out initial conditions for two fields
        quin0, phan0 = 10**mid, 0
        sol = self.solver(quin0 , 0.0, phan0, 0.0).T

        quin, dotq, phan, dotp = sol
        rho =  self.rhode(sol)[-1]
        Ode = rho*(self.h/self.hubble(0.0, [quin[-1], dotq[-1], phan[-1], dotp[-1]]))**2
        tol = (1- self.Ocb- self.Omrad) - Ode
        return sol, tol, Ode


    def calc_Ode2(self, phi0):
        "Select Quintess, Phantom or Quintom"
        #if (self.varymquin)   and (self.mphan == 0) and (self.beta ==0): quin0, phan0 = mid, 0
        #elif (self.varymphan) and (self.mquin == 0) and (self.beta ==0): quin0, phan0 = 0, mid
        #else : quin0, phan0 = mid, mid*self.iniphi
            #Still figuring out initial conditions for two fields

        quin0, phan0 = 10**phi0, 0
        sol = self.solver(quin0 , 0.0, phan0, 0.0).T

        quin, dotq, phan, dotp = sol
        rho = self.rhode(sol)[-1]
        Ode = rho*(self.h/self.hubble(0.0, [quin[-1], dotq[-1], phan[-1], dotp[-1]]))**2
        tol = (1- self.Ocb- self.Omrad) - Ode
        return tol


    def bisection(self):
        "Search for intial condition of phi such that \O_DE today is 0.7"
        lowphi, highphi = -2, 3 ##initial guess
        Ttol            = 5E-3
        mid = (lowphi + highphi)*0.5
        while (highphi - lowphi )*0.5 > Ttol:
            sol, tol_mid, Ode = self.calc_Ode(mid)
            tol_low = self.calc_Ode(lowphi)[1]
            if(np.abs(tol_mid) < Ttol):
                  if self.chatty: print('-done-, phi_0={}, error={}'.format(mid, tol_mid))
                  return mid, sol
            elif tol_low*tol_mid<0:
                highphi  = mid
            else:
                lowphi   = mid
            mid = (lowphi + highphi)*0.5

        if self.chatty: print ('No solution found!', mid, self.mquin, self.mphan, Ode)
        ##Check whether rho is constant or nearby
        grad = np.abs(np.gradient(self.rhode(sol))).max()
        mid  = -1 if grad < 1.0E-2 else 0
        return mid, sol


    def search_ini(self):
        #try:
        #    phi0 = newton(self.calc_Ode2, 1)
        #    print('phi0', phi0)
        #    sol, _,_ = self.calc_Ode(phi0)
        #except:
        #    phi0 = 0

        #It's slower than newton, but finds more solutions
        self.phi0, self.sol = self.bisection()

        if (self.phi0 == 0):
            if self.chatty: print('mass=', self.mquin, self.mphan, 'sol not found with', self.phi0)
            self.hub_SF   = interp1d(self.lna, np.zeros(len(self.lna)))
        elif (self.phi0 == -1):
            if self.chatty: print ('looks like LCDM', self.mquin, self.mphan)
            self.hub_SF   = interp1d(self.lna, self.hubble(self.lna, SF=False))
        else:
            self.hub_SF   = interp1d(self.lna, self.hubble(self.lna, self.sol))
        return True


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        lna = np.log(a)
        if (1./a-1 < self.zvals[-1]):
            hubble = (self.hub_SF(lna)/self.h)**2.
        else:
            hubble = self.hubble(self.lna, SF=False)
        return hubble


    #-------------------------

    def Hubble_a(self, a):
        lna = np.log(a)
        hubble = 100*self.hub_SF(lna)
        return hubble


    def tmp(self):
        if self.phi0==0:
            self.w_eos    = interp1d(self.lna, np.zeros(len(self.lna)))
        elif (self.phi0 == -1):
            self.w_eos    = interp1d(self.lna, -1*np.ones(len(self.lna)))
        else:
            self.w_eos    = interp1d(self.lna, self.eos(self.sol))



    def w_de(self, a):
        lna = np.log(a)
        return self.w_eos(lna)





