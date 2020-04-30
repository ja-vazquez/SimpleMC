## This is CDM cosmology with w, wa and Ok


import math as N
import numpy as np
import matplotlib.pyplot as plt
from LCDMCosmology import *

class QuintomCosmology(LCDMCosmology):
    def __init__(self, varymquin=False, varymphan=False, varyiniphi=False, varybeta=False, varydelta=False):
        ## two parameters: Om and h

        """Is better to start the chains at masses equal one, othewise
        may take much longer"""

        self.varymquin = varymquin
        self.varymphan = varymphan
        self.varyiniphi= varyiniphi
        self.varybeta = varybeta
        self.varydelta = varydelta

        self.mquin    = 0 if (varymphan and (not varymquin)) else mquin_par.value
        self.mphan    = 0 if (varymquin and (not varymphan)) else mphan_par.value
        self.beta     = 0 if (not varybeta) else beta_par.value
        self.iniphi   = iniphi_par.value
        self.delta    = delta_par.value

        LCDMCosmology.__init__(self, mnu=0)
        # it started from -15
        self.lna   = linspace(-15, 0, 300)
        self.z     = np.exp(-self.lna) - 1.
        self.zvals = np.linspace(0, 3, 100)

        self.cte   = 3.0*self.h**2
        self.n     = 2
        self.m     = 2
        self.comments = False

        self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varymquin)  : l.append(mquin_par)
        if (self.varymphan)  : l.append(mphan_par)
        if (self.varyiniphi) : l.append(iniphi_par)
        if (self.varybeta)   : l.append(beta_par)
        if (self.varydelta)  : l.append(delta_par)
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
            elif p.name=='delta':
                self.delta=p.value

        self.search_ini()
        print 'phi, phan, beta, delta, ini', self.mquin, self.mphan, self.beta, self.delta, self.iniphi
        if False:
        #if self.mid != 0:
            "Useful for testing purposes"
            quin, dquin, phan, dphan = self.solution

            w1 = 0.5*dquin**2 - self.Vtotal(quin, phan, 0)/self.cte - 0.5*dphan**2  - self.Vtotal(quin, phan, 0)/self.cte
            w2 = 0.5*dquin**2 + self.Vtotal(quin, phan, 0)/self.cte - 0.5*dphan**2  + self.Vtotal(quin, phan, 0)/self.cte

            plt.plot(self.zvals, self.logatoz(w1/w2))
            plt.title('mquin %f, mphan %f'%(self.mquin, self.mphan))
            plt.show()
        return True


    def Vtotal(self, x, y, select):
        """Cuadratic potential and its derivatives wrt phi or psi"""
        #a=self.delta
        #b=self.beta
        if select == 0:
            #Vtotal = 0.5*(x*self.mquin)**2  +  0.5*(y*self.mphan)**2 + self.beta*(x*y)**2.
            Vtotal = self.mquin**2*(1-np.cos(x*self.delta))  + self.mphan**2*(1-np.cos(y*self.beta))
            #Vtotal = 0.5*(x*self.mquin)**2  +  0.5*(y*self.mphan)**2 - 0.25*self.beta*y**4.
        elif select == 'phi':
            #Vtotal = x*self.mquin**2 + 2.0*self.beta*x*y**2
            Vtotal = self.mquin**2*self.delta*np.sin(x*self.delta)
        elif select == 'psi':
            #Vtotal = y*self.mphan**2 + 2.0*self.beta*y*x**2
            Vtotal = self.mphan**2*self.beta*np.sin(y*self.beta)
        return Vtotal


    def rhode(self, x_vec):
        quin, dotquin, phan, dotphan = x_vec
        Ode = 0.5*dotquin**2 - 0.5*dotphan**2 + self.Vtotal(quin, phan, 0) / self.cte
        return Ode


    def eos(self, x_vec):
        quin, dotquin, phan, dotphan = x_vec
        w1 = 0.5*dotquin**2 - 0.5*dotphan**2 - self.Vtotal(quin, phan, 0)/self.cte
        w2 = 0.5*dotquin**2 - 0.5*dotphan**2 + self.Vtotal(quin, phan, 0)/self.cte
        return w1/w2



    def hubble(self, lna, x_vec=None, SF = True):
        a = np.exp(lna)
        if SF:
            Ode = self.rhode(x_vec)
        else:
            Ode      = 1.0-self.Om

        return self.h*np.sqrt(self.Ocb/a**3 + self.Omrad/a**4 + Ode)


    def logatoz(self, func):
        # change function from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])



    def RHS(self, x_vec, lna):
        sqrt3H0 = np.sqrt(self.cte)
        quin, dotquin, phan, dotphan = x_vec
        hubble  = self.hubble(lna, x_vec)
        return [sqrt3H0*dotquin/hubble, -3*dotquin - self.Vtotal(quin, phan, 'phi')/(sqrt3H0*hubble),
                sqrt3H0*dotphan/hubble, -3*dotphan + self.Vtotal(quin, phan, 'psi')/(sqrt3H0*hubble)]



    def solver(self, quin0, dotquin0, phan_0, dotphan0):
        y0       = [quin0, dotquin0, phan_0, dotphan0]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result



    def calc_Ode(self, mid):
        "Select Quintess, Phantom or Quintom"
        if (self.varymquin) and (self.mphan == 0)   and (self.beta ==0): quin0, phan0 = mid, 0
        elif (self.varymphan) and (self.mquin == 0) and (self.beta ==0): quin0, phan0 = 0, mid
        else : quin0, phan0 = mid, mid*self.iniphi
            #Still figuring out these initial conditions.

        sol = self.solver(quin0 , 0.0, phan0, 0.0).T

        quin, dotq, phan, dotp = sol
        #rho_quin =   0.5*dotq[-1]**2 + self.Vquin(quin[-1], 0)/self.cte
        #rho_phan =  -0.5*dotp[-1]**2 + self.Vphan(phan[-1], 0)/self.cte
        rho =  self.rhode(sol)[-1]
        Ode = rho*(self.h/self.hubble(0.0, [quin[-1], dotq[-1], phan[-1], dotp[-1]]))**2
        #Ode = (rho_quin + rho_phan)*(self.h/self.hubble(0.0, [quin[-1], dotq[-1], phan[-1], dotp[-1]], SF=True))**2
        tol = (1- self.Ocb- self.Omrad) - Ode

        return sol, tol, Ode


    def bisection(self):
        "Search for intial condition of phi such that \O_DE today is 0.7"
        lowphi, highphi = 0, 3 #200
        Ttol            = 1E-2
        mid = (lowphi + highphi)/2.0
        while (highphi - lowphi )/2.0 > Ttol:
            sol, tol_mid, Ode = self.calc_Ode(mid)
            tol_low = self.calc_Ode(lowphi)[1]
            if(np.abs(tol_mid) < Ttol):
                  #print 'reach tolerance',  'phi_0=', mid, 'error=', tol_mid
                  return mid, sol
            elif tol_low*tol_mid<0:
                highphi  = mid
            else:
                lowphi   = mid
            mid = (lowphi + highphi)/2.0

        #print 'No solution found!', mid, self.mquin, self.mphan
        ##Check whether rho is constant or nearby
        grad = np.abs(np.gradient(self.rhode(sol))).max()
        self.mid = mid
        mid  = -1 if grad < 1.0E-2 else 0
        #mid = 0
        return mid, sol


    def search_ini(self):
        mid, sol = self.bisection()
        #self.solution = sol
        #self.mid      = mid
        if (mid == 0):
            print 'mass=',  self.mquin, self.mphan, 'solution not found with', self.mid
            self.hub_SF   = interp1d(self.lna, self.hubble(self.lna, SF=False))
            self.w_eos    = interp1d(self.lna, -1*np.ones(len(self.lna)))
            #self.hub_SF_z = self.logatoz(self.lna, np.ones(len(self.lna)))
        elif (mid == -1):
            #print 'looks like LCDM', self.mquin, self.mphan
            self.hub_SF   = interp1d(self.lna, self.hubble(self.lna, SF=False))
            self.w_eos    = interp1d(self.lna, -1*np.ones(len(self.lna)))
            #self.hub_SF_z = self.logatoz(self.hubble(self.lna, sol, SF=False))
        else:
            #(mid is not 0) and (mid is not -1):
            self.hub_SF   = interp1d(self.lna, self.hubble(self.lna, sol))
            self.w_eos    = interp1d(self.lna, self.eos(sol))
            #self.hub_SF_z = self.logatoz(self.hubble(self.lna, sol))
            #else:
            # self.hub_SF   = interp1d(self.lna, np.ones(len(self.lna)))
        return True




    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        #NuContrib=self.NuDensity.rho(a)/self.h**2
        #rhow= a**(-3*(1.0+self.w0+self.wa))*N.exp(-3*self.wa*(1-a))
        #return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+NuContrib+(1.0-self.Om-self.Ok)*rhow)
        lna = log(a)
        if (1./a-1 < 10):
            hubble = (self.hub_SF(lna)/self.h)**2.
        else:
            hubble = self.hubble(self.lna, SF=False)
        return hubble


    def Hubble_a(self, a):
        lna = log(a)
        hubble = 100*self.hub_SF(lna)
        return hubble

    def w_de(self, a):
        lna = log(a)
        return self.w_eos(lna)

    #def RHSquared_z(self, z):
    #    hubble = (self.hub_SF_z(z)/self.h)**2.
    #    return hubble



