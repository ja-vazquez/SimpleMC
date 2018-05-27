## This is CDM cosmology with w, wa and Ok


import math as N
import numpy as np
import matplotlib.pyplot as plt
from LCDMCosmology import *

class QuintomCosmology(LCDMCosmology):
    def __init__(self, varymquin=False, varymphan=False, varyiniphi=False):
        ## two parameters: Om and h

	"""Is better to start the chains at masses equal one, othewise
		may take much longer"""

	self.varymquin = varymquin
        self.varymphan = varymphan
	self.varyiniphi= varyiniphi

        self.mquin    = mquin_par.value   
        self.mphan    = mphan_par.value
	self.iniphi   = iniphi_par.value

	LCDMCosmology.__init__(self, mnu=0)

	self.lna   = linspace(-10, 0, 200)
	self.z     = np.exp(-self.lna) - 1.
	self.zvals = np.linspace(0, 3, 100)
	
	self.cte   = 3.0*self.h**2
	self.comments = False

	self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l=LCDMCosmology.freeParameters(self)
        if (self.varymquin) : l.append(mquin_par)
        if (self.varymphan) : l.append(mphan_par)
	if (self.varyiniphi) : l.append(iniphi_par)
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

	self.search_ini()

	if  False:
		"Useful for testing purposes"
		quin, dquin, phan, dphan = self.solution
		
		w1 = 0.5*dquin**2 - self.Vquin(quin, 0)/self.cte - 0.5*dphan**2  - self.Vphan(phan, 0)/self.cte
		w2 = 0.5*dquin**2 + self.Vquin(quin, 0)/self.cte - 0.5*dphan**2  + self.Vphan(phan, 0)/self.cte
		
		plt.plot(self.zvals, self.logatoz(w1/w2))
		plt.title('mquin %f, mphan %f'%(self.mquin, self.mphan))
		plt.show()
        return True


    def Vphan(self, x, select):
	#0-Potential, 1-Derivative
	mphan = self.mphan
        if select == 0:     return 0.5*(x*mphan)**2
        if select == 1:     return x*mphan**2


    def Vquin(self, x, select):
        #0-Potential, 1-Derivative
	mquin = self.mquin
        if select == 0:     return 0.5*(x*mquin)**2
        if select == 1:     return x*mquin**2


    def hubble(self, lna, x_vec=None, SF = False):
        a = np.exp(lna)
        if SF:
            quin, dotquin, phan, dotphan = x_vec
            Ode_quin =  0.5*dotquin**2 + self.Vquin(quin, 0)/self.cte
            Ode_phan = -0.5*dotphan**2 + self.Vphan(phan, 0)/self.cte
            Ode      = Ode_quin + Ode_phan
        else:
            Ode      = 1.0-self.Om 

        return self.h*np.sqrt(self.Ocb/a**3 + self.Omrad/a**4 + Ode)


    def logatoz(self, func):
        # change function from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return  np.interp(self.zvals, self.z[::-1], functmp[::-1])


    def RHS(self, x_vec, lna):
	quin, dotquin, phan, dotphan = x_vec
	sqrt3H0 = np.sqrt(self.cte)
	hubble  = self.hubble(lna, x_vec, SF=True)
        return [sqrt3H0*dotquin/hubble, -3*dotquin - self.Vquin(quin, 1)/(sqrt3H0*hubble),
                sqrt3H0*dotphan/hubble, -3*dotphan + self.Vphan(phan, 1)/(sqrt3H0*hubble)]



    def solver(self, quin0, dotquin0, phan_0, dotphan0):
        y0       = [10**quin0, dotquin0, phan_0, dotphan0]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result



    def calc_Ode(self, mid):
	"Select Quintess, Phantom or Quintom"
	if   self.varymquin and (not self.varyiniphi):   a, b = mid, 0
	elif self.varymphan and (not self.varyiniphi):   a, b = 0, mid
	elif self.varyiniphi:                            a, b = self.iniphi, mid 
	#np.sqrt(2*3*self.h**2*(1-self.Om) + (self.mphan*mid)**2)/self.mquin, mid

	sol = self.solver(a , 0.0, b, 0.0).T

        quin, dotq, phan, dotp = sol
        rho_quin =   0.5*dotq[-1]**2 + self.Vquin(quin[-1], 0)/self.cte
        rho_phan =  -0.5*dotp[-1]**2 + self.Vphan(phan[-1], 0)/self.cte

        Ode = (rho_quin + rho_phan)*(self.h/self.hubble(0.0, [quin[-1], dotq[-1], phan[-1], dotp[-1]], SF=True))**2
        tol = (1- self.Ocb- self.Omrad) - Ode

        return tol


    def bisection(self):
         "Search for intial condition of phi such that \O_DE today is 0.7"
         lowphi, highphi = -1, 4
         Ttol            = 1E-3
         mid = (lowphi + highphi)/2.0
         while (highphi - lowphi )/2.0 > Ttol*0.1:
	       ode_mid = self.calc_Ode(mid)
               ode_low = self.calc_Ode(lowphi)
               if(np.abs(ode_mid) < Ttol):
                  #print 'reach tolerance',  'phi_0=', mid, 'error=', ode_mid
                  return mid
               elif ode_low*ode_mid<0:
		  highphi  = mid
	       else:
		  lowphi   = mid
               #print mid, self.calc_Ode(mid) ,(1- self.Ocb- self.Omrad)
               mid = (lowphi + highphi)/2.0

	 print 'No solution found!', mid, self.mquin, self.mphan
	 mid = 0
	 return mid


    def search_ini(self):
        mid = self.bisection()
        if   self.varymquin and (not self.varyiniphi):   a, b = mid, 0
	elif self.varymphan and (not self.varyiniphi):   a, b = 0, mid
	elif self.varyiniphi:                            a, b = self.iniphi, mid

	#if   self.mquin == 0:   a, b = 0, mid
	#elif self.mphan == 0:   a, b = mid, 0
	#else: a,b = self.iniphi, mid #np.sqrt(2*3*self.h**2*(1-self.Om) + (self.mphan*mid)**2)/self.mquin, mid
	
        sol = self.solver(a , 0, b, 0).T
        self.hub_SF   = interp1d(self.lna, self.hubble(self.lna, sol, SF=True))
	self.hub_SF_z = self.logatoz(self.hubble(self.lna, sol, SF=True))
        self.solution = sol
        return mid




    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        #NuContrib=self.NuDensity.rho(a)/self.h**2
        #rhow= a**(-3*(1.0+self.w0+self.wa))*N.exp(-3*self.wa*(1-a))
        #return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+NuContrib+(1.0-self.Om-self.Ok)*rhow)
	lna = log(a)
	hubble = (self.hub_SF(lna)/self.h)**2.
	return hubble


    def RHSquared_z(self, z):
        hubble = self.hub_test(z)
        return hubble
    


