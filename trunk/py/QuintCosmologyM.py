## This is a CDM cosmology with \phi
from pylab import *
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from LCDMCosmology import *

class QuintCosmology2(LCDMCosmology):
    def __init__(self, varylam=True, varyV0=True):
        ## two parameters: Om and h

	self.varyV0= varyV0
	self.varylam= varylam

	self.V0= V0_par.value
	self.lam= lam_par.value

        LCDMCosmology.__init__(self,mnu=0)

	self.lna = linspace(-9, 0, 500)
	
	## force caching
	self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
	l=LCDMCosmology.freeParameters(self)
	if (self.varylam): l.append(lam_par)
	if (self.varyV0): l.append(V0_par)
	return l

    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="lam":
                self.lam=p.value
	for p in pars:
	    if p.name=="V0":
		self.V0=p.value	

	self.Ini_phi()
        return True


    def Pot(self, x, i):
	if i==0:
           return self.V0*exp(-self.lam*x)/self.h**2
        if i==1:
           return -self.lam*self.V0*exp(-self.lam*x)/self.h**2
        else:
           print 'wrong choice'
           stop()

				# x ->\phi, y->\dotphi    
    def RHS(self, x_vec, lna):
        a=exp(lna)
        x, y = x_vec
        return [y/self.hub(a,x_vec), -3*y -self.Pot(x,1)/self.hub(a,x_vec)]


    def solver(self, x0):
        y0 = [x0, 0]
        y_result = odeint(self.RHS, y0, self.lna)
        return y_result


    def Ini_phi(self):
	lowr, highr = -2.e1, 1e2
	tol, tol1 =100, 100
	Ttol= 1e-3
	count=0
				#search initial conditions
	while (abs(tol)>Ttol):
		mid= (lowr+highr)/2.0
		sol = self.solver(mid)
		sol1 = self.solver(lowr)
		
		Omegal= (0.5*sol[-1,1]**2+self.Pot(sol[-1,0],0))*self.h**2/self.hub(1.0, sol[-1])**2
		Omegal1= (0.5*sol1[-1,1]**2+self.Pot(sol1[-1,0],0))*self.h**2/self.hub(1.0, sol1[-1])**2
		tol = (1.0-self.Om) - Omegal
		tol1 = (1.0-self.Om) - Omegal1
		#print 'Omega_l',Omegal, '\phi_ini= ', mid

		if(abs(tol) < Ttol):
		   print 'reach tolerance', abs(tol)
		   break
		else:
		   if(tol*tol1>0):
 			lowr = mid
		   else:
			highr = mid
	 	   						
		count+=1
		if (count >50):
		   'No solution found!'
		   break     
	self.hubble=interp1d(self.lna, (self.hub(exp(self.lna), sol.T)/self.h)**2)	

				
    def hub(self, a, x_vec):
        x, y  = x_vec
        NuContrib=self.NuDensity.rho(a)/self.h**2
        return self.h*sqrt(0.5*y**2 + self.Pot(x,0) + self.Ocb/a**3 + self.Omrad/a**4 +NuContrib)

    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
	lna=log(a)
	return self.hubble(lna)
	

