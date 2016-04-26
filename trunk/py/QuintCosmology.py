## This is a CDM cosmology with \phi
from pylab import *
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from LCDMCosmology import *

class QuintCosmology(LCDMCosmology):
    def __init__(self, varylam=True, varyV0=False, varyA=True, varylB=True):
        ## two parameters: Om and h

	self.varyV0= varyV0
	self.varylam= varylam
	self.varyA= varyA
	self.varylB= varylB
	
	self.lam= lam_par.value
	self.lB = lB_par.value
	self.A	= A_par.value
	self.V0=  V0_par.value

#	self.V0i= 1.0*exp(self.lB)/self.A

#	self.oC=LCDMCosmology()
        LCDMCosmology.__init__(self,mnu=0)

	self.lna = linspace(-35, 10, 300)
	self.Cte = sqrt(3.0)*self.h

	## force caching
	self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
	l=LCDMCosmology.freeParameters(self)
	if (self.varylam): l.append(lam_par)
	if (self.varyV0): l.append(V0_par)
	if (self.varyA): l.append(A_par)
	if (self.varylB): l.append(lB_par)
	return l

    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
#	self.oC.updateParams(pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="lam":
                self.lam=p.value
	for p in pars:
	    if p.name=="A":
		self.A=p.value
	for p in pars:
	    if p.name=="lB":
		self.lB=p.value	
	for p in pars:
	    if p.name=="V0":
		self.V0=p.value	

	self.Ini_phi()
        return True


    def Pot(self, x, i):
	B= self.lB	
        AA= self.A/self.lam**2
	funct1 = (x-B)**2 + AA
	funct2 = 2.0*(x-B) - self.lam*funct1
	if i==0:
	    return funct1*exp(-self.lam*(x-self.lB))/AA*self.vpo
#           return funct1*self.V0i*exp(-self.lam*x+self.vpo)
        if i==1:
	    return funct2*exp(-self.lam*(x-self.lB))/AA*self.vpo
#           return funct2*self.V0i*exp(-self.lam*x+self.vpo)
        else:
           print 'wrong choice'
           stop()

				# x ->\phi, y->\dotphi    
    def RHS(self, x_vec, lna):
        x, y = x_vec
        return [sqrt(3.0)*y/self.hub(lna,x_vec), -3*y -self.Pot(x,1)/self.h/(self.Cte*self.hub(lna,x_vec))]


    def solver(self, x0, vpo):
	self.vpo= vpo
	self.x0 = x0
        y0 = [self.x0, 0]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result


    def Ini_phi(self):
	lowr, highr = 0, 2
	tol, tol1 =101, 100
	Ttol= 5e-3
	count=0
 					#search initial conditions
	if True:
	   while (abs(tol)>Ttol):
		mid= (lowr+highr)/2.0
		sol = self.solver(0.0, mid)
		
		Omegal= (0.5*sol[-1,1]**2+self.Pot(sol[-1,0],0)/self.Cte**2)/self.hub(0.0, sol[-1])**2
		tol = (1.0-self.Om) - Omegal

		if(abs(tol) < Ttol):
		   #print 'reach tolerance', abs(tol), count
		   break
		else:
		   if(tol<0):
 			highr = mid
		   else:
			lowr = mid
	 	   						
		count+=1
		if (count >10):
		   print 'No solution found!'
		   break    

	#print 'mid', self.lB, self.lam, self.A
	#sol =self.solver(0., self.V0)
 	#print 'O_L calc', (0.5*sol[-1,1]**2+self.Pot(sol[-1,0],0)/self.Cte**2)/self.hub(0.0, sol[-1])**2
	self.sol =sol
	xx, yy = sol.T
	self.Ophi = interp1d(self.lna, (0.5*yy**2+self.Pot(xx,0)/self.Cte**2))
	self.hubble=interp1d(self.lna, (self.hub(self.lna, sol.T))**2)	
	return self.sol	
	
				
    def hub(self, lna, x_vec):
	a=exp(lna)
        x, y  = x_vec
        NuContrib=0 #self.NuDensity.rho(a)/self.h**2
        return sqrt(0.5*y**2 + self.Pot(x,0)/self.Cte**2 + self.Ocb/a**3 + self.Omrad/a**4 +NuContrib)
 
 
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
	lna=log(a)
	return self.hubble(lna)

    def RHSquared_lna(self,lna):
        return self.hubble(lna)
 
    def O_phi(self, lna):
	return self.Ophi(lna)*(1.0/self.hub((lna), self.sol.T))**2
    

#    def prefactor(self):
        #if self.userd_DE:
#        self.rd = self.oC.rd*((1.-self.O_phi(-7.0)[-1])**(0.5))
        #else:
        #self.rd = self.oC.rd
	#print 'rd ', self.c_/(self.rd*self.h*100)
#        return  self.c_/(self.rd*self.h*100)
	
