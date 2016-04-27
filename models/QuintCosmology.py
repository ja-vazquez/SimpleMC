## This is a CDM cosmology with \phi
from pylab import *
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from LCDMCosmology import *
import scipy.optimize as optimize

class QuintCosmology(LCDMCosmology):
    def __init__(self, varylam=True, varyV0=True, varyA=True, varylB=True):
        ## two parameters: Om and h

	self.varyV0= varyV0
	self.varylam= varylam
	self.varyA= varyA
	self.varylB= varylB
	
	self.lam= lam_par.value
	self.lB = lB_par.value
	self.A	= A_par.value
	self.V0=  V0_par.value

	self.V0i= 1.0 #E25 #*exp(self.lB)/self.A

        LCDMCosmology.__init__(self,mnu=0)

	self.lna = linspace(-20, 5, 300)

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

	self.Ini_phi(num=0)
        return True


    def Pot(self, x, i):
	B  = self.lB	
        AA = self.A/self.lam**2
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
        return [sqrt(3.0)*y/self.hub(lna,x_vec), -3*y - sqrt(3.)*self.Pot(x,1)/self.hub(lna,x_vec)]


    def phidot(self, x0):
	return sqrt(4*self.Pot(x0, 0))


    def solver(self, x0, vpo):
	self.vpo= vpo
	self.x0 = x0
        y0 = [self.x0, self.phidot(self.x0)]
        y_result = odeint(self.RHS, y0, self.lna, h0=1E-10)
        return y_result


    def find_phi0(self, xx):
	#print self.Omrad
	return (self.Pot(xx,1)**2/self.Pot(xx,0)/4. - self.Pot(xx,0))*exp(4*self.lna[0])*(3./self.Omrad) - 1.
	#return (self.Pot(xx,1)**2/self.Pot(xx,0)/4. - 3.*self.Pot(xx,0))*exp(4*self.lna[0])*(3./self.Omrad) - 1.

    def Ini_phi(self, num=1):
	lowr, highr = 0, 2
	tol, tol1   = 101, 100
	Ttol  = 5e-3
	count = 0
 					#search initial conditions
	if True:
	   phi_0 = optimize.bisect(self.find_phi0, -20, 20, xtol= 0.001) if num != 0 else 0
	   print 'phi_0 = ', phi_0
	   while (abs(tol)>Ttol):
		mid= (lowr+highr)/2.0
		sol = self.solver(phi_0, mid)
		
		Omegal= (0.5*sol[-1,1]**2 + self.Pot(sol[-1,0],0))/self.hub(0.0, sol[-1])**2
		tol = (1.0-self.Om) - Omegal

		if(abs(tol) < Ttol):
		   print 'reach tolerance', abs(tol), count
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

#	sol =self.solver(-4.6, self.V0)
#	if num !=0 : print 'solll', optimize.bisect(self.find_phi0, -10, 0)
	print '***',self.V0
	self.sol =sol
	xx, yy   = sol.T
	self.Ophi  = interp1d( self.lna, 0.5*yy**2+self.Pot(xx,0) )
	self.hubble= interp1d( self.lna, self.hub(self.lna, sol.T)**2 )	
	return self.sol	
	
				
    def hub(self, lna, x_vec):
	a=exp(lna)
        x, y  = x_vec
        NuContrib = self.NuDensity.rho(a)/self.h**2
        return sqrt(0.5*y**2 + self.Pot(x,0) + self.Ocb/a**3 + self.Omrad/a**4)
 
 
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
	lna=log(a)
	return self.hubble(lna)


    def RHSquared_lna(self,lna):
        return self.hubble(lna)
 

    def O_phi(self, lna):
	return self.Ophi(lna)/self.hub((lna), self.sol.T)**2
   
 
    def w_ede(self, lna):
	xx, yy = self.sol.T
	return (0.5*yy**2 - self.Pot(xx,0))/(0.5*yy**2+self.Pot(xx,0))

	
