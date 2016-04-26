## This is a CDM cosmology with \phi

from LCDMCosmology import *

class QuintCosmology2(LCDMCosmology):
    def __init__(self, varylam=True):
        ## two parameters: Om and h
	self.varylam=varylam
	self.lam=lam_par.value
        LCDMCosmology.__init__(self,mnu=0)

	## force caching
	self.updateParams([])

    ## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
	l=LCDMCosmology.freeParameters(self)
	if (self.varylam): l.append(lam_par)
	return l

    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        for p in pars:
            if p.name=="lam":
                self.lam=p.value
	self.Solve_Quint()
        return True

			# x ->\dotphi, y->V(\phi), z->DM, h-> Hubble
    def RHS(self, x_vec, t):
        x, y, z, h  = x_vec
	lam= 10**(self.lam)
        func = 1.5*(2*x**2 + z**2+ 1.33*(1 -x**2 -y**2 -z**2))
        return [-3*x+ lam*sqrt(1.5)*y**2+ x*func, -lam*sqrt(1.5)*x*y + y*func, func*z-1.5*z, -h*func]


    def Solve_Quint(self):
	ini = [0.5, sqrt(1.0-self.Om), sqrt(self.Ocb), self.h]
	self.logar= linspace(0.0,-7.1,100)
	self.ilogar= self.logar[::-1]

	y_result = odeint(self.RHS, ini, self.logar)
	xx, yy, zz, hh= y_result.T

	self.hubble=interp1d(self.ilogar, hh[::-1])


    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
	lna=log(a)
	#NuContrib=self.NuDensity.rho(a)/self.h**2
	#return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om))
	return (self.hubble(lna)/self.h)**2
	

