## This is CDM cosmology with w, wa and Ok

import timeit

import numpy as np
import matplotlib.pyplot as plt
from LCDMCosmology import *
#from scipy.interpolate import interp1d
#from scipy.misc import derivative
#from scipy.interpolate import UnivariateSpline

from   scipy      import interpolate, integrate
from   scipy.misc import derivative

class HornFcoCosmology(LCDMCosmology):
	def __init__(self, varyOk=False):
		## two parameters: Om and h

		self.varyOk=varyOk

		vary = False
		self.varyf1 = vary
		self.varyf2 = vary
		self.varyf3 = vary
		self.varyf4 = vary

		self.varyc1 = True
		self.varyc2 = True
		self.varyc3 = True
		self.varyc4 = True


		self.Ok = Ok_par.value
		self.f1 = f1_par.value
		self.f2 = f2_par.value
		self.f3 = f3_par.value
		self.f4 = f4_par.value

		self.c1 = c1_par.value
		self.c2 = c2_par.value
		self.c3 = c3_par.value
		self.c4 = c4_par.value

		LCDMCosmology.__init__(self)

		self.name = 'Alberto'

		self.npoints =  20

		self.zstart  =  0.0
		self.zend    =  3.0
		self.epsilon =  1e-5

		self.zvals   = np.linspace(self.zstart, self.zend, self.npoints)
		#eta from negative to zero, means from z decreasing
		self.etavals = np.linspace(self.eta(self.zstart), self.eta(self.zend), self.npoints )[::-1]

		self.z       = np.array([self.zstart-self.epsilon, 0.8, 1.6, 2.4, self.zend+self.epsilon])
		#eta from negative to zero, means from z decreasing
		self.eta     = np.log( 1.0/(1.0+ self.z))[::-1]

		self.itype   = 'cubic'

		self.updateParams([])

	## my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
	def freeParameters(self):
		l=LCDMCosmology.freeParameters(self)
		if (self.varyf1): l.append(f1_par)
		if (self.varyf2): l.append(f2_par)
		if (self.varyf3): l.append(f3_par)
		if (self.varyf4): l.append(f4_par)
		if (self.varyc1): l.append(c1_par)
		if (self.varyc2): l.append(c2_par)
		if (self.varyc3): l.append(c3_par)
		if (self.varyc4): l.append(c4_par)
		if (self.varyOk): l.append(Ok_par)
		return l


	def updateParams(self,pars):
		ok=LCDMCosmology.updateParams(self,pars)
		if not ok:
			return False
		for p in pars:
			if p.name=="f1":
				self.f1=p.value
			elif p.name=="f2":
				self.f2=p.value
			elif p.name=="f3":
				self.f3=p.value
			elif p.name=="f4":
				self.f4=p.value
			elif p.name=="c1":
				self.c1=p.value
			elif p.name=="c2":
				self.c2=p.value
			elif p.name=="c3":
				self.c3=p.value
			elif p.name=="c4":
				self.c4=p.value
			elif p.name=="Ok":
				self.Ok=p.value
				self.setCurvature(self.Ok)
				if (abs(self.Ok)>1.0):
					return False

		self.ini_f_c()
		self.ini_lambda()
		#plt.plot(self.zvals, self.ffunction(self.etavals))
		#plt.show()
		return True


	def something(self):
		print self.name


	def eta(self, z):
		return np.log( 1.0/(1.0+z) )


	def interp(self, x, y):
		return interpolate.interp1d(x, y, kind=self.itype) #, copy=False)


	def ffunction(self, N):
		# decreasing in z
		x = self.eta
		y = [1.0, self.f1, self.f2, self.f3, self.f4]
		#print 'test', x, y
		g = self.interp(x, y)
		f = g(N)
		return f


	def cfunction(self, N):
		# increasing in z
		x = self.eta
		y = [0.0, self.c1, self.c2, self.c3, self.c4]
		#print 'test', x, y
		c = self.interp(x, y)
		return c(N)


	def ini_f_c(self):
		etavals  = self.etavals

		#print 'eta',etavals
		self.f   = np.array(self.ffunction(etavals))
		self.df  = np.array([self.deriv(self.ffunction, i) for i in etavals])
		self.ddf = np.array([self.dderiv(self.ffunction, i) for i in etavals])

		#print 'f,df,ddf',self.f, self.df, self.ddf
		self.c   = np.array(self.cfunction(etavals)) #self.cfunction(etavals)
		self.dc  = np.array([self.deriv(self.cfunction, i) for i in etavals])

		#print 'c, dc', self.c, self.dc
		self.interc  = self.interp(etavals, self.c)
		self.interf  = self.interp(etavals, self.f)
		self.interdf = self.interp(etavals, self.df)
		return True


	def Aterm(self):
		num = self.df*( self.ddf - 3.0*self.df - 4.0*self.f)
		den = (self.df + self.f)*( self.df + 2.0*self.f)
		return num/den


	def Bterm(self):
		rhom = self.Om*np.exp(self.etavals)**(-3)
		pm   = 0.0
		num1 = self.dc*(  self.f*( 2.*self.f + 3.*self.df ) + self.df**2)
		num2 = self.c*( 4.*self.f*(3.*self.f + 5.*self.df ) + self.df*( self.ddf + 9.*self.df))
		num3 = self.df*( (self.ddf - self.f)*rhom + 3.*( self.df + self.f)*pm)

		num  = num1 + num2 + num3
		den  = (self.df + self.f)*( self.df + 2.*self.f)
		return num/den


	def deriv(self, g, inter):
		return derivative(g, inter, dx=1e-6)


	def dderiv(self, g, inter):
		return derivative(g, inter, n=2, dx=1e-6)


	def my_Int(self, h, N):
		return integrate.quad(lambda x: h(x), 0, N)[0]


	def first_term(self):
		Atmp = self.interp(self.etavals, self.Aterm())
		II   = np.array([self.my_Int(Atmp, i) for i in self.etavals])
		return II


	def second_term(self, first):
		etavals   = self.etavals

		Btmp      = self.interp(etavals, self.Bterm())
		argtmp    = self.interp(etavals, np.exp(first*Btmp(etavals)))
		II_second = np.array([self.my_Int(argtmp, i) for i in etavals])
		return II_second


	def ini_lambda(self):
		first  = self.first_term()
		lam    = np.exp(-first)*((1-self.Om) - self.second_term(first))
		self.lambda_term = self.interp(self.etavals, lam)
		#print '**', lam
		return True


	def hub_horn(self, N):
		rhom = self.Om*(np.exp(N)**(-3))
		num  = rhom + self.interc(N) + self.lambda_term(N)
		den =  self.interf(N) + self.interdf(N)
		return 100*self.h*np.sqrt(num/den)


	def hubble(self, N):
		rhom = self.Om*(np.exp(N)**(-3))
		num = rhom + self.c + self.lambda_term(N)
		den = 1.0*( self.f + self.df )
		return 100*self.h*np.sqrt(num/den)[::-1]


		# noinspection PyUnreachableCode
	def RHSquared_a(self, a):
		z = 1./a - 1.0
		N = log(a)

		if (z>= 0.0 and z<= self.zend):
			H = self. hub_horn(N)
		else:
			H = (self.Ocb/a**3 + self.Ok/a**2 + self.Omrad/a**4 + (1-self.Om))
		return H


##=================Plotting ==============

	def plot_func(self):
		self.ini_f_c()
		plt.plot(self.etavals, self.c)
		#plt.plot(self.etavals, self.df)
		#plt.plot(self.etavals, self.ddf)
		plt.show()



	def plot_A(self):
		self.ini_f_c()
		#plt.plot(self.etavals, map(self.Aterm, self.etavals), 'ro')
		plt.plot(self.etavals, self.first_term(), 'ro')
		#plt.plot(self.etavals, self.second_term())
		#plt.plot(self.etavals, self.second_term(), 'ro')
		plt.show()

	def plot_lam(self):
		self.ini_f_c()
		lt = self.lambda_term(self.etavals)
		print lt
		plt.plot(self.etavals, lt)
		plt.show()


	def plot_hubble(self):
		self.ini_f_c()
		fig = plt.figure(figsize=(9,7))
		ax1 = fig.add_subplot(111)
		plt.plot(self.zvals, self.hubble(self.etavals))
		dataHz = np.loadtxt('data/Hz_all.dat')
		redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
		ax1.errorbar(redshifts, obs, errors, xerr=None,
				color='purple', marker='o', ls='None',
				elinewidth =2, capsize=3, capthick = 1)
		ax1.legend(loc='lower right', frameon=False)
		plt.ylabel('$H(z)$', fontsize=20)
		plt.xlabel('redshift $z$', fontsize=20)
		plt.show()
