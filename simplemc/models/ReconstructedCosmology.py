
#TODO: incorporate linear and spline interpolation
#TODO: steps, bins or tanh function

from scipy.interpolate import InterpolatedUnivariateSpline
from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
from scipy.interpolate import interp1d
from scipy.integrate import quad
import numpy as np

class ReconstructedCosmology(LCDMCosmology):
    def __init__(self,
                 select='rho_de',   # select which function wants to reconstruct i.e. [hubble, rho_de, w_de]
                 kspline= 1,            # k=0 tanh, =1 linear, >=2 splines
                 nodes= 5,              # numer of nodes used in the interpolation
                 eta= 100               # smoothing params for tannh
                 ):
        """
        This file implements several types of reconstruction techniques
        for different functions, i.e. H(z), rho(z), w(z), Q(z)
        """

        self.select = select
        self.kspline = kspline
        self.nnodes = nodes
        self.eta = eta

        # zend is the last point in linear-interpolation
        # however is the first point of the last step in the binning/tanh version
        self.zini = 0.0
        self.zend = 2.0

        # range to perform the interpolation
        self.zvals = np.linspace(self.zini, 3.0, 50)

        #this is useful when using MCMC but not for nested
        mean = -1 if self.select == 'w_de' else 1
        priors = (-2, 0) if self.select == 'w_de' else (0.2, 1.5)
        sigma = 0.2

        self.pname = 'amp_'
        self.params = [Parameter(self.pname+str(i), mean, sigma, priors, self.pname+str(i)) for i in range(nodes)]
        self.pvals = [i.value for i in self.params]

        fixOm = True if self.select == 'hubble' else False
        LCDMCosmology.__init__(self, mnu=0, fixOm=fixOm)
        self.updateParams([])


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        l+= self.params
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for i, p in enumerate(pars):
            if p.name == (self.pname+str(i)):
                self.pvals[i] = p.value

        self.initialize()
        return True


    def bin_tanh(self, z, zi, fi, fi1):
        return 0.5*(fi1-fi)*( 1+ np.tanh( (z-zi)*self.eta ) )


    def fun_tanh(self, z):
        f = self.y[0]
        for i in range(self.mnodes):
            f += self.bin_tanh(z, self.x[i+1], self.y[i], self.y[i+1])
        return f


    def initialize(self):


        if self.select in ('hubble', 'rho_de'):
            # the first node-amplitude is fixed by today's conditions
            self.y = [1.0] + self.pvals
            self.mnodes = self.nnodes
        else:
            self.y = self.pvals
            self.mnodes = self.nnodes-1

        delta = (self.zend - self.zini)/(self.mnodes)
        self.x = [self.zini + delta*i for i in range(self.mnodes+1)]

        if self.kspline == 0:
            ftanh = [self.fun_tanh(i) for i in self.zvals]
            self.spline = interp1d(self.zvals, ftanh)
        else:
            self.spline = InterpolatedUnivariateSpline(self.x, self.y, k=self.kspline)

        if self.select== 'w_de':
            expr = lambda x: (1+self.spline(x))/(1+x)
            integ = lambda x: quad(expr, 0, x)[0]
            rho = np.exp(3 * np.array(list(map(integ, self.zvals))))
            self.rhow = interp1d(self.zvals, rho)

        return True



    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        z= 1./a - 1
        NuContrib = 0 #self.NuDensity.rho(a)/self.h**2
        if z< 3.0:
            if self.select == 'hubble':
                return (self.spline(z))**2
            else:
                if self.select == 'rho_de':
                    return self.Ocb/a**3 + self.Omrad/a**4 + NuContrib + (1.0-self.Om)*self.spline(z)
                elif self.select == 'w_de':
                    return self.Ocb/a**3 + self.Omrad/a**4 + NuContrib + (1.0-self.Om)*self.rhow(z)
        else:
            return self.Ocb/a**3 + self.Omrad/a**4 + NuContrib + (1.0-self.Om)

        return True


    def rho_de(self, z):
        if self.select == 'w_de':
            return self.rhow(z)
        else:
            return self.spline(z)