
#TODO: incorporate linear and spline interpolation
#TODO: steps, bins or tanh function

from scipy.interpolate import InterpolatedUnivariateSpline
from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
import numpy as np

class ReconstructedCosmology(LCDMCosmology):
    def __init__(self, reconstuct='hub', kspline=1, nodes=7, fixOm=True):
        """
        This file implements serveral types of reconstruction techniques
        for different functions, i.e. H(z), rho(z), w(z)
        """

        # select which function wants to reconstruct
        self.select = 'hubble'

        self.kspline = kspline

        self.nnodes = nodes
        self.zini = 0.0
        self.zend = 3.0

        #this is useful when using MCMC
        mean = 1
        sigma = 0.2
        priors = (0.5, 15)

        self.pname = 'amp_'

        self.params = [Parameter(self.pname+str(i), mean, sigma, priors, self.pname+str(i)) for i in range(self.nnodes)]
        self.pvals = [i.value for i in self.params]
        self.zvals = np.linspace(self.zini, self.zend, self.nnodes+1)

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
        for p in pars:
            for i in range(self.nnodes):
                if p.name == (self.pname+str(i)):
                    self.pvals[i] = p.value


        self.initialize()
        return True


    def initialize(self):
        x = x = [self.zini + (self.zend - self.zini) / (self.nnodes - 1) * i
             for i in range(self.nnodes)]
        y = self.pvals

        k = min(3, len(x) - 1)
        self.spline = InterpolatedUnivariateSpline(x, y, k=self.kspline)
        return True


    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    # return (self.Ocb/a**3+self.Omrad/a**4 +(1.0-self.Om)*f(z)  )
    def RHSquared_a(self,a):
        z= 1./a - 1
        if self.select == 'hubble':
            return (self.spline(z)/self.h)**2

        else:
            pass
