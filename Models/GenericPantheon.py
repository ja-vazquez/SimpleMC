
#from __future__ import print_function
from ParamDefs import *
import sys
import numpy as np
sys.path = ["../Cosmo"] + sys.path
import matplotlib.pyplot as plt

class BaseGenericPModel:
    def __init__(self):
        #This class may be used to add your own model
        BaseGenericPModel.updateParams(self,[])

    def freeParameters(self):
        l = []
        return l

    def printFreeParameters(self):
        print("Free parameters and values currently accepted:")
        self.printParameters(self.freeParameters())

    def printParameters(self, params):
        for p in params:
            print(p.name, '=' , p.value , '+/-' , p.error)

    def updateParams(self, pars):
        return True


    def prior_loglike(self):
        return 0 

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def genericModel(self,a):
        print("You should not instatiate BaseGenericModel")
        sys.exit(1)


class GenericPantheon(BaseGenericPModel):
    def __init__(self):
        self.zini    = 0 # 13.9
        self.parvals = zbin_par

        self.zs = [i.value for i in self.parvals]
        names   = [i.name  for i in self.parvals]
        self.index = dict((j,i) for i,j in enumerate(names))

        self.zvals = np.logspace(np.log10(0.01),np.log10(2.261), len(self.parvals)+1)

        BaseGenericPModel.__init__(self)



    # my free params (parameters/priors see ParamDefs.py)
    def freeParameters(self):
        l = BaseGenericPModel.freeParameters(self)
        l+= self.parvals
        return l


    def updateParams(self, pars):
        #BaseGenericModel.updateParams(self,pars)
        for p in pars:
            i = self.index[p.name]
            self.zs[i] = p.value

        return True


    def mu_bar(self,n,z,m1,m2):
        alpha = np.log10(z/self.zvals[n])/np.log10(self.zvals[n+1]/self.zvals[n])
        return (1- alpha)*m1 + alpha*m2


    def genericPModel(self, z):
        if z>=self.zvals[0] and z< self.zvals[len(self.parvals)]:
            for i in range(len(self.parvals)):
                if z>=self.zvals[i] and z< self.zvals[i+1]:
                    if i ==0:
                        y = self.mu_bar(0,  z, self.zini, self.zs[0])
                    else:
                        y = self.mu_bar(i,  z, self.zs[i-1], self.zs[i])
        else:
            y = 0
        return y

    def prior_loglike(self):
        return 0
