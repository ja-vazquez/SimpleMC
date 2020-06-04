#from __future__ import print_function
from simplemc.cosmo.paramDefs import Anfw_par, rs_par
import numpy as np
import sys


class BaseGenericModel:
    def __init__(self):
        #This class may be used to add your own model
        BaseGenericModel.updateParams(self,[])

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




class RotationCurves(BaseGenericModel):
    def __init__(self, varya= True, varyb= True):
        ## Example used to fit a straigth line two parameters: a, b
        self.varya = varya
        self.varyb = varyb

        self.Anfw = Anfw_par.value
        self.rs = rs_par.value
        BaseGenericModel.__init__(self)


    # my free params (parameters/priors see ParamDefs.py)
    def freeParameters(self):
        l = BaseGenericModel.freeParameters(self)
        if (self.varya): l.append(Anfw_par)
        if (self.varyb): l.append(rs_par)
        return l


    def updateParams(self, pars):
        for p in pars:
            if p.name == "Anfw":
                self.Anfw = p.value
            elif p.name == "rs":
                self.rs = p.value
        return True


    def rotation(self,x):
        #def nfw(r, rs =1, A=0.05):
        A  = self.Anfw
        rs = self.rs
        return (A**2*(rs**3)/x)*(np.log((rs+x)/rs) - x/(rs+x))


    def prior_loglike(self):
        return 0
