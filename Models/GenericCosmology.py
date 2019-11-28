#from __future__ import print_function
from ParamDefs import a_par, b_par
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




class GenericCosmology(BaseGenericModel):
    def __init__(self, varya= True, varyb= True):
        ## Example used to fit a straigth line two parameters: a, b
        self.varya = varya
        self.varyb = varyb

        self.a = a_par.value
        self.b = b_par.value
        BaseGenericModel.__init__(self)


    # my free params (parameters/priors see ParamDefs.py)
    def freeParameters(self):
        l = BaseGenericModel.freeParameters(self)
        if (self.varya): l.append(a_par)
        if (self.varyb): l.append(b_par)
        return l


    def updateParams(self, pars):
        #BaseGenericModel.updateParams(self,pars)
        for p in pars:
            if p.name == "a":
                self.a = p.value
            elif p.name == "b":
                self.b = p.value
        return True


    def genericModel(self,x):
        return (self.a*x + self.b)


    def prior_loglike(self):
        return 0
