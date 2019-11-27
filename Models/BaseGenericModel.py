from __future__ import print_function
from scipy import *
from Parameter import *
from GenericParamDefs import *

class BaseGenericModel:
    # speed of light

    def __init__(self):
        #self.a=0
        #self.b=149.50
        BaseGenericModel.updateParams(self,[])

    def freeParameters(self):
        ##a_par.setValue(self.a)
        l=[a_par]
        return l
        
    def printFreeParameters(self):
        print("Free parameters and values currently accepted:")
        self.printParameters(self.freeParameters())

    def printParameters(self,params):
        for p in params:
            print(p.name,'=',p.value,'+/-',p.error)

    def updateParams(self,pars):       
        return True
        
    def prior_loglike(self):
        return 0 

    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def genericModel(self,a):
        print("You should not instatiate BaseGenericModel")
        error("BAD")
        


class GenericModel(BaseGenericModel):

    def __init__(self,a=a_par.value, b=b_par.value):
        ## two parameters: a, b
        self.a=a
        self.b=b
        BaseGenericModel.__init__(self)
        
        # Force rd update
        GenericModel.updateParams(self,[])
        
    # to change parameters/priors see GenericParamDefs.py    
    def freeParameters(self):
        l=[]
        l+=[b_par]  
        l+=BaseGenericModel.freeParameters(self)
        return l

    def updateParams(self,pars):
        BaseGenericModel.updateParams(self,pars)        
        for p in pars:
            if p.name=="a":                
                self.a=p.value
            elif p.name=="b":
                self.b=p.value
        return True

    def genericModel(self,x):
        return (self.a*x+self.b)

    def prior_loglike(self):
        return 0
