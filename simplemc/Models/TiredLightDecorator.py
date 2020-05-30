##
# This helper class adds a "tired light" parameter
# to any class
##
from ParamDefs import beta_par
import math as N


def TiredLightDecorator(LikeInstance):
    class TiredLightLikelihood(LikeInstance.__class__):
        def __init__(self, LikeInstance):
            # copy  instance into my class
            self.__dict__.update(LikeInstance.__dict__)
            # save type since we need to call parent type
            self.LikeType = LikeInstance.__class__
            self.beta = beta_par.value


        def freeParameters(self):
            return self.LikeType.freeParameters(self) + [beta_par]


        def updateParams(self, pars):
            for p in pars:
                if p.name == "beta": self.beta = p.value
            self.LikeType.updateParams(self, pars)


       # distance modulus
        def distance_modulus(self, z):
            assert(not self.varyPrefactor)
            return 5*N.log10(self.Da_z(z)*(1+z)**(1+self.beta))

    return TiredLightLikelihood(LikeInstance)
