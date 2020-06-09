

from simplemc.cosmo.paramDefs import a_par, b_par

class GenericModel():
    def __init__(self, varya= True, varyb= True):
        """
        This class may be used to add your own model -- cosmology independent --
        For the moment contains a straight line with two parameters: a, b
        Returns
        -------

        Parameters
        ----------
        varya
        varyb

        Returns
        -------

        """
        self.varya = varya
        self.varyb = varyb

        self.a = a_par.value
        self.b = b_par.value


    # my free params (parameters/priors see ParamDefs.py)
    def freeParameters(self):
        l = []
        if (self.varya): l.append(a_par)
        if (self.varyb): l.append(b_par)
        return l


    def printFreeParameters(self):
        print("Free parameters and values currently accepted:")
        self.printParameters(self.freeParameters())


    def printParameters(self, params):
        for p in params:
            print(p.name, '=' , p.value , '+/-' , p.error)


    def updateParams(self, pars):
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
