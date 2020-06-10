
from simplemc.cosmo.paramDefs import Anfw_par, rs_par
import numpy as np


class RotationCurves():
    def __init__(self, varya= True, varyb= True):
        """
        Class to constrain rotational curves profiles,
            Here we assume a NFW profile
        Parameters
        ----------
        varya
        varyb

        Returns
        -------

        """
        ## Example used to fit a straigth line two parameters: a, b
        self.varya = varya
        self.varyb = varyb

        self.Anfw = Anfw_par.value
        self.rs = rs_par.value


    # my free params (parameters/priors see ParamDefs.py)
    def freeParameters(self):
        l = []
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
