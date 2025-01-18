from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
from simplemc.cosmo.paramDefs import Ok_par
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import numpy as np



class HolographicCosmology(LCDMCosmology):
    """
        This is Holographic cosmology.
        This class inherits LCDMCosmology class as the rest of the cosmological
        models already included in SimpleMC.

        :param varyc: variable w0 parameter

    """
    def __init__(self, varyOk=False, varyc=True, varyD=False):
        # Holographic parameter
        # Notation: capital C = 3c**2*M_p, here, we use little c
        self.c_par = Parameter("c", 0.7, 0.1, (0.5, 2.0), "c")
        self.D_par = Parameter("D", 0.0, 0.1, (0., 3.0), "D")

        self.varyOk = varyOk
        self.varyc  = varyc
        self.varyD  = varyD

        self.Ok = Ok_par.value
        self.c_hde  = self.c_par.value
        self.D_hde = self.D_par.value

        # This value is quite related to the initial z
        self.zini = 3
        self.zvals = np.linspace(0, self.zini, 100)

        LCDMCosmology.__init__(self, disable_radiation=True)
        self.updateParams([])


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varyOk): l.append(Ok_par)
        if (self.varyc):  l.append(self.c_par)
        if (self.varyD): l.append(self.D_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name  == "c":
                self.c_hde = p.value
            if p.name == "D":
                self.D_hde = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False

        self.initialize()
        return True

    def Q_term(self, z):
        D = self.D_hde
        term1 = 1./(D - 2)
        Q = (2-D)*(self.c_hde)**(2*term1)*(100*self.h*np.sqrt(self.Om*(1 + z)**3))**(-D*term1)
        return Q

    def extra_term(self, z, Ode):
        D = self.D_hde
        term1 = 1./(D - 2)
        return (1- Ode)**(0.5*D*term1)*Ode**(-term1)


    # Right hand side of the equations
    def RHS_hde(self, vals, z):
        Ode = vals
        fact = self.D_hde + 1 + self.Q_term(z)*self.extra_term(z, Ode)
        dOmega = -Ode*(1 - Ode)*fact/(1 + z) + 0
        return dOmega


    def initialize(self):
        """
        Main method.
        """
        Ode0 = (1 - self.Om - self.Ok)

        result_E = odeint(self.RHS_hde, Ode0, self.zvals)
        self.Ode = interp1d(self.zvals, result_E[:, 0])

        return True



    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        z = 1./a-1
        if z>self.zini:
            Ode =  (1-self.Om-self.Ok)
            hubble = self.Omrad/a**4 + self.Ocb/a**3 + self.Ok/a**2 + Ode
        else:
            hubble = self.Om*(1+z)**3/(1-self.Ode(z))

        return hubble


    def EoS(self, z):
        Ode = self.Ode(z)
        w = -(1+self.D_hde)/3. - (1./3)*self.Q_term(z)*self.extra_term(z, Ode)
        return w
