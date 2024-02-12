from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter, Ok_par
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from scipy import optimize
import numpy as np



class HolographicCosmology(LCDMCosmology):
    """
        This is Holographic cosmology.
        This class inherits LCDMCosmology class as the rest of the cosmological
        models already included in SimpleMC.

        :param varyc: variable w0 parameter

    """


    def __init__(self, varyOk=False, varyc=True):
        # Holographic parameter
        self.c_par = Parameter("c", 0.7, 0.05, (0.5, 1.), "c")

        self.varyOk = varyOk
        self.varyc  = varyc

        self.Ok = Ok_par.value
        self.c  = self.c_par.value

        # This value is quite related to the initial z
        self.zini = 10
        self.scale = 10**(-2)
        self.avals = np.linspace(1./(1+self.zini), 1, 300)

        LCDMCosmology.__init__(self)
        self.updateParams([])


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varyOk): l.append(Ok_par)
        if (self.varyc):  l.append(self.c_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name   == "c":
                self.c = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False

        self.initialize()
        return True


    def RHS_a_hde(self, Omega, a):
        dOmega = -Omega*(Omega - 1)*(1 + 2*np.sqrt(Omega)/self.c)/a + 0
        return dOmega


    def compute_Ode(self, Ode_ini):
        Ode = Ode_ini*self.scale
        solution = odeint(self.RHS_a_hde, Ode, self.avals, h0=1E-5)
        return solution


    def ini_sol(self, Ode_ini):
        diference = self.compute_Ode(Ode_ini)[-1] - (1-self.Om-self.Ok)
        return diference


    def initialize(self):
        """
        Main method that searches the initial conditions for a given model.
        """
        ini_val = optimize.newton(self.ini_sol, 1)
        Ode = self.compute_Ode(ini_val)
        self.Omega_hde = interp1d(self.avals, Ode[:, 0])

        ## Add a flag in case the ini condition isn't found, i.e. c<0.4
        return True




    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self, a):
        if a< 1./(1+self.zini):
            Ode =  (1-self.Om-self.Ok)
        else:
            Ode = self.Omega_hde(a)

        return self.Omrad/a**4 + self.Ocb/a**3 + self.Ok/a**2 + Ode