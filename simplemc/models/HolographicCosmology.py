from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
from simplemc.cosmo.paramDefs import Ok_par
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
        self.c_par = Parameter("c", 0.7, 0.1, (0.5, 1.1), "c")

        self.varyOk = varyOk
        self.varyc  = varyc

        self.Ok = Ok_par.value
        self.c  = self.c_par.value

        # This value is quite related to the initial z
        self.zini = 3
        #self.scale = 10**(-2)
        #self.avals = np.linspace(1./(1+self.zini), 1, 300)
        self.zvals = np.linspace(0, self.zini, 50)
        #self.zhde = np.linspace(0, 10, 500)

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


    def RHS_hde(self, vals, z):
        Omega = vals
        fact = (1 + 2*np.sqrt(Omega)/self.c)
        dOmega = -Omega*(1 - Omega)*fact/(1 + z) + 0
        #dlogE = -0.5*Omega*(-3/Omega + fact)/(1 + z)
        #dOmega = -Omega*(Omega - 1)*(1 + 2*np.sqrt(Omega)/self.c)/a + 0
        return dOmega


    #def compute_Ode(self, ini_vals):
    #    solution = odeint(self.RHS_hde, ini_vals, self.zvals)
        #Ode = Ode_ini*self.scale
        #solution = odeint(self.RHS_a_hde, ini_vals, self.avals, h0=1E-5)
    #    return solution

    #For shooting
    #def ini_sol(self, Ode_ini):
    #    diference = self.compute_Ode(Ode_ini)[-1] - (1-self.Om-self.Ok)
    #    return diference


    def initialize(self):
        """
        Main method that searches the initial conditions for a given model.
        """
        #ini_val = optimize.newton(self.ini_sol, 1)
        Ode0 = (1 - self.Om - self.Ok)
        #logE0 = 0

        #ini_vals = Ode0
        result_E = odeint(self.RHS_hde, Ode0, self.zvals)
        #Ode = np.exp(result_E[:, 0])
        self.Ode = interp1d(self.zvals, result_E[:, 0])
        #self.Omega_hde = interp1d(self.avals, Ode[:, 0])

        ## Add a flag in case the ini condition isn't found, i.e. c<0.4
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