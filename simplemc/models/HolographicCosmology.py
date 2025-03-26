from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline
from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
from simplemc.cosmo.paramDefs import Ok_par
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import numpy as np
from scipy import constants


class HolographicCosmology(LCDMCosmology):
    """
        This is Holographic cosmology.
        This class inherits LCDMCosmology class as the rest of the cosmological
        models already included in SimpleMC.

        :param varyc: variable w0 parameter

    """
    def __init__(self, varyOk=False, varychde=True, mean=-2,
                 nodes=3,   # numer of nodes used in the interpolation
                 interp= 'lineal'
                 ):
        # Holographic parameter
        # Notation: capital C = 3c**2*M_p, here, we use little c
        # Nodes: =0 delta constant, = 1 varying delta, = 2 two nodes, etc.

        self.nnodes = nodes
        # =1 lineal, 3 = cubic
        self.interp = 1 if interp == 'lineal' else 3

        self.varyOk = varyOk
        self.varychde = varychde

        self.c_par = Parameter("c", 1.0, 0.1, (0.5, 2.0), "c")

        self.Ok = Ok_par.value
        self.c_hde = self.c_par.value

        # zend is the last point in linear-interpolation
        # however is the first point of the last step in the binning/tanh version
        self.zini = 0.0
        self.zend = 2.5

        # range to perform the interpolation
        self.zvals = np.linspace(self.zini, self.zend, 100)

        # = 0 for CC, = -2 for Barrow
        mean = mean
        priors = (-3, 1)
        sigma = 0.2

        self.pname = 'amp_'
        amps = 1 if nodes <= 1 else nodes
        self.params = [Parameter(self.pname+str(i), mean, sigma, priors, self.pname+str(i)) for i in range(amps)]
        self.pvals = [i.value for i in self.params]

        LCDMCosmology.__init__(self, disable_radiation=True)
        self.updateParams([])


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varychde): l.append(self.c_par)
        if (self.varyOk): l.append(Ok_par)
        if self.nnodes > 0:
            l += self.params
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False

        for p in pars:
            if self.pname in p.name:
                for i in range(self.nnodes):
                    if p.name == (self.pname+str(i)):
                        self.pvals[i] = p.value
            if p.name == "c":
                self.c_hde = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False

        self.ini_function()
        self.initialize()
        return True



    def ini_function(self):
        if self.nnodes <= 1:
            x = [self.zini, self.zend]
            y = self.pvals[0]*np.ones(2)
            self.ffunc = InterpolatedUnivariateSpline(x, y, k=1)
        else:
            y = self.pvals
            delta = (self.zend - self.zini)/(self.nnodes-1)
            x = [self.zini + delta*i for i in range(self.nnodes)]

            self.ffunc = InterpolatedUnivariateSpline(x, y, k=self.interp)
            self.dffunc = self.ffunc.derivative(n=1)
        return True



    def extra_term(self, z, Ode):
        # Bear in mind that f=D-2
        q_term = self.c_hde**2/(self.Om*(1 + z)**3)
        return q_term*(1 - Ode)/Ode



    # Right hand side of the equations
    def RHS_hde(self, Ode, z):
        func = self.ffunc(z)
        fact = 3 + func

        if np.abs(func) > 0.005:
            ex_term = self.extra_term(z, Ode)
            fact += -1*func*np.sqrt(Ode)/self.c_hde*ex_term**(0.5*(func+2)/func)

            # derivative term
            if self.nnodes > 1:
                fact += (1+z)/func*self.dffunc(z)*np.log(ex_term)
        dOmega = -Ode*(1 - Ode)*fact/(1 + z)
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
        if z>self.zend:
            Ode = (1 - self.Om - self.Ok)
            hubble = self.Omrad/a**4 + self.Ocb/a**3 + self.Ok/a**2 + Ode
        else:
            hubble = self.Om*(1+z)**3/(1 - self.Ode(z))

        return hubble


    def EoS(self, z):
        Ode = self.Ode(z)
        func = self.ffunc(z)
        fact = 3 + func

        if np.abs(self.ffunc(z)) > 0.005:
            ex_term = self.extra_term(z, Ode)
            fact += -1*func*np.sqrt(Ode)/self.c_hde*ex_term**(0.5*(func+2)/func)

            if self.nnodes > 1:
                fact += (1+z)/func*self.dffunc(z)*np.log(ex_term)
        return -fact/3.
