from scipy.interpolate import InterpolatedUnivariateSpline
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
    def __init__(self, varyOk=False, varychde=True,
                 nodes=3  # numer of nodes used in the interpolation
                 ):
        # Holographic parameter
        # Notation: capital C = 3c**2*M_p, here, we use little c

        self.nnodes = nodes

        self.varyOk = varyOk
        self.varychde = varychde

        self.c_par = Parameter("c", 0.7, 0.1, (0.5, 2.0), "c")

        self.Ok = Ok_par.value
        self.c_hde = self.c_par.value

        # zend is the last point in linear-interpolation
        # however is the first point of the last step in the binning/tanh version
        self.zini = 0.0
        self.zend = 3.0

        # range to perform the interpolation
        self.zvals = np.linspace(self.zini, self.zend, 100)


        mean = -2
        priors = (-1, 3)
        sigma = 0.2

        self.pname = 'amp_'
        self.params = [Parameter(self.pname+str(i), mean, sigma, priors, self.pname+str(i)) for i in range(nodes)]
        self.pvals = [i.value for i in self.params]

        LCDMCosmology.__init__(self, disable_radiation=True)
        self.updateParams([])


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        l += self.params
        if (self.varychde):  l.append(self.c_par)
        if (self.varyOk): l.append(Ok_par)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False

        for i, p in enumerate(pars):
            if p.name == (self.pname+str(i)):
                self.pvals[i] = p.value
            if p.name == "c":
                self.c_hde = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False

        if self.nnodes > 1:
            self.ini_function()
        self.initialize()
        return True




    def ini_function(self):
        self.y = self.pvals
        self.mnodes = self.nnodes-1

        delta = (self.zend - self.zini)/(self.mnodes)
        self.x = [self.zini + delta*i for i in range(self.mnodes+1)]

        self.spline = InterpolatedUnivariateSpline(self.x, self.y, k=1)

        return True



    def ffunc(self, z):
        # = 0 for CC, = 2 for Barrow
        # f = \Delta-2
        # function = self.alpha_hde + self.beta_hde*z
        if self.nnodes==1:
            funct = self.pvals[0]
        else:
            funct = self.spline(z)
        return funct


    def dffunc(self, z):
        # dfunction = self.beta_hde
        if self.nnodes==1:
            dfunc = 0
        else:
            dfucn = x
        return dfunc


    def Q_term(self, z):
        delta = self.ffunc(z) + 2
        term1 = 1./(delta - 2)

        Q = (2-delta)*(self.c_hde)**(2*term1)*(100*self.h*np.sqrt(self.Om))**(-delta*term1)
        return Q


    def extra_term(self, z, Ode):
        delta = self.ffunc(z) + 2
        term1 = 1./(delta - 2)

        return (1 - Ode)**(0.5*delta*term1)*Ode**(-term1)


    def deriv_func(self, z, Ode):
        delta = self.ffunc(z) + 2
        term1 = 1. / (delta - 2)

        qterm= self.c_hde**2*(1+z)**(-3)/(100*self.h)**2/self.Om
        # this is where the sign changes compared to previous works
        value= (1 + z)*self.dffunc(z)*term1*np.log(qterm*(1 - Ode)/Ode)
        return value


    # Right hand side of the equations
    def RHS_hde(self, Ode, z):
        delta = self.ffunc(z) + 2
        term1 = 1./(delta - 2)

        term_cte = 1 + delta + self.Q_term(z)*self.extra_term(z, Ode)*(1 + z)**(-1.5*delta*term1)
        fact = term_cte + self.deriv_func(z, Ode)
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
        if z>self.zini:
            Ode =  (1-self.Om-self.Ok)
            hubble = self.Omrad/a**4 + self.Ocb/a**3 + self.Ok/a**2 + Ode
        else:
            hubble = self.Om*(1+z)**3/(1-self.Ode(z))

        return hubble


    def EoS(self, z):
        Ode = self.Ode(z)
        delta = self.ffunc(z)

        tmp = (1+z)**(-1.5*delta/(2-delta))
        w = -(1+delta) - self.Q_term(z)*self.extra_term(z, Ode)*tmp \
            + self.deriv_func(z, Ode)
        return w/3.
