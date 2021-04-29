
## This is Quintom cosmology
# with the interaction term

from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import *
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import numpy as np


class QuintomCosmology(LCDMCosmology):
    def __init__(self, vary_mquin=False, vary_mphan=False, vary_iniphi=False, vary_coupling=False, vary_Ok=False):

        """
        It better to start the chains at masses equal one, othewise may take much longer.
        mphi - mass of Quintessence.
        mpsi - mass of Phantom.
        """

        self.vary_mquin = vary_mquin
        self.vary_mphan = vary_mphan
        self.vary_coupling = vary_coupling
        self.vary_iniphi = vary_iniphi
        self.vary_Ok = vary_Ok

        self.mquin = 0 if (vary_mphan and (not vary_mquin)) else mquin_par.value
        self.mphan = 0 if (vary_mquin and (not vary_mphan)) else mphan_par.value
        self.coupling = 0 if (not vary_coupling) else coupling_par.value
        self.iniphi = 0
        self.Ok = Ok_par.value

        self.zvals = np.linspace(0, 3, 100)
        self.lna = np.linspace(-10, 0, 500)
        self.z = np.exp(-self.lna) - 1.

        #whether we rather printing all
        self.chatty = False

        self.Ol = 0

        LCDMCosmology.__init__(self, mnu=0)

        self.updateParams([])



    # Free parameters. We add Ok on top of LCDM ones (we inherit LCDM).
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if self.vary_mquin: l.append(mquin_par)
        if self.vary_mphan: l.append(mphan_par)
        if self.vary_iniphi: l.append(iniphi_par)
        if self.vary_coupling: l.append(coupling_par)
        if self.vary_Ok: l.append(Ok_par)
        return l



    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False

        for p in pars:
            if p.name  == "mquin":
                self.mquin = p.value
            elif p.name == "mphan":
                self.mphan = p.value
            elif p.name == "iniphi":
                self.iniphi = p.value
            elif p.name == 'beta':
                self.coupling = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False


        # Searches initial conditions and stores the Hubble function.
        if self.chatty:
            print('-'*10)
            print ('mphi={}, mphan={}, coupling={}'.format(self.mquin, self.mphan, self.coupling))

        # Main method.
        self.initialize()

        return True


    def sf_potential(self, phi, psi, select):
        """ The quadratic potential.
        Select:
            (0) the potential,
            ('phi') derivative wrt phi,
            ('psi') derivative wrt psi.
        """

        if select == 0:
            potential = 0.5*(self.mquin*phi)**2 + 0.5*(self.mphan*psi)**2 + self.coupling*(phi*psi)**2
        elif select == 'phi':
            potential = self.mquin**2*phi + 2*self.coupling*phi*psi**2
        elif select == 'psi':
            potential = self.mphan**2*psi + 2*self.coupling*psi*phi**2
        return potential




    def sf_rho(self, variables):
        """ Computes the density for the scalar field. """
        phi, x_phi, psi, x_psi = variables
        sf_rho = x_phi**2 - x_psi**2 + self.sf_potential(phi, psi, 0)/(3*self.h**2)
        return sf_rho



    def hubble(self, lna, variables=None, use_sf=True):
        """
        Computes the Hubble function for either a scalar-field or cosmological constant.

        Parameters
        ----------
        lna : array
            Interval of scale factor.
        variables : list , optional
            Array [phi, dotPhi, psi, dotPhi]. Defaults to ``None`` or LCDM model.
        use_sf : Boolean , optional
            Defaults to ``True`` or SFDE model, otherwise use LCDM.

        Returns
        -------
        hubble : array
            The ``hubble`` function value at the range ``lna``.

        """
        if use_sf:
            omega_de = self.sf_rho(variables)
        else:
            omega_de = 1.0 - self.Om - self.Ok - self.Ol

        a = np.exp(lna)
        hubble = self.h*np.sqrt(self.Omrad/a**4 + self.Ocb/a**3 + self.Ok/a**2 + self.Ol + omega_de)
        return hubble



    def right_hand_side(self, variables, lna):
        """ Right hand side of the dynamical system. """

        factor = np.sqrt(6)*self.h

        # Compute hubble function.
        hubble  = self.hubble(lna, variables=variables)

        phi, x_phi, psi, x_psi = variables

        # Right hand side of the dynamical system.
        rhs = [factor*x_phi/hubble,
               -3*x_phi - self.sf_potential(phi, psi, 'phi')/(factor*hubble),
               factor*x_psi/hubble,
               -3*x_psi + self.sf_potential(phi, psi, 'psi')/(factor*hubble)]
        return rhs


    def solve_eqns(self, *y0):
        """
        Input initial conditions for [phi, dotPhi, psi, dotPhi] and returns
        the solution of the Klein-Gordon equations.
        """
        solution = odeint(self.right_hand_side, y0, self.lna, h0=1E-10)
        return solution




    def compute_omega_de(self, ini_guess):
        """
        Given an guess for the initial condition of the field
        return the solution (if any) of the Klein-Gordon equation.
        """

        # Still figuring out initial conditions for two fields.
        if self.vary_mquin and (self.mphan == 0) and (self.coupling == 0):
            phi_ini, psi_ini = 10**ini_guess, 0
        elif self.vary_mphan and (self.mquin == 0) and (self.coupling == 0):
            phi_ini, psi_ini = 0, 10**ini_guess
        else:
            if self.vary_iniphi:
                phi_ini, psi_ini = 10**ini_guess, self.iniphi
            else:
                phi_ini, psi_ini = 10**ini_guess, 10**ini_guess

        # Find the solution for such a guess.
        solution = self.solve_eqns(phi_ini, 0.0, psi_ini, 0.0)

        # Solution at a=1.
        sf_rho = self.sf_rho(solution[-1])
        omega_de = sf_rho*(self.h/self.hubble(self.lna[-1], solution[-1]))**2

        tolerance = (1 - self.Omrad - self.Ocb - self.Ok - self.Ol) - omega_de
        return solution.T, tolerance



    def find_initial_phi(self, low_guess=-2, high_guess=3, tolerance=5E-3):
        """
        Using the bisection method, searches for the initial condition of the
        field such that \Omega_DE today is the value given by observations.

        Parameters
        ----------
        low_guess : float , optional
            Minimum value of the field.
        high_guess : float , optional
            Maximum value of the field.
        tolerance : float , optional
            Tolerance or stopping criteria, abs(Omega_observable - Omega_theory).

        Returns
        -------
        mid_guess : float
            Initial condition for the field: ``mid_guess``.

        solution : list
            Found solution (if any): ``solution`` which contains [phi, xPhi, psi, xPsi].

        """

        mid_guess = 0.5*(low_guess + high_guess)

        while (high_guess - low_guess) > tolerance:
            # Compute the solution for an initial guess.
            solution, current_tolerance = self.compute_omega_de(mid_guess)
            low_tolerance = self.compute_omega_de(low_guess)[1]

            if self.chatty:
                print('-done- phi_0={}, error={}'.format(mid_guess, current_tolerance))

            if np.abs(current_tolerance) < tolerance:
                return mid_guess, solution

            elif low_tolerance*current_tolerance < 0:
                high_guess = mid_guess
            else:
                low_guess  = mid_guess
            mid_guess = (low_guess + high_guess)*0.5

        # Check whether rho is constant or nearby.
        grad = np.abs(np.gradient(self.sf_rho(solution))).max()
        mid_guess = -1 if grad < 1.0E-2 else 0

        if self.chatty:
            if mid_guess == -1:
                print('looks a lot like LCDM', mid_guess, self.mquin, self.mphan, self.iniphi)
        if mid_guess == 0:
            if self.chatty: print('-- No solution found!', mid_guess, self.mquin, self.mphan, self.iniphi)
        return mid_guess, solution



    def initialize(self):
        """
        Main method that searches the initial conditions and computes [phi, xPhi, psi, xPsi]
        for a given scalar field potential.
        """

        high_guess = 2 if self.coupling >0 else 1

        #It's slower than newton, but finds more solutions
        self.phi_ini, self.solution = self.find_initial_phi(-2, high_guess)

        if self.phi_ini == 0:
            # Solution couldn't be found, then reject that point.
            hubble = np.zeros(len(self.lna))

        elif self.phi_ini == -1:
            # Solution couldn't be found, but it is very close to the cosmological constant.
            hubble = self.hubble(self.lna, use_sf=False)

        else:
            # Solution found.
            hubble = self.hubble(self.lna, self.solution)

        self.sf_hubble = interp1d(self.lna, hubble)
        return True





    def RHSquared_a(self, a):
        """ This is relative hsquared as a function of a, i.e. H(z)^2/H(z=0)^2. """

        if (1./a-1 < self.zvals[-1]):
            hubble = (self.sf_hubble(np.log(a))/self.h)**2.
        else:
            hubble = (self.hubble(np.log(a), use_sf=False)/self.h)**2
        return hubble


    #-------------------------
    # External methods useful for plotting or testing.


    def logatoz(self, func):
        # change function from lna to z
        tmp     = interp1d(self.lna, func)
        functmp = tmp(self.lna)
        return np.interp(self.zvals, self.z[::-1], functmp[::-1])


    def Hubble_a(self, a):
        hubble = 100*self.sf_hubble(np.log(a))
        return hubble


    def sf_eos(self, variables):
        phi, dot_phi, psi, dot_psi = variables

        w1 = dot_phi**2 - dot_psi**2 - self.sf_potential(phi, psi, 0)/(3*self.h**2)
        w2 = dot_phi**2 - dot_psi**2 + self.sf_potential(phi, psi, 0)/(3*self.h**2)
        return  w1/w2


    def call_functions(self):
        if self.phi_ini == 0:
            w_eos = np.zeros(len(self.lna))
        elif self.phi_ini == -1:
            w_eos = -1*np.ones(len(self.lna))
        else:
            w_eos = self.sf_eos(self.solution)

        self.w_eos = interp1d(self.lna, w_eos)
        return True



    def w_de(self, a):
        lna = np.log(a)
        return self.w_eos(lna)





