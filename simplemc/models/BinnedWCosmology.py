

from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.Parameter import Parameter
from scipy.interpolate import interp1d
from scipy.integrate import quad
import numpy as np

## Binned cosmology, where the DE eqn of state is assumed to be a set of bins
#  with varying amplitudes and fix positions.
#  Nevertheless, this class may be deprecated and substituted by the code use
#  for writing the paper https://arxiv.org/abs/2111.10457
##

class BinnedWCosmology(LCDMCosmology):
    def __init__(self, dz=0.2, zmax=1.0):
        """
        This class corresponds to a CDM cosmology with binned w.
        Still testing the file but it seems to be working just fine
        Parameters
        ----------
        dz : float
            Step size for the position of the bins.

        zmax : float
            Maximum redshift to use for the reconstruction.

        Returns
        -------

        """
        # Bunch of parameters for amplitudes and positions of the bins.
        self.zbins = np.arange(0, zmax, dz)
        self.Nb = len(self.zbins)
        self.wvals = np.ones(self.Nb)*-1.0
        self.pnames = ["w%i" % i for i in range(self.Nb)]

        LCDMCosmology.__init__(self)
        self.integrateOmega()


    # My free parameters. We use a flat cosmology.
    def freeParameters(self):
        wpars = [Parameter(name, self.wvals[i], err=0.05)
                 for i, name in enumerate(self.pnames)]
        return LCDMCosmology.freeParameters(self) + wpars


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        gotone = False
        for p in pars:
            ## Something's happening here, check it later.
            print ('**', p.name, self.pnames)
            i = self.pnames.index(p.name)

            if i > 0:
                self.wvals[i] = p.value
                gotone = True
        if gotone:
            self.integrateOmega()
        return True


    def integrateOmega(self):
        abins = np.hstack((1./(1+self.zbins), [1e-4]))
        w = np.hstack((self.wvals, [self.wvals[-1]]))
        itg = interp1d(np.log(abins), 3*(1 + w))
        oabins = np.hstack((np.logspace(-4, -1, 10), np.linspace(0.1, 1, 100)))
        olnrho = [quad(itg, np.log(a), 0)[0] for a in oabins]
        print(1/oabins**4)
        print(np.exp(olnrho))
        self.DEomega = interp1d(oabins, np.exp(olnrho))



    # This is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2.
    def RHSquared_a(self, a):
        NuContrib = self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om)*self.DEomega(a))
