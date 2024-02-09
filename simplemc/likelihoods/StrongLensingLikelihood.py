

from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import scipy.linalg as la
import scipy as sp
from simplemc.setup_logger import cdir
import numpy as np

class StrongLensingLikelihood(BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename):
        """
        This module calculates likelihood for Strong Lensing.
        based on the CompressHD file
        Parameters
        ----------
        name
        values_filename
        cov_filename

        Returns
        -------

        """
        BaseLikelihood.__init__(self,name)
        print("Loading ",values_filename)
        da = sp.loadtxt(values_filename)
        self.zl = da[:,0]
        self.zs = da[:, 1]
        self.theta_E = da[:,2]
        self.sigma = da[:,4]

        print("Loading ", cov_filename)
        cov = sp.loadtxt(cov_filename, skiprows=1)
        assert(len(cov) == len(self.zs))
        cov += 3 ** 2

        vals, vecs = la.eig(cov)
        vals = sorted(sp.real(vals))
        print("Eigenvalues of cov matrix:", vals[0:3],'...',vals[-1])
        print("Adding marginalising constant")

        self.icov = la.inv(cov)


    def loglike(self):
        c = 3*10**5
        fact = (c**2)*(np.pi/(3600*180))/(4*np.pi)
        tvec = sp.array([self.theory_.distance_ratio(zl0, zs0) for zl0, zs0 in zip(self.zl, self.zs)])
        tobs = sp.array([fact*theta/sigma**2 for theta, sigma in zip(self.theta_E, self.sigma)])

        tvec += 0
        delta = tvec - tobs
        loglike= -sp.dot(delta, sp.dot(self.icov, delta))/2.0
        return loglike

class StrongLensing(StrongLensingLikelihood):
    # data from https://arxiv.org/abs/1906.04107
    def __init__(self):
        StrongLensingLikelihood.__init__(self, "SL", cdir+"/data/Strong_Lensing_196.txt",
                                            cdir+"/data/Strong_Lensing_cov_196.txt")

