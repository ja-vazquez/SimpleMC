from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d
from simplemc import cdir


class PantheonSNLikelihood(BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename, ninterp=150):
        """
        This module calculates likelihood for Pantheon datasets.
        Parameters
        ----------
        name: str
            name of the likelihood
        values_filename: str
            directory and name of the data file
        cov_filename: str
            directory and name of the covariance matrix file
        ninterp: int
        """
        # first read data file
        self.name_ = name
        BaseLikelihood.__init__(self, name)
        print("Loading", values_filename)
        da = np.loadtxt(values_filename, skiprows=1, usecols=(1, 2, 4, 5))
        self.zcmb = da[:, 0]
        self.zhelio = da[:, 1]
        self.mag = da[:, 2]
        self.dmag = da[:, 3]
        self.N = len(self.mag)
        self.syscov = np.loadtxt(cov_filename, skiprows=1).reshape((self.N, self.N))
        self.cov = np.copy(self.syscov)
        self.cov[np.diag_indices_from(self.cov)] += self.dmag**2
        self.xdiag = 1/self.cov.diagonal()  # diagonal before marginalising constant
        # add marginalising over a constant
        self.cov += 3**2
        self.zmin = self.zcmb.min()
        self.zmax = self.zcmb.max()
        self.zmaxi = 1.1 ## we interpolate to 1.1 beyond that exact calc
        print("Pantheon SN: zmin=%f zmax=%f N=%i" % (self.zmin, self.zmax, self.N))
        self.zinter = np.linspace(1e-3, self.zmaxi, ninterp)
        self.icov = la.inv(self.cov)

    def loglike(self):
        # we will interpolate distance
        dist = interp1d(self.zinter, [self.theory_.distance_modulus(z) for z in self.zinter],
                        kind='cubic', bounds_error=False)(self.zcmb)
        who = np.where(self.zcmb > self.zmaxi)
        dist[who] = np.array([self.theory_.distance_modulus(z) for z in self.zcmb[who]])
        tvec = self.mag-dist

        # tvec = self.mag-np.array([self.theory_.distance_modulus(z) for z in self.zcmb])
        # print (tvec[:10])
        # first subtract a rought constant to stabilize marginaliztion of
        # intrinsic mag.
        tvec -= (tvec*self.xdiag).sum() / (self.xdiag.sum())
        # print(tvec[:10])
        chi2 = np.einsum('i,ij,j', tvec, self.icov, tvec)
        # print("chi2=",chi2)
        return -chi2/2


class PantheonSN(PantheonSNLikelihood):
    """
    Likelihood to full Pantheon SNIa compilation.
    """
    def __init__(self):
        PantheonSNLikelihood.__init__(self, "Pantheon", cdir+"/data/pantheon_lcparam_full_long_zhel.txt",
                                      cdir+"/data/pantheon_sys_full_long.txt")


class BinnedPantheon(PantheonSNLikelihood):
    """
    Likelihood to binned Pantheon dataset.
    """
    def __init__(self):
        PantheonSNLikelihood.__init__(self, "BPantheon", cdir+"/data/binned_pantheon.txt",
                                      cdir+"/data/binned_cov_pantheon.txt")