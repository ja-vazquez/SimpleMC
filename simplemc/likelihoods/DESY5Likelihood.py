from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d
from simplemc.setup_logger import cdir
import pandas as pd

class DESY5Likelihood(BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename, ninterp=150):
        """
        This module calculates likelihood for DESY5 datasets.
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
        da = pd.read_csv(values_filename)
        self.zcmb = da['zHD']
        self.zhelio = da['zHEL']
        self.mag = da['MU']
        self.dmag = da['MUERR_FINAL'] 
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
        print("DESY5: zmin=%f zmax=%f N=%i" % (self.zmin, self.zmax, self.N))
        self.zinter = np.linspace(1e-3, self.zmaxi, ninterp)
        self.icov = la.inv(self.cov)

    def loglike(self):
        # we will interpolate distance
        dist = interp1d(self.zinter, [self.theory_.distance_modulus(z) for z in self.zinter],
                        kind='cubic', bounds_error=False)(self.zcmb)
        who = np.where(self.zcmb > self.zmaxi)
        dist[who] = np.array([self.theory_.distance_modulus(z) for z in self.zcmb.loc[who]])
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


class DESY5(DESY5Likelihood):
    """
    Likelihood to full DESY5 compilation.
    """
    def __init__(self):
        DESY5Likelihood.__init__(self, "DESY5", cdir+"/data/DES-SN5YR_HD.csv",
                                      cdir+"/data/covsys_000.txt")


