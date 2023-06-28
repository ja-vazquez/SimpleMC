from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d
from simplemc import cdir
import pandas as pd


class PantheonPlusSNLikelihood(BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename, ninterp=150):
        """
        This module calculates likelihood for Pantheon plus dataset without SH0ES.

        Likelihood based on:
        Pantheon+_only_cosmosis_likelihood.py from the
        https://github.com/PantheonPlusSH0ES repository.

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
        self.name_ = name
        BaseLikelihood.__init__(self, name)
        print("Loading", values_filename)
        data = pd.read_csv(values_filename, delim_whitespace=True)
        self.origlen = len(data)
        self.ww = (data['zHD'] > 0.01)
        self.zcmb = data['zHD'][self.ww].values  # use the vpec corrected redshift for zCMB
        self.zhelio = data['zHEL'][self.ww].values
        self.mag = data['m_b_corr'][self.ww].values
        self.N = len(self.mag)

        filename = cov_filename
        print("Loading covariance from {}".format(filename))
        f = open(filename)
        line = f.readline()
        n = int(len(self.zcmb))
        C = np.zeros((n,n))
        ii = -1
        jj = -1
        mine = 999
        maxe = -999
        for i in range(self.origlen):
            jj = -1
            if self.ww[i]:
                ii += 1
            for j in range(self.origlen):
                if self.ww[j]:
                    jj += 1
                val = float(f.readline())
                if self.ww[i]:
                    if self.ww[j]:
                        C[ii,jj] = val
        f.close()
        print('Done')
        self.cov = C
        self.xdiag = 1/self.cov.diagonal()  # diagonal before marginalising constant
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
        tvec = self.mag - dist

        # tvec = self.mag-np.array([self.theory_.distance_modulus(z) for z in self.zcmb])
        # print (tvec[:10])
        # first subtract a rought constant to stabilize marginaliztion of
        # intrinsic mag.
        tvec -= (tvec * self.xdiag).sum() / (self.xdiag.sum())
        # print(tvec[:10])
        chi2 = np.einsum('i,ij,j', tvec, self.icov, tvec)
        # print("chi2=",chi2)
        return -chi2 / 2
class PantheonPlus(PantheonPlusSNLikelihood):
    """
    Likelihood to full Pantheon+ SNIa compilation.

    Data from the 4_DISTANCES_AND_COVAR folder on the
    https://github.com/PantheonPlusSH0ES repo.
    """
    def __init__(self):
        PantheonPlusSNLikelihood.__init__(self, "PantheonPlus", cdir+"/data/Pantheon+SH0ES.dat",
                                          cdir+"/data/Pantheon+SH0ES_STAT+SYS.cov")
