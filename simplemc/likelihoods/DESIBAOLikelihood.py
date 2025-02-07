
import numpy as np
from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import scipy.linalg as la
import numpy as sp



class DESIBAOLikelihood(BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename, fidtheory):
        """
        This module calculates likelihood for the consensus DESI-BAO
        ----------
        name
        values_filename
        cov_filename
        fidtheory

        Returns
        -------

        """
        BaseLikelihood.__init__(self ,name)

        self.rd = fidtheory.rd
        print("Loading ", values_filename)
        da = np.loadtxt(values_filename, usecols = (0 ,1 ,2))
        self.zs    = da[:, 0]
        self.DM_DH = da[:, 1]
        self.type  = da[:, 2]

        print("Loading covariance DESIBAO")
        cov = np.loadtxt(cov_filename)
        assert(len(cov) == len(self.zs))
        print("Adding marginalising constant")
        cov += 3**2
        self.icov = la.inv(cov)

    def loglike(self):
        tvec = []
        for i, z in enumerate(self.zs):
            if self.type[i] == 4:
                tvec.append(self.theory_.DaOverrd(z))
            elif self.type[i] == 5:
                tvec.append(self.theory_.HIOverrd(z))
            elif self.type[i] == 3:
                tvec.append(self.theory_.DVOverrd(z))
        tvec = np.array(tvec)
        tvec += 0
        delta = tvec - self.DM_DH
        return -np.dot(delta, np.dot(self.icov, delta))/2.0