import numpy as np
from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
from simplemc.models.LCDMCosmology import LCDMCosmology
import scipy.linalg as la
import numpy as sp
from simplemc.setup_logger import cdir


class DESIDR2BAOLikelihood(BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename, fidtheory):
        """
        This module calculates likelihood for the consensus DESIDR2-BAO
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
        da = sp.loadtxt(values_filename, usecols = (0 ,1 ,2))
        self.zs    = da[:, 0]
        self.DM_DH = da[:, 1]
        self.type  = da[:, 2]

        print("Loading covariance DESIDR2BAO")
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
        tvec = sp.array(tvec)
        tvec += 0
        delta = tvec - self.DM_DH
        return -sp.dot(delta, sp.dot(self.icov, delta))/2.0

class DESIDR2BAO(DESIDR2BAOLikelihood):
    """
    Likelihood to full DESIBAO compilation.
    """
    def __init__(self):
        obh2 = 0.022
        Om   = 0.31
        h    = 0.676
        mnu  = 0.06
        fidtheory = LCDMCosmology(obh2, Om, h, mnu)
        DESIDR2BAOLikelihood.__init__(self, "DESIDR2BAO", cdir+"/data/desi_gaussian_bao_ALL_GCcomb_mean.txt",
                                          cdir+"/data/desi_gaussian_bao_ALL_GCcomb_cov.txt", fidtheory)
