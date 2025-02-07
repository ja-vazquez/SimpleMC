#
# This module calculates likelihood for a Generic DATA.
#
from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
from simplemc.setup_logger import logger
import scipy.linalg as la
import scipy as sp
from simplemc.setup_logger import cdir
import numpy as np


# uncoment lines for use a covariance matrix
class SimpleLikelihood(BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename, fn='generic'):
        BaseLikelihood.__init__(self, name)
        # print("Loading ", values_filename)
        logger.info("Loading {}".format(values_filename))
        data = np.loadtxt(values_filename)
        self.xx  = data[:,0]
        self.yy  = data[:,1]
        self.cov = np.loadtxt(cov_filename,skiprows=0)
        assert(len(self.cov) == len(self.xx))
        self.icov = la.inv(self.cov)
        self.fn = fn

    def loglike(self):
        #delta is the difference between theory and data
        if self.fn == "generic":
            tvec  = np.array([self.theory_.genericModel(z) for z in self.xx])
        elif self.fn == "h":
            tvec = np.array([100.0 * self.theory_.h * np.sqrt(self.theory_.RHSquared_a(1.0 / (1 + z))) for z in self.xx])
        elif self.fn == "fs8":
            tvec = np.array([self.theory_.fs8(z) for z in self.xx])
        elif self.fn == "distance_mod":
            tvec = np.array([self.theory_.distance_modulus(z) for z in self.xx])

        delta = self.yy - tvec
        return -0.5*np.dot(delta, np.dot(self.icov, delta))



class StraightLine(SimpleLikelihood):
     def __init__(self):
         SimpleLikelihood.__init__(self,"GenericData", cdir+"/data/line_data.txt",
                                   cdir+"/data/line_cov.txt")


class GenericLikelihood(SimpleLikelihood):
    def __init__(self, path_to_data, path_to_cov, fn):
        SimpleLikelihood.__init__(self, "GenericData", path_to_data,
                                  path_to_cov, fn)
