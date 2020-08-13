#
# This module calculates likelihood for a Generic DATA.
#
from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
from simplemc import logger
import scipy.linalg as la
import scipy as sp

# uncoment lines for use a covariance matrix
class SimpleLikelihood (BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename, fn='generic'):
        BaseLikelihood.__init__(self, name)
        # print("Loading ", values_filename)
        logger.info("Loading {}".format(values_filename))
        data = sp.loadtxt(values_filename)
        self.xx  = data[:,0]
        self.yy  = data[:,1]
        self.cov = sp.loadtxt(cov_filename,skiprows=0)
        assert(len(self.cov) == len(self.xx))
        self.icov = la.inv(self.cov)
        self.fn = fn

    def loglike(self):
        #delta is the difference between theory and data
        if self.fn == "generic":
            tvec  = sp.array([self.theory_.genericModel(z) for z in self.xx])
        elif self.fn == "h":
            tvec = sp.array([100.0 * self.theory_.h * sp.sqrt(self.theory_.RHSquared_a(1.0 / (1 + z))) for z in self.zs])
        elif self.fn == "fs8":
            tvec = sp.array([self.theory_.fs8(z) for z in self.zs])
        elif self.fn == "distance_mod":
            tvec = sp.array([self.theory_.distance_modulus(z) for z in self.zs])

        delta = self.yy - tvec
        return -0.5*sp.dot(delta, sp.dot(self.icov, delta))



class StraightLine(SimpleLikelihood):
     def __init__(self):
         SimpleLikelihood.__init__(self,"GenericData","simplemc/data/line_data.txt",
            "simplemc/data/line_cov.txt")


class GenericLikelihood(SimpleLikelihood):
    def __init__(self, path_to_data, path_to_cov, fn):
        SimpleLikelihood.__init__(self, "GenericData", path_to_data,
                                  path_to_cov, fn)
