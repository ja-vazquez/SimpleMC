#
# This module calculates likelihood for a Generic DATA.
#
from BaseLikelihood import BaseLikelihood
import scipy.linalg as la
import scipy as sp

#uncoment lines for use a covariance matrix
class CompressedGenericLikelihood (BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename):
        BaseLikelihood.__init__(self,name)
        print("Loading ", values_filename)
        data = sp.loadtxt(values_filename)
        self.xx  = data[:,0]
        self.yy  = data[:,1]
        self.cov = sp.loadtxt(cov_filename,skiprows=0)
        assert(len(self.cov) == len(self.xx))
        self.icov = la.inv(self.cov)



    def loglike(self):
        #delta is the difference between theory and data
        tvec  = sp.array([self.theory_.genericModel(z) for z in self.xx])

        delta = self.yy - tvec
        return -0.5*sp.dot(delta, sp.dot(self.icov,delta))



class StraightLine(CompressedGenericLikelihood):
    def __init__(self):
        CompressedGenericLikelihood.__init__(self,"GenericData","data/line_data.txt",
            "data/line_cov.txt")
