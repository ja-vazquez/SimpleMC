#
# This module calculates likelihood for a Generic DATA.
#
from BaseLikelihood import *
from scipy import *
import scipy.linalg as la
import numpy as np
#uncoment lines for use a covariance matrix
class CompressedGenericLikelihood (BaseLikelihood):
    def __init__(self,name,values_filename, cov_filename):
        BaseLikelihood.__init__(self,name)
        print("Loading ",values_filename)
        data=np.loadtxt(values_filename)
        self.xx=data[:,0]
        self.yy=data[:,1]
        self.cov=np.loadtxt(cov_filename,skiprows=0)
        assert(len(self.cov)==len(self.xx))
        self.icov=la.inv(self.cov)

    def loglike(self):
        #delta is the difference between theory and data
        tvec=array([self.theory_.genericModel(z) for z in self.xx])
        #norm = -0.5*len(self.xx)*np.log(2*np.pi)-0.5*np.log(la.det(cov))
        norm = -0.5*len(self.xx)*np.log(np.pi)-0.5*np.log(la.det(self.cov))
        ## This is the factor that we need to correct
        ## note that in principle this shouldn't matter too much, we will marginalise over this
        #delta=tvec-self.yy
        delta=self.yy-tvec
        return -0.5*dot(delta,dot(self.icov,delta))



class GenericLikelihood(CompressedGenericLikelihood):
    def __init__(self):
        CompressedGenericLikelihood.__init__(self,"GenericData","data/linedata.txt",
            "data/cov_line.txt")
