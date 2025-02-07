

from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import scipy.linalg as la
import scipy as sp
import numpy as np
from simplemc.setup_logger import cdir


#uncoment lines for use a covariance matrix
class RotationCurvesLikelihood(BaseLikelihood):
    def __init__(self, name, values_filename, cov_filename):
        """
        Class that includes the likelihood for a  Particular rotation curve: NGC2403
        Parameters
        ----------
        name
        values_filename
        cov_filename

        Returns
        -------

        """
        BaseLikelihood.__init__(self,name)
        print("Loading ", values_filename)
        data = np.loadtxt(values_filename)
        self.xx  = data[:,0]
        self.yy  = data[:,1]
        self.cov = np.loadtxt(cov_filename,skiprows=0)
        assert(len(self.cov) == len(self.xx))
        self.icov = la.inv(self.cov)



    def loglike(self):
        #delta is the difference between theory and data
        tvec  = np.array([self.theory_.rotation(r) for r in self.xx])

        delta = self.yy - tvec
        return -0.5*np.dot(delta, np.dot(self.icov,delta))



class RotationCurvesLike(RotationCurvesLikelihood):
    def __init__(self):
        RotationCurvesLikelihood.__init__(self, "RotCurves", cdir+"/data/NGC2403.curve.02",
                                          cdir+"/data/NGC2403.curve.02-cov.txt")
