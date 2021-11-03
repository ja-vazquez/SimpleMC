

from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import scipy.linalg as la
import scipy as sp
from simplemc import cdir


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
        data = sp.loadtxt(values_filename)
        self.xx  = data[:,0]
        self.yy  = data[:,1]
        self.cov = sp.loadtxt(cov_filename,skiprows=0)
        assert(len(self.cov) == len(self.xx))
        self.icov = la.inv(self.cov)



    def loglike(self):
        #delta is the difference between theory and data
        tvec  = sp.array([self.theory_.rotation(r) for r in self.xx])

        delta = self.yy - tvec
        return -0.5*sp.dot(delta, sp.dot(self.icov,delta))



class RotationCurvesLike(RotationCurvesLikelihood):
    def __init__(self):
        RotationCurvesLikelihood.__init__(self, "RotCurves", cdir+"/data/NGC2403.curve.02",
                                          cdir+"/data/NGC2403.curve.02-cov.txt")
