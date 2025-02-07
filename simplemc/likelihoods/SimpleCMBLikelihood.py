

from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import scipy.linalg as la
import numpy as np
import scipy as sp

class SimpleCMBLikelihood(BaseLikelihood):
    def __init__(self, name, mean, cov, kill_Da=False, kill_rd=False):
        """
        This is inspired by WangWang CMB, but ultimately simpler
        where we compress to just obh2, obh2+och2 and thetas
        Parameters
        ----------
        name
        mean
        cov
        kill_Da
        kill_rd

        Returns
        -------

        """
        print("Initializing CMB likelihood:", name)
        cov = np.array(cov)
        if kill_Da:
            cov[2, :] = 0.0
            cov[:, 2] = 0.0
            cov[2, 2] = 1e10
            name += "_noDa"
        if kill_rd:
            cov[0:2, 0:2] = 0.0
            cov[0, 0] = 1e10
            cov[1, 1] = 1e10
            name += "_noRd"
        BaseLikelihood.__init__(self, name)
        self.mean = np.array(mean)
        self.cov  = cov
        self.icov = la.inv(cov)


    def setTheory(self, theory):
        self.theory_ = theory
        self.theory_.setNoObh2prior()


    def loglike(self):
        delt = self.theory_.CMBSimpleVec() - self.mean
        return -np.dot(delt, np.dot(self.icov, delt))/2.0



class PLK(SimpleCMBLikelihood):
    def __init__(self, kill_Da= False, kill_rd= False):
        mean = np.array([2.24519776e-02,   1.38572404e-01,   9.43303000e+01])
        cov  = np.array([[1.28572727e-07,  -6.03323687e-07,   1.44305285e-05],
                        [-6.03323687e-07,   7.54205794e-06,  -3.60547663e-05],
                        [1.44305285e-05,  -3.60547663e-05,   4.26414740e-03]])
        name = "SPlanck"
        SimpleCMBLikelihood.__init__(self, name, mean, cov, kill_Da, kill_rd)



                #Calibrated with plikHM_TTTEEE_lowTEB
class PLK15(SimpleCMBLikelihood):
    def __init__(self, kill_Da= False, kill_rd= False):
        mean = np.array([2.24001583e-02,   1.40200580e-01,   9.44043640e+01 ])
        cov = np.array([[3.02751758e-08,  -1.54495460e-07,   4.26868164e-06],
                 [ -1.54495460e-07,   2.16079050e-06,  -1.49955437e-05 ],
                 [ 4.26868164e-06,  -1.49955437e-05,   1.30349464e-03  ]])
        name = "SPlanck_15"
        SimpleCMBLikelihood.__init__(self,name,mean,cov, kill_Da, kill_rd)


#Calibrated with plikHM_TTTEEE_lowTEB
class PLK18(SimpleCMBLikelihood):
    def __init__(self, kill_Da= False, kill_rd= False):
        mean = np.array([2.23619584e-02,   1.425557737e-01,   9.43342292580e+01])
        cov = np.array([[ 2.2280476e-08,  -9.5339119e-08,  -1.5060900e-06],
                 [ -9.5339119e-08,   1.6369370e-06,   1.1826342e-05],
                 [ -1.5060900e-06,   1.1826342e-05,   7.9308203e-04 ]])
        name = "SPlanck_18"
        SimpleCMBLikelihood.__init__(self,name,mean,cov, kill_Da, kill_rd)


class WMAP9(SimpleCMBLikelihood):
    def __init__(self, kill_Da=False, kill_rd=False):
        mean = np.array([2.25946978e-02,   1.35359318e-01,   9.45118918e+01])
        cov  = np.array([[2.86459327e-07,  -4.80929954e-07,  -1.11081266e-05],
                        [-4.80929954e-07,   1.90757225e-05,   7.49542945e-06],
                        [-1.11081266e-05,   7.49542945e-06,   2.54207102e-02]])
        name = "SWMAP"
        SimpleCMBLikelihood.__init__(self, name, mean, cov, kill_Da, kill_rd)



if __name__ == "__main__":
    D = PlanckLikelihood()
    for v in D.mean:
        print("%5.4g" % v)
    for l in D.cov:
        print("%5.4g &  %5.4g &  %5.4g \\" % tuple(l))
