

from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import scipy.linalg as la
import scipy as sp


class Compressedfs8Likelihood(BaseLikelihood):
    def __init__(self,name,values_filename, cov_filename):
        """
        This module calculates likelihood for Hubble Diagram.
        based on the CompressSN file
        Parameters
        ----------
        name
        values_filename
        cov_filename

        Returns
        -------

        """
        BaseLikelihood.__init__(self,name)
        print("Loading ",values_filename)
        da = sp.loadtxt(values_filename)
        self.zs  = da[:,0]
        self.fs8 = da[:,1]
        print("Loading ",cov_filename)
        cov = sp.loadtxt(cov_filename,skiprows=1)
        assert(len(cov) == len(self.zs))
        vals, vecs = la.eig(cov)
        vals = sorted(sp.real(vals))
        print("Eigenvalues of cov matrix:", vals[0:3],'...',vals[-1])
        print("Adding marginalising constant")
        cov += 3**2
        self.icov = la.inv(cov)


    def loglike(self):
        tvec = sp.array([self.theory_.fs8(z) for z in self.zs])
        delta = tvec - self.fs8
        return -sp.dot(delta, sp.dot(self.icov, delta))/2.0



class fs8Diagram(Compressedfs8Likelihood):
    # data from arXiv:1806.10822
    def __init__(self):
        Compressedfs8Likelihood.__init__(self,"fs8", "simplemc/data/fs8Diagram.txt",
                                             "simplemc/data/fs8Diagram-cov.txt")
