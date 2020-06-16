

from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import scipy.linalg as la
import scipy as sp

class CompressedHDLikelihood(BaseLikelihood):
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
        self.zs = da[:,0]
        self.Hs = da[:,1]
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
        tvec = sp.array([100.0*self.theory_.h*sp.sqrt(self.theory_.RHSquared_a(1.0/(1+z))) for z in self.zs])
        #print tvec, self.Hs
        ## This is the factor that we need to correct
        ## note that in principle this shouldn't matter too much, we will marginalise over this
        tvec += 0
        delta = tvec - self.Hs
        return -sp.dot(delta, sp.dot(self.icov, delta))/2.0


class HubbleDiagram(CompressedHDLikelihood):
    # data from https://arxiv.org/abs/1802.01505
    def __init__(self):
        CompressedHDLikelihood.__init__(self,"HD","simplemc/data/HDiagramCompilacion-data_31.txt",
                                             "simplemc/data/HDiagramCompilacion-cov_31.txt")
