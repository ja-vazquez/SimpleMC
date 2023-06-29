

from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import scipy.linalg as la
import scipy as sp
from simplemc import cdir
import numpy as np

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
        self.errHz = da[:,2]
        print("Loading ", cov_filename)
        cov = sp.loadtxt(cov_filename, skiprows=1)
        assert(len(cov) == len(self.zs))
        if len(cov) == 15:
            filename = cdir+'/data/data_MM20.dat'
            zmod, imf, slib, sps, spsooo = np.genfromtxt(filename, comments='#', usecols=(0, 1, 2, 3, 4), unpack=True)
            cov_mat_diag = np.zeros((len(self.zs), len(self.zs)), dtype='float64')

            for i in range(len(self.zs)):
                cov_mat_diag[i,i] = self.errHz[i]**2

            imf_intp = np.interp(self.zs, zmod, imf) / 100
            slib_intp = np.interp(self.zs, zmod, slib) / 100
            sps_intp = np.interp(self.zs, zmod, sps) / 100
            spsooo_intp = np.interp(self.zs, zmod, spsooo) / 100

            cov_mat_imf = np.zeros((len(self.zs), len(self.zs)), dtype='float64')
            cov_mat_slib = np.zeros((len(self.zs), len(self.zs)), dtype='float64')
            cov_mat_sps = np.zeros((len(self.zs), len(self.zs)), dtype='float64')
            cov_mat_spsooo = np.zeros((len(self.zs), len(self.zs)), dtype='float64')

            for i in range(len(self.zs)):
                for j in range(len(self.zs)):
                    cov_mat_imf[i, j] = self.Hs[i] * imf_intp[i] * self.Hs[j] * imf_intp[j]
                    cov_mat_slib[i, j] = self.Hs[i] * slib_intp[i] * self.Hs[j] * slib_intp[j]
                    cov_mat_sps[i, j] = self.Hs[i] * sps_intp[i] * self.Hs[j] * sps_intp[j]
                    cov_mat_spsooo[i, j] = self.Hs[i] * spsooo_intp[i] * self.Hs[j] * spsooo_intp[j]

            cov = cov_mat_spsooo + cov_mat_imf + cov_mat_diag
        else:
            cov += 3 ** 2

        vals, vecs = la.eig(cov)
        vals = sorted(sp.real(vals))
        print("Eigenvalues of cov matrix:", vals[0:3],'...',vals[-1])
        print("Adding marginalising constant")

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
        CompressedHDLikelihood.__init__(self,"HD", cdir+"/data/HDiagramCompilacion-data_31.txt",
                                            cdir+"/data/HDiagramCompilacion-cov_31.txt")

class HD23(CompressedHDLikelihood):
    # Updated data (2020) from https://arxiv.org/abs/2003.07362
    # with full covariance matrix
    def __init__(self):
        CompressedHDLikelihood.__init__(self,"HD23", cdir+"/data/HzTable_MM_BC03.dat",
                                            cdir+"/data/cov_cc_MM_BC03.txt")
