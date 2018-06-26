#
# This module calculates likelihood for Hubble Diagram.
#

from BaseLikelihood import *
from scipy import *
import scipy.linalg as la

class CompressedHDLikelihood (BaseLikelihood):
    def __init__(self,name,values_filename, cov_filename):
        BaseLikelihood.__init__(self,name)
        print "Loading ",values_filename
        da=loadtxt(values_filename)
        self.zs=da[:,0]
        self.Hs=da[:,1]
        print "Loading ",cov_filename
        cov=loadtxt(cov_filename,skiprows=1)
        assert(len(cov)==len(self.zs))
        vals, vecs=la.eig(cov)
        vals=sorted(real(vals))
        print "Eigenvalues of cov matrix:", vals[0:3],'...',vals[-1]
        print "Adding marginalising constant"
        cov+=3**2
        vals, vecs=la.eig(cov)
        vals=sorted(real(vals))
        print "Eigenvalues of cov matrix:", vals[0:3],'...',vals[-1]
        self.icov=la.inv(cov)


    def loglike(self):
        tvec=array([self.theory_.RHSquared_a(1.0/(1+z)) for z in self.zs])
	print '*****',tvec, self.Hs
	## This is the factor that we need to correct
        ## note that in principle this shouldn't matter too much, we will marginalise over this
        tvec+=70
        delta=tvec-self.Hs
        return -dot(delta,dot(self.icov,delta))/2.0


class HD(CompressedHDLikelihood):
    def __init__(self):
        CompressedHDLikelihood.__init__(self,"HD","data/HDiagramCompilacion-data.txt",
                                             "data/HDiagramCompilacion-cov.txt")
