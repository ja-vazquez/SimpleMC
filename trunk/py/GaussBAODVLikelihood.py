#
# This is a BAO likelihood, where we give DV and its error and optionally a value
# at which chi2 is cut.
#

from BaseLikelihood import *
from scipy import *
from scipy.interpolate import RectBivariateSpline

class GaussBAODVLikelihood (BaseLikelihood):
    def __init__(self,name,z,DV,DVErr, fidtheory, maxchi2=1e30):
        BaseLikelihood.__init__(self,name)
        self.z=z
        rd=fidtheory.rd
        DV/=rd
        DVErr/=rd
        self.maxchi2=maxchi2
        print name,"measurement in ", fidtheory.rd_approx,":",DV,"+-",DVErr
        self.setData(DV,DVErr)

    def setData(self,DV,DVErr):
        self.DV=DV
        self.DVErr2=DVErr**2
        

    def loglike(self):
        DVT=self.theory_.DVOverrd(self.z)
        chi2=min(self.maxchi2,(DVT-self.DV)**2/(self.DVErr2))
        return -chi2/2.0


    
