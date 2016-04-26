##
## This is inspired by WangWang CMB, but ultimately simpler
## where we compress to just obh2, obh2+och2 and thetas
## 

from BaseLikelihood import *
from scipy import *
import scipy.linalg as la

class SimpleCMB (BaseLikelihood):
    def __init__(self, name, mean, cov, kill_Da=False, kill_rd=False):
        print "Initializing CMB likelihood:",name
        cov=array(cov)
        if kill_Da:
            cov[2,:]=0.0
            cov[:,2]=0.0
            cov[2,2]=1e10
            name+="_noDa"
        if kill_rd:
            cov[0:2,0:2]=0.0
            cov[0,0]=1e10
            cov[1,1]=1e10
            name+="_noRd"
        BaseLikelihood.__init__(self,name)
        self.mean=array(mean)
        self.cov=cov
        self.icov = la.inv(cov)

    def setTheory(self,theory):
        self.theory_=theory
        self.theory_.setNoObh2prior()

    def loglike(self):
        delt=self.theory_.CMBSimpleVec()-self.mean
        return -dot(delt,dot(self.icov,delt))/2.0

    

class PlanckLikelihood(SimpleCMB):
    def __init__(self, kill_Da=False, kill_rd=False):
        mean = array([   2.24519776e-02,   1.38572404e-01,   9.43303000e+01])
        cov = array([[  1.28572727e-07,  -6.03323687e-07,   1.44305285e-05 ],
               [       -6.03323687e-07,   7.54205794e-06,  -3.60547663e-05 ],
               [        1.44305285e-05,  -3.60547663e-05,   4.26414740e-03]])
        name="SPlanck"
        SimpleCMB.__init__(self,name,mean,cov, kill_Da, kill_rd)



class WMAP9Likelihood(SimpleCMB):
    def __init__(self, kill_Da=False, kill_rd=False):
        mean = array([   2.25946978e-02,   1.35359318e-01,   9.45118918e+01])
        cov = array([[  2.86459327e-07,  -4.80929954e-07,  -1.11081266e-05],
                     [ -4.80929954e-07,   1.90757225e-05,   7.49542945e-06],
                     [ -1.11081266e-05,   7.49542945e-06,   2.54207102e-02]])
        name="SWMAP"
        SimpleCMB.__init__(self,name,mean,cov, kill_Da, kill_rd)


if __name__=="__main__":
    D=PlanckLikelihood()
    for v in D.mean:
        print "%5.4g"%v
    for l in D.cov:
        print "%5.4g &  %5.4g &  %5.4g \\"%tuple(l)

    
