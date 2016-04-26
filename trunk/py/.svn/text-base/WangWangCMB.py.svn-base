##
##
##  !!!!!THIS IS DEPRECATED!!!!!
##  
##  SEE SimpleCMB.py
##
## This module implements Wang+Wang CMB
##

from BaseLikelihood import *
from scipy import *
import scipy.linalg as la

print "DEPRECATED!"
stop()


class WangWangCMB (BaseLikelihood):
    def __init__(self, name, mean, err, cor):
        BaseLikelihood.__init__(self,name)
        self.mean=array(mean)
        ## wang and wang give correlation matrix
        ## and errors and presumably we need to multipy
        ## them back
        err=array(err)
        cov=array(cor)*outer(err,err)
        self.icov = la.inv(cov)

    def setTheory(self,theory):
        self.theory_=theory
        self.theory_.setNoObh2prior()

    def loglike(self):
        delt=self.theory_.WangWangVec()-self.mean
        return -dot(delt,dot(self.icov,delt))/2.0

    

class PlanckLikelihood(WangWangCMB):
    def __init__(self, matrices="PLA3"):
        if matrices=="WW":
            mean=[301.57,1.7407,0.02228]
            err=[0.18,0.0094,0.00030]
            cov=[[ 1.0    , 0.5250  , -0.4235],
                 [ 0.5250 , 1.0    , -0.6925],
                 [-0.4235 , -0.6925,  1.0   ]]
        elif matrices=="PLA1":
            # base_omegak
            mean=[3.01510344e+02,   1.74340204e+00,   2.23128506e-02]
            err=[ 0.18550629,  0.00928648,  0.00031125]
            cov=[[ 1.,          0.57063287, -0.48438085],
                 [ 0.57063287,  1.,         -0.67948316],
                 [-0.48438085, -0.67948316,  1.        ]]
        elif matrices=="PLA2":
            # base_nrun_r_omegak
            mean=[3.01412183e+02, 1.73538851e+00, 2.27716863e-02]
            err=[0.20675659,  0.01203558,  0.00044825]
            cov=[[ 1.,          0.66889481, -0.59554585],
                 [ 0.66889481,  1.,         -0.78180126],
                 [-0.59554585, -0.78180126,  1.        ]]
        elif matrices=="PLA3":
            # base_nrun_r_omegak
            mean= [3.01539311e+02,   1.74066358e+00,   2.25879526e-02]
            err= [ 0.19676676,  0.01056994,  0.00037733]
            cov= [[ 1.        ,  0.62040041, -0.53470863],
                  [ 0.62040041,  1.        , -0.71007217],
                  [-0.53470863, -0.71007217,  1.        ]]





        else:
            print "Basd mtrices param"
            stop()

        WangWangCMB.__init__(self,"CMB_WW_"+matrices,mean,err,cov)
            

class WMAP9Likelihood(WangWangCMB):
    def __init__(self, matrices='PLA'):
        if matrices=="WW":
            mean=[302.02,1.7327,0.02260]
            err=[0.66,0.0164,0.00053]
            cov=[[ 1.0    , 0.3883  , -0.6089],
                 [ 0.3883 , 1.0    , -0.5239],
                 [-0.6089 , -0.5239,  1.0   ]]
        elif matrices=="PLA":
            mean= [ 3.01969812e+02,   1.73081611e+00,   2.26651563e-02]
            err= [  6.61780337e-01,   1.60432146e-02,   5.18589026e-04]
            cov= [[ 1.         , 0.37679434 ,-0.61134328],
                  [ 0.37679434 , 1.         ,-0.51194784],
                  [-0.61134328 ,-0.51194784 , 1.        ]]


        else:
            print "Basd mtrices param"
            stop()


        WangWangCMB.__init__(self,"CMB_WW_WMAP9",mean,err,cov)


        
