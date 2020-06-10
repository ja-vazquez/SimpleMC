#
# This module calculates likelihood for the full SN Pantheon
#

from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
from scipy.interpolate import interp1d
import scipy.linalg as la
import numpy as np


class PantheonLikelihood(BaseLikelihood):
    def __init__ (self, ninterp=150):
        ## first read data file
        self.name_="PantheonSN"
        da=[x.split() for x in open('simplemc/data/pantheon_lcparam_full_long_zhel.txt').readlines()[1:]]
        self.zcmb = np.array([float(line[1]) for line in da])
        self.zhelio = np.array([float(line[2]) for line in da])
        self.mag = np.array([float(line[4]) for line in da])
        self.dmag = np.array([float(line[5]) for line in da])
        self.N=len(self.mag)
        self.syscov=np.loadtxt('simplemc/data/pantheon_sys_full_long.txt',skiprows=1).reshape((self.N,self.N))
        self.cov=np.copy(self.syscov)
        self.cov[np.diag_indices_from(self.cov)]+=self.dmag**2
        self.xdiag=1/self.cov.diagonal() ## diagonal before marginalising constant
        ## add marginalising over a constant
        self.cov+=3**2
        self.zmin=self.zcmb.min()
        self.zmax=self.zcmb.max()
        self.zmaxi=1.1 ## we interpolate to 1.1 beyond that exact calc
        print ("Pantheon SN: zmin=%f zmax=%f N=%i"%(self.zmin,self.zmax, self.N))
        self.zinter=np.linspace(1e-3,self.zmaxi,ninterp) 
        self.icov=la.inv(self.cov)


    def loglike(self):
        ## we will interpolate distance
        dist=interp1d(self.zinter,[self.theory_.genericPModel(z) for z in self.zinter],
                       kind='cubic',bounds_error=False)(self.zcmb)

        who=np.where(self.zcmb>self.zmaxi)
        dist[who]=np.array([self.theory_.genericPModel(z) for z in self.zcmb[who]])

        tvec = self.mag-dist
        tvec-= (tvec*self.xdiag).sum() / (self.xdiag.sum())

        chi2 = np.einsum('i,ij,j',tvec,self.icov,tvec)

        return -chi2/2

    
    
        
        
        
    
