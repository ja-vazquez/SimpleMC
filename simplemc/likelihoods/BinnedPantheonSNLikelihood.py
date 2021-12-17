
from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d
from simplemc import cdir


class BinnedPantheonSNLikelihood(BaseLikelihood):
    def __init__ (self, ninterp=150):
        """
        This module calculates likelihood for the binned version of Pantheon Compilation.
        Data from https://github.com/dscolnic/Pantheon and arXiv:1710.00845.

        Parameters
        ----------
        ninterp

        Returns
        -------

        """
        ## first read data file
        self.name_="BinnedPantheonSN"
        da= np.loadtxt(cdir+'/data/binned_pantheon.txt', skiprows=1)
        self.zcmb = da[:, 1]
        self.zhelio = da[:, 2]
        self.mag = da[:, 4]
        self.dmag = da[:, 5]
        self.N=len(self.mag)
        self.syscov=np.loadtxt(cdir+'/data/binned_cov_pantheon.txt',skiprows=1).reshape((self.N,self.N))
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
        dist=interp1d(self.zinter,[self.theory_.distance_modulus(z) for z in self.zinter],
                       kind='cubic',bounds_error=False)(self.zcmb)
        who=np.where(self.zcmb>self.zmaxi)
        dist[who]=np.array([self.theory_.distance_modulus(z) for z in self.zcmb[who]])
        tvec = self.mag-dist

        #tvec = self.mag-np.array([self.theory_.distance_modulus(z) for z in self.zcmb])
        #print (tvec[:10])
        ## first subtract a rought constant to stabilize marginaliztion of
        ## intrinsic mag.
        tvec-= (tvec*self.xdiag).sum() / (self.xdiag.sum())
        #print (tvec[:10])
        chi2 = np.einsum('i,ij,j',tvec,self.icov,tvec)
        #print ("chi2=",chi2)
        return -chi2/2
        
    
