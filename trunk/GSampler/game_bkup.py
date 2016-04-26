from scipy import *
import random
import scipy.linalg as la

class Sample:
    def __init__ (self,pars, like, glikes):
        self.pars=pars
        self.like=like
        self.glikes=glikes
        
class Gaussian:
    def __init__(self,mean,cov):
        self.cov=cov
        self.mean=mean
        self.chol=la.cholesky(cov)
        self.lndet=log(self.chol.diagonal()).sum()*2.0
        self.icov=la.inv(cov)
        self.N=len(cov)

    def sample(self):
        da=array([random.gauss(0.,1.) for x in range(self.N)])
        glike = -(da**2).sum()/2.0-self.lndet/2.0
        sa=dot(da,self.chol)
        if (self.mean!=None):
            sa+=self.mean
        return sa,glike
    
    def chi2(self,vec):
        if mean!=None:
            delta=vec-self.mean
        else:
            delta=vec
        return dot(dot(delta,self.icov),delta)
        
    def like(self,vec):
        return -self.chi2(vec)/2-self.lndet/2.0
        

class Game:
    def __init__ (self, likefuncmany, par0, sigreg=0.0):
        random.seed(10)
        self.like=likefuncmany ## returns log like
        self.sigreg=array(sigreg)
        self.N=len(par0)
        self.N1=1000
        self.blow=2.0 ## factor by which to increase the enveloping Gauss
        self.tweight=2.00
        self.wemin=0.00
        self.mineffsamp=self.N1*1
        self.fixedcov=False
        self.toexplore=array(par0)
        self.maxiter=30

    def run(self):
        done=False
        toexplore=self.toexplore
        badlist=[]
        self.Gausses=[]
        self.SamList=[]
        while not done:
            sample_list, G=self.isample (toexplore)
            self.Gausses.append(G)
            self.SamList+=sample_list

            toexplore=self.rebuild_samples(self.SamList, self.Gausses)
            
            if (self.wemax<self.tweight):
                done=True
            if (len(self.Gausses)>=self.maxiter):
                print "Max iter exceeded"
                done=True
            if (self.effsamp<self.mineffsamp):
                done=False
            
    def gausses_eval(self,sam):
        if len(sam.glikes)!=len(self.Gausses):
            stop("SHIT")
        probi=(exp(array(sam.glikes))).sum()
        return probi

    def rebuild_samples(self, SamList,Gausses):
        maxlike=-1e30
        gmaxlike=-1e30
        for sa in SamList:
            if (sa.like>maxlike):
                maxlike=sa.like
                maxlikesa=sa
            sa.glike=self.gausses_eval(sa) 
            if (sa.glike>gmaxlike):
                gmaxlike=sa.glike
                
        gmaxlike2=self.gausses_eval(maxlikesa)
        print gmaxlike, gmaxlike2,'AA'
        #gmaxlike=gmaxlike2
        wemax=0.0
        flist=[]
        wemax=0.0
        parmaxw=None
        effsamp=0
        for sa in SamList:
            rellike=exp(sa.like-maxlike)
            glike=sa.glike/gmaxlike
            we=rellike/glike
            sa.we=we
            effsamp+=min(1.0,we)
            if we>wemax:
                wemax=we
                parmaxw=sa.pars
            if we>self.wemin:
                flist.append(sa)
                
        self.sample_list=flist
        print "#G=",len(Gausses), "maxlike=",maxlike,"wemax=",wemax,"effsamp=",effsamp
        self.effsamp=effsamp
        self.wemax=wemax
        return parmaxw

                        
    def getcov(self, around):
        N=self.N

        if (self.fixedcov):
            cov=zeros((N,N))
            for i in range(N):
                cov[i,i]=self.sigreg[i]**2
            print cov
            G=Gaussian(around,cov)    
            return G

        icov=zeros((N,N))
        delta=self.sigreg/1000.0
        toget=[]
        toget.append(around)
        
        ### This is a kinda ugly hack
        ### We repeat the exactly the same loop twice.
        ### first populating where to evaluate like 
        ### and the popping hoping for perfect sync


        for i in range(N):
            parspi=around*1.0
            parsmi=around*1.0
            parspi[i]+=delta[i]
            parsmi[i]-=delta[i]
            for j in range(N):
                if (i==j):
                    toget.append(parspi)
                    toget.append(parsmi)
                else:
                    parspp=parspi*1.0
                    parspm=parspi*1.0
                    parsmp=parsmi*1.0
                    parsmm=parsmi*1.0
                    parspp[j]+=delta[j]
                    parspm[j]-=delta[j]
                    parsmp[j]+=delta[j]
                    parsmm[j]-=delta[j]
                    toget.append(parspp)
                    toget.append(parsmm)
                    toget.append(parspm)
                    toget.append(parsmp)

        likes=self.like(toget)

        like0=likes.pop(0)
        for i in range(N):
            for j in range(N):
                if (i==j):
                    der=(likes.pop(0)+likes.pop(0)-2*like0)/(delta[i]**2)
                else:
                    der=(likes.pop(0)+likes.pop(0)-likes.pop(0)-likes.pop(0))/(4*delta[i]*delta[j])
                icov[i,j]=-der
                icov[j,i]=-der

        while True:
            print "Regularizing cholesky"
            for i in range(N):
                icov[i,i]+=1/self.sigreg[i]**2
            try:
                ch=la.cholesky(icov)
                break
            except:
                pass

        cov=la.inv(icov)
        print cov
        G=Gaussian(around,self.blow*cov)    
        return G


    def isample (self, zeropar):
        
        ## Get local covariance matrix
        G=self.getcov(zeropar)
        
        ### first update the old samples
        for i,s in enumerate(self.SamList):
            self.SamList[i].glikes.append(G.like(s.pars))
            
        slist=[]
        lmany=[]
        ## now sample around this 
        for i in range(self.N1):
            par,glike=G.sample()
            glikel=[g.like(par) for g in self.Gausses] + [glike]
            lmany.append(par)
            like=None
            slist.append(Sample(par,like, glikel))
        likes=self.like(lmany)
        for like,sa in zip(likes,slist):
            sa.like=like

        return slist,G
