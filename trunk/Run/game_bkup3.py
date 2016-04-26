from scipy import *
import random
import scipy.linalg as la
import math as M
import time
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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
        
    def mean_g(self):
        return self.mean

    def converg(self,vec):
	return dot(dot(vec,self.icov),vec)



class Game:
    def __init__ (self, likefuncmany, outfile, par0, sigreg=0.0, like=None, verbose=0):
        random.seed(100)
        self.like=likefuncmany ## returns log like
	self.likeL=like
        self.sigreg=array(sigreg)
        self.N=len(par0)
        self.N1=1000
        self.blow=1.0 ## factor by which to increase the enveloping Gauss
        self.tweight=2.00
        self.wemin=0.00
        self.mineffsamp=self.N1*1
        self.fixedcov=False
        self.toexplore=array(par0)
        self.maxiter=1000
	self.outfile=outfile
	self.verbose=verbose
	if like!=None:
	 self.cpars=like.freeParameters()
	 if (like.name()=="Composite"):
	    self.sublikenames=like.compositeNames()
	    self.composite=True

    def openFiles(self):
	outfile=self.outfile

	fpar=open(outfile+".paramnames",'w')
	for p in self.cpars:
               fpar.write(p.name+"\t\t\t"+p.Ltxname+"\n")
	if self.composite:
           for name in self.sublikenames:
               fpar.write(name+"_like \t\t\t"+name+"\n")
           fpar.write("theory_prior \t\t\t None \n")
        fpar.close()

	
	self.cout=open(outfile+".convergence",'w')
	#self.cout.write('hi there')
	#self.cout.close()
	formn ='%g '*(3) + '\n'
	self.formn=formn

        self.fout=open(outfile+"_1.txt",'w')
        self.mlfout=open(outfile+"_1.maxlike",'w')
	
        formstr = '%g '+'%g '*(self.N+1)
	if (self.composite):
           formstr+='%g '*(len(self.sublikenames)+1)
        formstr+='\n'

	self.formstr=formstr

    def closeFiles(self):
	self.fout.close()
	self.cout.close()

    def run(self):
        comm.Barrier()
        self.initial= time.time()
        toexplore=self.toexplore
	done=False
	self.Gausses=[]
	self.SamList=[]
        if rank==0: self.openFiles()
	weight=0
	mean=0
	cova=0
	b =[]
        i=0
	comm.Barrier()
        while not done:
          self.tsample=time.time()
          sample_list, G=self.isample (toexplore)
          #if rank==0: print '*** end isample', i, time.time()-self.tsample, rank
          #if rank==1: 
          #print '*** end isample', i, time.time()-self.tsample, rank
	  comm.Barrier()
	  if rank==0: 
					#It stores ALL information, check it later
            self.Gausses.append(G)
            self.SamList+=sample_list

            toexplore=self.rebuild_samples(self.SamList, self.Gausses)

	 
                #Convergence
            for sa in self.SamList:
                weight+=sa.we
                mean+=(sa.we*sa.pars)

            a =array(mean)
            b.append(a/weight)

	    #convergt=self.convergt(self.SamList)
	    #mm=1.0
	    #mm2=self.formn%tuple([mm])	
	    #self.cout.write('hi there\n')
	    #print convergt	
            for sa in self.SamList:
                cova+= sa.we*(sa.pars-a)*array([sa.pars-a]).T
		#convp= dot(dot(b[i]-b[i-1],cova)/weight,b[i]-b[i-1])
		#convpp= self.formn%tuple([convp])
		#print convpp
		#self.Cout.write(convpp)		

            #convergence= dot(dot(b[i]-b[i-1],cova)/weight,b[i]-b[i-1])
	    #mm=self.formn%tuple([convergence]) 
	    #self.cout.write(mm)
            if i>0:
		  convergence= dot(dot(b[i]-b[i-1],cova)/weight,b[i]-b[i-1])
	    	  cvg=self.formn%tuple([i, self.N1*i, convergence])
	    	  self.cout.write(cvg)	
            #   if Convergence < self.N/1000.0:
                  print '---Convergence= ', convergence, ', #params/100= ', self.N/100.0
                  print 'for Gaussians= ', i, ' and Nsamples= ', self.N1
           #       done=True
                

            if (self.wemax<self.tweight):
                done=True
            if (len(self.Gausses)>=self.maxiter):
                print "Max iter exceeded"
                done=True
            if (self.effsamp<self.mineffsamp):
                done=False


            if (i==50):
                done=True
	    #if i>0:
            #   vect = self.Gausses[i].mean_g() - self.Gausses[i-1].mean_g()
            #   print '---------------------', self.Gausses[i].converg(vect)
          else: 
            toexplore=[]
            G=None
	  
	  i+=1
	  comm.Barrier()
          #print '***end gaussian', i, time.time()-self.tsample, rank
          #comm.Barrier()
        print '***--end run', time.time()-self.initial, rank

        if rank==0: self.closeFiles()
        comm.Barrier()     
    #def convergt(self, SamList):
	#self.cout.write('hi there\n')
	#weight=0
        #mean=0
        #cova=0
        #b =[]
        #i=0
	#for sa in self.SamList:
	#    weight+=sa.we
	#    mean+=(sa.we*sa.pars)
	#a =array(mean)
	#b.append(a/weight)

	#for sa in self.SamList:
	#    cova+= sa.we*(sa.pars-a)*array([sa.pars-a]).T

	#convergence= dot(dot(b[i]-b[i-1],cova)/weight,b[i]-b[i-1])
	#print convergence
	#return convergence


    def gausses_eval(self,sam):
        if len(sam.glikes)!=len(self.Gausses):
            stop("SHIT")
        probi=(exp(array(sam.glikes))).sum()
        return probi

    def rebuild_samples(self, SamList,Gausses):
        maxlike=-1e30
        gmaxlike=-1e30
        j=0
	params=self.cpars

        for sa in SamList:
            if (sa.like>maxlike):
                maxlike=sa.like
                maxlikesa=sa
            sa.glike=self.gausses_eval(sa) 
            if (sa.glike>gmaxlike):
                gmaxlike=sa.glike
                
        gmaxlike2=self.gausses_eval(maxlikesa)
        if(self.verbose>0):
	   print 'gmaxlike, gmaxlike2' , gmaxlike, gmaxlike2
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

	    #print sa.we, j
  	    for i in range(len(params)):
            	params[i].setValue(sa.pars[i])

	    self.likeL.updateParams(params) 	    
	    self.ploglikes=self.getLikes()
	    vect=[sa.pars[i] for i in range(0,len(sa.pars))]

	    #print values
	    if (self.composite):
 	       outstr=self.formstr%tuple([sa.we,rellike]+vect+self.ploglikes.tolist())
	    else:
	       outstr=self.formstr%tuple([sa.we,rellike]+vect)

	    if(not M.isnan(sa.we)):
		 self.fout.write(outstr)

	    if (rellike>maxlike):	
		maxlike=rellike
		#print "New maxloglike", maxlike
		self.mlfout.seek(0)
		self.mlfout.write(outstr)
		self.mlfout.flush()

	    j+=1
        self.sample_list=flist
        #print "#G=",len(Gausses), "maxlike=",maxlike,"wemax=",wemax,"effsamp=",effsamp
        #print 'time so far= ', time.time()-self.initial
        self.effsamp=effsamp
        self.wemax=wemax
        return parmaxw

    def getLikes(self):
	if (self.composite):
  	   cloglikes=self.likeL.compositeLogLikes_wprior()
	else:
	   cloglikes = []

	return cloglikes
                        
    def getcov(self, around):

      toget=[]
      if rank==0:
         toget.append(around)
         N=self.N

         if (self.fixedcov):
            cov=zeros((N,N))
            for i in range(N):
                cov[i,i]=self.sigreg[i]**2
            if(self.verbose>0):
	       print cov
            G=Gaussian(around,cov)    
            return G

         icov=zeros((N,N))
         delta=self.sigreg/1000.0

        
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

      else: toget=[]
      comm.Barrier()
	
      likes=self.like(toget, Npar=True)
      comm.Barrier()    

      if rank==0:
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
            if(self.verbose>0):
	       print "Regularizing cholesky"
            for i in range(N):
                icov[i,i]+=1/self.sigreg[i]**2
            try:
                ch=la.cholesky(icov)
                break
            except:
                pass

         cov=la.inv(icov)
         if(self.verbose>0):
	   print cov
         G=Gaussian(around,self.blow*cov)    

      else: G=None 
      return G


    def isample (self, zeropar):
       
      ## Get local covariance matrix 
      init2=time.time()
      if True: #rank==0:
         G=self.getcov(zeropar) 
      #else:
      #   G=[]
      comm.Barrier()
      #print '**--- getcov',time.time()-init2, rank
      slist=[]
      lmany=[]
      likes=[]      
 
      if rank==0:
        ### first update the old samples
        for i,s in enumerate(self.SamList):
            self.SamList[i].glikes.append(G.like(s.pars))
            
        ## now sample around this 
        for i in range(self.N1):
            par,glike=G.sample()
            glikel=[g.like(par) for g in self.Gausses] + [glike]
            lmany.append(par)
            like=None
            slist.append(Sample(par,like, glikel))
      else:
        slist=[]
        lmany=[]
        likes=[]
 
      comm.Barrier()
      init=time.time()
      likes=self.like(lmany, Npar=True)
      #print '*** time parallel', time.time()-init, rank 
      
      if rank==0:
        for like,sa in zip(likes,slist):
            sa.like=like
      else:
        for like,sa in zip(likes,slist):
            sa.like=None

      comm.Barrier()
      return slist,G
