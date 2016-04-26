#!/usr/bin/env python
from RunBase import *
from game import *
from scipy import *
import pylab, sys
import time
from mpi4py import MPI

initime=time.time()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


if (len(sys.argv)<4):
    print """Usage: 
Run/driver_game.py pre/phy model dataset [#Samples]
where model is one of 
%s 
and datasets is one of 
%s
Example:
Run/driver.py phy LCDM BBAO+CMBP 
"""%(model_list, data_list)
    sys.exit(1)

prefact,model,datasets=sys.argv[1:4]
if len(sys.argv)>4:
    Nsampling=int(sys.argv[4])
else:
    Nsampling=100


verbose =2
random.seed(100)
chainsdir = 'chains/'
name = chainsdir + "/" + model+"_"+prefact+"_"+datasets+"_GS_"+str(Nsampling)
print "Running ", name
print 'with processors:',size, 'in the rank, ', rank 

###---------------------------Set up Theory and Likelihood ------------

			#load likelihood and datasets
if True:
  T=ParseModel(model)
  L=ParseDataset(datasets)
  if prefact=="pre":
     T.setVaryPrefactor()
  T.printFreeParameters()
  L.setTheory(T)

  params= L.freeParameters()
  values= [p.value for p in params]
  errors= [p.error for p in params]


###----------------------------------------Game Sampler------------
if True:
			#Likelihood
  def negloglike(x):
      for i in range(len(params)):
          params[i].setValue(x[i])
      L.updateParams(params)
      loglike=L.loglike_wprior()	
      return loglike

			#Parallelization for the Lik evaluations
  def likemany(x, Npar=True):	
        comm.Barrier()  
        
        if rank==0:
          chunks = [[] for _ in range(size)]
          n= ceil(float(len(x))/size)
          for i, chunk in enumerate(x):
              chunks[int(i//n)].append(chunk)
        else:
          chunks = []
          
        y2=comm.scatter(chunks, root=0)
        result2=map(negloglike,y2)
        comm.Barrier()
        z=comm.gather(result2, root=0)

        if rank==0:
            zr=reduce(lambda xi,yi: xi+yi,z)
        else:
            zr=None
        result= comm.bcast(zr, root=0)	
        
        comm.Barrier()
        return result



  comm.Barrier() 
  ga=Game(likemany, name , values, errors, like=L, verbose=0)
  ga.N1=Nsampling
  ga.tweight=1.5
  ga.mineffsamp=5000
  
  ga.run()
  
  print 'total time =',time.time()-initime
  comm.Abort()


### Useful for simple tests, i.e. ring, gauss...
###----------------------------------------Plotting ------------

if False:
 print 'making plots'

 def plotel(G):
    cov=G.cov
    print G.cov
    val,vec=linalg.eig(cov)
    vec=vec.T

    vec[0]*=sqrt(real(val[0]))
    vec[2]*=sqrt(real(val[2]))
    print vec[0],'A'
    print vec[1],'B'
    pylab.plot(G.mean[0],G.mean[2],'bo')
    pylab.plot([G.mean[0]-vec[0][0],G.mean[0]+vec[0][0]],
               [G.mean[2]-vec[0][2],G.mean[2]+vec[0][2]],'r-')

    print 'vectors', val, vec[2][0], vec[2][2], vec[0][2], vec[0][0]
    pylab.plot([G.mean[0]-vec[2][0],G.mean[0]+vec[2][0]],
               [G.mean[2]-vec[2][2],G.mean[2]+vec[2][2]],'r-')


 sname='gauss.pdf'
 xx=array([sa.pars[0] for sa in ga.sample_list])
 yy=array([sa.pars[2] for sa in ga.sample_list])
 ww=array([sa.we for sa in ga.sample_list])

 ## now we plot
 Np=50
 cxmin=0.0
 cxmax=0.9
 cymin=0.4
 cymax=0.9	

 cxstep=(cxmax-cxmin)/(1.0*Np)
 cystep=(cymax-cymin)/(1.0*Np)

 sums=zeros((Np,Np))
 wsums=zeros((Np,Np))
 trvals=zeros((Np,Np))

 for x,y,w in zip(xx,yy,ww):
    if (x<cxmin) or (x>cxmax) or (y<cymin) or (y>cymax):
        continue
    ix=int((x-cxmin)/cxstep)
    iy=int((y-cymin)/cystep)
    sums[iy,ix]+=1.0
    wsums[iy,ix]+=w

 def like(x):
   return negloglike(x)

 for i in range(Np):
    x=cxmin+(i+0.5)*cxstep
    for j in range(Np):
        y=cymin+(j+0.5)*cystep
	#trvals[j,i]=exp(like([x,y]))
        trvals[j,i]=exp(like([x,0.02234,y]))

 trvalsa=trvals/trvals.sum()
 wsumsa=wsums/wsums.sum()
 diffp=wsumsa-trvalsa
 vmax=trvalsa.max()*1.1

 pylab.figure(figsize=(8,8))

 #pylab.subplot(2,2,1)
 pylab.imshow(sums, interpolation='nearest', origin='lower left',extent=[cxmin,cxmax,cymin,cymax], aspect='auto')
 for G in ga.Gausses:
    plotel(G)
 pylab.xlim(cxmin,cxmax)
 pylab.ylim(cymin,cymax)

 #pylab.subplot(2,2,2)
 #pylab.imshow(wsumsa, interpolation='nearest', origin='lower left',extent=[cxmin,cxmax,cymin,cymax],vmin=0, vmax=vmax, aspect='auto')

 #pylab.subplot(2,2,3)
 #pylab.imshow(trvalsa, interpolation='nearest', origin='lower left',extent=[cxmin,cxmax,cymin,cymax],vmin=0, vmax=vmax, aspect='auto')

 #pylab.subplot(2,2,4)
 #pylab.imshow(diffp, interpolation='nearest', origin='lower left',extent=[cxmin,cxmax,cymin,cymax],aspect='auto')
 #pylab.colorbar()


 mx=(xx*ww).sum()/(ww.sum())
 vx=sqrt((xx**2*ww).sum()/(ww.sum())-mx*mx)
 my=(yy*ww).sum()/(ww.sum())
 vy=sqrt((yy**2*ww).sum()/(ww.sum())-my*my)
 rr=xx**2+yy**2
 mr=(rr*ww).sum()/(ww.sum())
 vr=sqrt((rr**2*ww).sum()/(ww.sum())-mr*mr)

 print 'xmean,xvar=',mx,vx
 print 'ymean,yvar=',my,vy
 print 'rmean,rvar=',mr,vr

 pylab.tight_layout()
 pylab.savefig(sname)
 pylab.show()

