#!/usr/bin/env python
from RunBase import *
from game import *
from scipy import *
import pylab, sys
import time
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if (len(sys.argv)<4):
    print """Usage: 
Run/driver_game.py pre/phy model dataset
where model is one of 
%s 
and datasets is one of 
%s
Example:
Run/driver.py phy LCDM BBAO+CMBP 
"""%(model_list, data_list)
    sys.exit(1)

prefact,model,datasets=sys.argv[1:4]
chainsdir = 'chains/'
random.seed(100)
verbose =2
print "Running ", model, prefact,datasets,'with processors:',size 

###---------------------------Set up Theory and Likelihood ------------

#print rank
if rank==0:
  T=ParseModel(model)
  L=ParseDataset(datasets)
  if prefact=="pre":
     T.setVaryPrefactor()
  T.printFreeParameters()
  Ll=L.setTheory(T)
else:
  T=None 
  Ll=None

T =comm.bcast(T, root=0)
L =comm.bcast(Ll, root=0)
if rank==1: T.printFreeParameters()

params=L.freeParameters()
values= [p.value for p in params]
errors= [p.error for p in params]
 

###----------------------------------------Game Sampler------------
if False:
  def negloglike(x):
      for i in range(len(params)):
          params[i].setValue(x[i])
      L.updateParams(params)
      loglike=L.loglike_wprior()
      return loglike


  def likemany(x):	
      if rank==0:
        chunks = [[] for _ in range(size)]
        n= ceil(float(len(x))/size)
        for i, chunk in enumerate(x):
            chunks[int(i//n)].append(chunk)
      else:
        y = None
        chunks = None

      y=comm.scatter(chunks, root=0)
      result2=map(negloglike,y)
      comm.Barrier()
      z=comm.gather(result2, root=0)
      if rank==0:
          zr=reduce(lambda x,y: x+y,z)
      else:
          zr=None
      result=comm.bcast(zr, root=0)	
      return result

  #params =comm.bcast(params, root=0)
  #values =comm.bcast(values, root=0)
  #errors =comm.bcast(errors, root=0)

  name = chainsdir + "/" + model+"_"+prefact+"_"+datasets 
  ga=Game(likemany,name , values, errors, like=L, verbose=0)
  ga.N1=5
  ga.tweight=1.5
  ga.mineffsamp=5000
  ga.run()



### Useful now only for simple tests, i.e. ring, gauss...
###----------------------------------------Plotting ------------

if False:
 sname='gauss.pdf'
 xx=array([sa.pars[0] for sa in ga.sample_list])
 yy=array([sa.pars[2] for sa in ga.sample_list])
 ww=array([sa.we for sa in ga.sample_list])

 ## now we plot
 Np=100
 cxmin=0.15 
 cxmax=0.65
 cymin=0.5
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

 pylab.subplot(2,2,1)
 pylab.imshow(sums, interpolation='nearest', origin='lower left',extent=[cxmin,cxmax,cymin,cymax], aspect='auto')
 pylab.xlim(cxmin,cxmax)
 pylab.ylim(cymin,cymax)

 pylab.subplot(2,2,2)
 pylab.imshow(wsumsa, interpolation='nearest', origin='lower left',extent=[cxmin,cxmax,cymin,cymax],vmin=0, vmax=vmax, aspect='auto')

 pylab.subplot(2,2,3)
 pylab.imshow(trvalsa, interpolation='nearest', origin='lower left',extent=[cxmin,cxmax,cymin,cymax],vmin=0, vmax=vmax, aspect='auto')

 pylab.subplot(2,2,4)
 pylab.imshow(diffp, interpolation='nearest', origin='lower left',extent=[cxmin,cxmax,cymin,cymax],aspect='auto')
 pylab.colorbar()


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

