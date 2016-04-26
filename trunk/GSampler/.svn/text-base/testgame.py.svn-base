#!/usr/bin/env python
from game import *
from scipy import *
import pylab, sys
import scipy.linalg as linalg

random.seed(100)
def likemany(x):
    return map(like,x)

if sys.argv[1]=='gauss':
    def like(x):
        return -((x[0])**2+(x[1])**2/1.0- 0.0*x[0]*x[1])/2.0
    ga=Game(likemany,[0.0,0.0],[0.5,0.5])
    ga.N1=1000
    ga.tweight=1.50
    ga.mineffsamp=5000
    sname='gauss.pdf'
    ga.run()


elif sys.argv[1]=='ring':
    def like(x):
        r2=x[0]**2+x[1]**2
        return -(r2-4.0)**2/(2*0.5**2)
    ga=Game(likemany,[3.5,0.0],[0.5,0.9])
    ga.blow=2.0
    ga.tweight=1.50
    sname='ring.pdf'
    ga.run()

elif sys.argv[1]=='box':
    def like(x):
        if (abs(x[0])>1) or (abs(x[1])>1):
            return -30
        else:
            return 0
    ga=Game(likemany,[0.5,0.0],[0.4,0.4])
    ga.tweight=1.5
    ga.N1=1000
    ga.run()
    sname='box.pdf'
else:
    stop ("define")

def plotel(G):
    cov=G.cov
    print G.cov
    val,vec=linalg.eig(cov)
    vec=vec.T

    vec[0]*=sqrt(real(val[0]))
    vec[1]*=sqrt(real(val[1]))
    print vec[0],'A'
    print vec[1],'B'
    pylab.plot(G.mean[0],G.mean[1],'bo')
    pylab.plot([G.mean[0]-vec[0][0],G.mean[0]+vec[0][0]],
               [G.mean[1]-vec[0][1],G.mean[1]+vec[0][1]],'r-')

    pylab.plot([G.mean[0]-vec[1][0],G.mean[0]+vec[1][0]],
               [G.mean[1]-vec[1][1],G.mean[1]+vec[1][1]],'r-')



xx=array([sa.pars[0] for sa in ga.sample_list])
yy=array([sa.pars[1] for sa in ga.sample_list])
ww=array([sa.we for sa in ga.sample_list])

## now we plot
Np=100
cmin=-5.
cmax=5.
cstep=(cmax-cmin)/(1.0*Np)

sums=zeros((Np,Np))
wsums=zeros((Np,Np))
trvals=zeros((Np,Np))

for x,y,w in zip(xx,yy,ww):
    if (x<cmin) or (x>cmax) or (y<cmin) or (y>cmax):
        continue
    ix=int((x-cmin)/cstep)
    iy=int((y-cmin)/cstep)
    sums[iy,ix]+=1.0
    wsums[iy,ix]+=w

for i in range(Np):
    x=cmin+(i+0.5)*cstep
    for j in range(Np):
        y=cmin+(j+0.5)*cstep
        trvals[j,i]=exp(like([x,y]))

trvalsa=trvals/trvals.sum()
wsumsa=wsums/wsums.sum()
diffp=wsumsa-trvalsa
vmax=trvalsa.max()*1.1

pylab.subplot(2,2,1)
pylab.imshow(sums, interpolation='nearest', origin='lower left',extent=[cmin,cmax,cmin,cmax])
for G in ga.Gausses:
    plotel(G)
pylab.xlim(cmin,cmax)
pylab.ylim(cmin,cmax)
#pylab.colorbar()

pylab.subplot(2,2,2)
pylab.imshow(wsumsa, interpolation='nearest', origin='lower left',extent=[cmin,cmax,cmin,cmax],vmin=0, vmax=vmax)
pylab.colorbar()

pylab.subplot(2,2,3)
pylab.imshow(trvalsa, interpolation='nearest', origin='lower left',extent=[cmin,cmax,cmin,cmax],vmin=0, vmax=vmax)



pylab.subplot(2,2,4)
pylab.imshow(diffp, interpolation='nearest', origin='lower left',extent=[cmin,cmax,cmin,cmax])
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


pylab.savefig(sname)
pylab.show()



