#!/usr/bin/env python
from RunBase import *
import pylab

D=DecayLCDMCosmology()


color=['red','blue','green','cyan','black']

for lw,lam in enumerate([0,0.05, 0.1, 0.2, 0.3]):
    #h_par.setValue(0.68-lam/2.0)
    lambda_par.setValue(lam)
    print "Updaring pars..."
    D.updateParams([lambda_par,h_par])
    print "done, plotting"
    pylab.subplot(2,1,1)
    #zl=arange(0,5.01,0.1)
    zl=arange(0,5.01,2.5)
    print zl
    y1=[D.rx(-log(1+z)) for z in zl]
    y2=[D.rr(-log(1+z)) for z in zl]
    #pylab.plot(D.logar, D.sol[:,0],'b-',lw=2,color=color[lw])
    #pylab.plot(D.logar, D.sol[:,1],'r-',lw=2,color=color[lw])
    pylab.plot(zl,y1,lw=2,color=color[lw])
    pylab.plot(zl,y2,'r-',lw=2,color=color[lw])

    pylab.subplot(2,1,2)

    print "DA..."
    y1=[D.DaOverrd(z)/sqrt(z)   for z in zl]
    print "HI..."
    y2=[D.HIOverrd(z)*z/sqrt(z) for z in zl]
    print "DV..."
    y3=[D.DVOverrd(z)/sqrt(z)   for z in zl]

    pylab.plot(zl,y1,'-',lw=2,color=color[lw])
    pylab.plot(zl,y2,':',lw=2,color=color[lw])
    pylab.plot(zl,y3,'--',lw=2,color=color[lw])

    print D.RHSquared_a(1)


pylab.subplot(2,1,1)
pylab.ylabel('r_x, r_r')

pylab.subplot(2,1,2)
pylab.ylabel('D_x/z^{1/2}')
pylab.xlabel('z')
plt.tight_layout()
pylab.savefig('ddmplot.pdf')
pylab.show()
