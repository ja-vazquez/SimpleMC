#!/usr/bin/env python
from RunBase import *
import pylab 
T=LCDMCosmology(mnu=0)


zLOWZ  = 0.32 
zCMASS = 0.57
zLyaA  = 2.34
zLyaC  = 2.36

zl=arange(0,3,0.1)
pylab.figure(figsize=(8,10))



pylab.subplot(3,1,1)
y1=[T.DaOverrd(z) for z in zl]
pylab.errorbar(zCMASS,9.519*(1+zCMASS), yerr=0.134*(1+zCMASS),color ='red', fmt='-o')
pylab.errorbar(zLyaA, 11.28*(1+zLyaA),  yerr=0.65*(1+ zLyaA),  color ='blue', fmt='-o')
pylab.errorbar(zLyaC, 10.8*(1+zLyaC),   yerr=0.4*(1+zLyaC),   color ='magenta', fmt='-o')
pylab.plot(zl,y1,'k-')
pylab.ylabel("$D_a(z)/r_d$")

pylab.subplot(3,1,2)
y1=[T.HIOverrd(z) for z in zl]
pylab.errorbar(zCMASS,20.75, yerr=0.73,color ='red', fmt='-o')
pylab.errorbar(zLyaA, 9.18,  yerr=0.28,color ='blue', fmt='-o')
pylab.errorbar(zLyaC, 9.0,   yerr=0.3,color ='magenta', fmt='-o')
pylab.plot(zl,y1,'k-')
pylab.ylabel("$H^{-1}(z)/r_d$")

pylab.subplot(3,1,3)
y1=[T.DVOverrd(z) for z in zl]
pylab.errorbar(zLOWZ,8.467,yerr=0.167,color ='green', fmt='-o')
pylab.plot(zl,y1,'k-')
pylab.ylabel("$D_v(z)/r_d$")

pylab.plot([],[],'g-',label='LOWZ')
pylab.plot([],[],'r-',label='CMASS')
pylab.plot([],[],'b-',label='Lyman-$\\alpha$ auto')
pylab.plot([],[],'magenta',label='Lyman-$\\alpha$ cross')
pylab.legend(loc='lower right')

pylab.xlabel("z")
pylab.savefig("Fig1.pdf")
pylab.show()
