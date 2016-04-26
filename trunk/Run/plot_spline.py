#!/usr/bin/env python
from RunBase import *
import matplotlib.pyplot as plt
import matplotlib.ticker
import pylab 
T=LCDMCosmology()
T2=SplineLCDMCosmology()


zLOWZ  = 0.32
zCMASS = 0.57
zLyaA  = 2.34
zLyaC  = 2.36


zl=arange(0,2.6,0.05)
pylab.figure(figsize=(8,8))

pylab.subplot(3,1,1)
y1=[T.DaOverrd(z) for z in zl]
y2=[T2.DaOverrd(z) for z in zl]

pylab.errorbar(zCMASS, 9.519*(1+zCMASS), yerr=0.134*(1+zCMASS), color ='red', fmt='-o')
pylab.errorbar(zLyaA,  11.28*(1+zLyaA),  yerr=0.65*(1+ zLyaA),  color ='blue', fmt='-o')
pylab.errorbar(zLyaC,  10.8*(1+zLyaC),   yerr=0.4*(1+zLyaC),    color ='magenta', fmt='-o')

pylab.plot(zl,y1,'r-')
pylab.plot(zl,y2,'b-')
pylab.ylabel("$D_a(z)/r_d$")


pylab.subplot(3,1,2)
y1=[T.HIOverrd(z) for z in zl]
y2=[T2.HIOverrd(z) for z in zl]

pylab.errorbar(zCMASS,20.75, yerr=0.73,  color ='red', fmt='-o')
pylab.errorbar(zLyaA, 9.18,   yerr=0.28, color ='blue', fmt='-o')
pylab.errorbar(zLyaC, 9.0,    yerr=0.3,  color ='magenta', fmt='-o')

pylab.plot(zl,y1,'r-')
pylab.plot(zl,y2,'b-')
pylab.ylabel("$H^{-1}(z)/r_d$")


ax = plt.subplot(3,1,3)
ax.set_xscale('log')
y1=[1.0 for z in zl]
y2=[T2.Rho_de(1.0/(1.0+z)) for z in zl]
plt.plot(zl,y1,'r-')
plt.plot(zl,y2,'b-', label = "$\\rho_{de}(z)/\\rho_{de}$")
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.xaxis.set_minor_formatter(plt.ScalarFormatter())
ax.xaxis.set_minor_locator(plt.FixedLocator([0.2,0.5,2]))
pylab.xlim(0.09,2.5)


pylab.xlabel("$z$")
pylab.savefig("Spline_de.pdf")
pylab.show()
