#!/usr/bin/env python
from RunBase import *
import pylab
import numpy as np 
import math as N
import matplotlib.pyplot  as pyplot

T =LCDMCosmology(mnu=0)

T2=DecayLCDMCosmology()

Lambda_=  Lambda_par
Omx_= Omx_par
Omr_= Omr_par

Lambda_.setValue(0.01)
Omx_.setValue(0.2558)
Omr_.setValue(0.0001)

T2.updateParams([Lambda_, Omx_, Omr_])


zl=arange(0,10,0.1)

fig = pyplot.figure()


ax = fig.add_subplot(3,1,1)
y1=[T.RHSquared_a(1.0/(1.0+z)) for z in zl]
y2=[T2.RHSquared_a(1.0/(1.0+z)) for z in zl]
ax.plot(zl,y1,'b-')
ax.plot(zl,y2,'r-')
ax.set_ylabel("$H^2(z)$")

ax2= fig.add_subplot(3,1,2)
y1=[T.HIOverrd(z) for z in zl]
y2=[T2.HIOverrd(z) for z in zl]
ax2.plot(zl,y1,'b-')
ax2.plot(zl,y2,'r-')
ax2.set_ylabel("$H^{-1}(z)/r_d$")

ax3= fig.add_subplot(3,1,3)
y1=[T.DaOverrd(z) for z in zl]
y2=[T2.DaOverrd(z) for z in zl]
ax3.plot(zl,y1,'b-')
ax3.plot(zl,y2,'r-')
ax3.set_ylabel("$Da(z)/r_d$")


#T.printA()
#zl=arange(0,3,0.1)
#pylab.figure(figsize=(8,10))
#al=arange(0.0001, 0.9, 0.001)

#pylab.subplot(1,1,1)
#y1=[T.DaOverrd(1/(z+1)) for z in zl]
#y1=[T.Spline(a) for a in al]

#xl = np.logspace(0.0,-1.0,1000, base=N.exp(1))
#xnew  =np.logspace(0.0,-1.0,1000, base=N.exp(1))
#fig = pyplot.figure()
#ax = fig.add_subplot(1,1,1)
#ax.set_yscale('log')
#ax.set_xscale('log')
##a,b = T2.Spline(xnew)
#y=T2.RHSquared_a2(xnew)
#y2=T2.RHSquared_a(xl)
#ax.plot((xnew),y,'b-')
#ax.plot((xl),y2,'r-')
##ax.plot(xnew,a,label='$O_r$')
##ax.plot(xnew,b,label='$O_r$')


#al=np.logspace(-5.0,0.0,1000) #arange(0.0001, 0.9, 0.001)
#y1=[T.RHSquared_a(a) for a in al]
#y2=[T.RHSquared_a2(a) for a in al]
#fig = pyplot.figure()
#ax = fig.add_subplot(1,1,1)
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.plot(al,y1, al, y2)

pylab.xlabel("z")
show()
