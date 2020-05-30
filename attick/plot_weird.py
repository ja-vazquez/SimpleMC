#!/usr/bin/env python
from RunBase import *
import pylab 
T=LCDMCosmology()


T2=WeirdCDMCosmology()
#T2=LCDMCosmology()

mu_=mu_par
amp_=Amp_par
sig_=sig_par

mu_.setValue(3.0)
amp_.setValue(0.05)
sig_.setValue(0.7)
T2.updateParams([amp_, sig_,mu_])

zl=exp(arange(0,8,0.1))
pylab.figure(figsize=(8,8))

pylab.subplot(3,1,3)
#y1=[T.RHSquared_a(1.0/(1.0+z)) for z in zl]
y2=[T2.RHSquared_a(1.0/(1.0+z)) for z in zl]
y1=[T2.Om*(1+z)**3 for z in zl]
print y2,y1
pylab.plot(zl,array(y2)-array(y1),'r-')
pylab.ylabel("\Delta H^2(z)")
pylab.loglog()
pylab.subplot(3,1,1)
y1=[T.DaOverrd(z) for z in zl]
y2=[T2.DaOverrd(z) for z in zl]
pylab.plot(zl,y1,'r-')
pylab.plot(zl,y2,'b-')
pylab.ylabel("$Da(z)/r_d$")
pylab.loglog()

pylab.subplot(3,1,2)
y1=[T.HIOverrd(z) for z in zl]
y2=[T2.HIOverrd(z) for z in zl]
pylab.plot(zl,y1,'r-')
pylab.plot(zl,y2,'b-')
pylab.ylabel("$H^{-1}(z)/r_d$")
pylab.loglog()
#pylab.subplot(3,1,3)
#y1=[T.RHSqua(z) for z in zl]
#y2=[T2.DVOverrd(z) for z in zl]
#pylab.plot(zl,y1,'r-')
#pylab.plot(zl,y2,'b-')
#pylab.ylabel("Dv(z)/rd")


pylab.xlabel("z")
pylab.savefig("compensated.pdf")
pylab.show()
