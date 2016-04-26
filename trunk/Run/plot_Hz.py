#!/usr/bin/env python
from RunBase import *
import math as N
import pylab 


params1 = {'backend': 'pdf',
               'axes.labelsize':28,
               'text.fontsize': 15,
               'xtick.labelsize': 15,
               'ytick.labelsize': 15,
               'legend.draw_frame': False,
               'legend.fontsize': 13,
               'lines.markersize': 6,
               'font.size': 18,
               'text.usetex': True}#
pylab.rcParams.update(params1)


T=LCDMCosmology()
T2=EarlyDECosmology() 
#T3=owa0CDMCosmology()
#T4=owa0CDMCosmology()
#T5=owa0CDMCosmology()

#w_=w_par
#h_=h_par

#w_.setValue(-1)
#h_.setValue(0.72)
#T.updateParams([w_,h_])

#w_.setValue(0.0)
#h_.setValue(0.72)
#T2.updateParams([w_,h_])

#w_.setValue(0.02)
#h_.setValue(0.72)
#T3.updateParams([w_,h_])

#w_.setValue(0.04)
#h_.setValue(0.72)
#T4.updateParams([w_,h_])

#w_.setValue(0.06)
#h_.setValue(0.72)
#T5.updateParams([w_,h_])

fig = plt.figure(figsize=(9,8))

zl= arange(0.01,1000.,0.5)
fig.add_subplot(1,1,1)

for z in zl:
  print z, T.RHSquared_a(1/(z+1))**0.5, T2.RHSquared_a(1/(z+1))**0.5

y1=[T.RHSquared_a(1/(1+z))**0.5 for z in zl]
#y2=[T2.RHSquared_a(1/(1+z))**0.5 for z in zl]
#y3=[T3.RHSquared_a(1/(1+z))**0.5 for z in zl]
#y4=[T4.RHSquared_a(1/(1+z))**0.5 for z in zl]
#y5=[T5.RHSquared_a(1/(1+z))**0.5 for z in zl]

pylab.plot(zl,y1,'k-',label="LCDM")
#pylab.plot(zl,y2,'r-',label="n=0")
#pylab.plot(zl,y3,'b-',label="n=0.02")
#pylab.plot(zl,y4,'g-',label="n=0.04")
#pylab.plot(zl,y5,'y-',label="n=0.06")

#pylab.errorbar(0,(70.6/69), yerr=3.3*(1.0/69))
#pylab.errorbar(0.57,(92.4/69), yerr=4.5*(1.0/69))
#pylab.errorbar(2.34,(222.0/69), yerr=7*(1.0/69))

#zCMASS = 0.57
#zLyaA  = 2.34
#zLyaC  = 2.36

#pylab.errorbar(zCMASS,20.75, yerr=0.73,color ='red', fmt='-o')
#pylab.errorbar(zLyaA, 9.18,  yerr=0.28,color ='blue', fmt='-o')
#pylab.errorbar(zLyaC, 9.0,   yerr=0.3,color ='magenta', fmt='-o')

pylab.ylabel("$H(z)/H_0$")
pylab.legend(loc="upper left")
pylab.xlabel("z")
#pylab.loglog()

#fig.add_subplot(2,1,2)

#w1=[-1 for z in zl]
#w2=[T2.w_a(1/(1+z)) for z in zl]
#w3=[T3.w_a(1/(1+z)) for z in zl]
#w4=[T4.w_a(1/(1+z)) for z in zl]
#w5=[T5.w_a(1/(1+z)) for z in zl]

#pylab.plot(zl,w1,'k-',label="lcdm")
#pylab.plot(zl,w2,'r-',label="a=-6")
#pylab.plot(zl,w3,'b-',label="a=-3")
#pylab.plot(zl,w4,'g-',label="a=3")
#pylab.plot(zl,w5,'y-',label="a=6")

#pylab.ylabel("$w_{\\rm eff}(z)$")
#pylab.legend(loc="upper left")
#pylab.xlabel("z")
#pylab.ylim(-2,1.5)
#pylab.semilogx()


pylab.savefig("Hz.pdf")
pylab.show()
