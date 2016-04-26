#!/usr/bin/env python
from RunBase import *
import math as N
import pylab 

T2=EarlyDECosmology()
T3=EarlyDECosmology()
T4=EarlyDECosmology()
T5=EarlyDECosmology()

T=LCDMCosmology() 
#Obh2=0.022032,Om=0.3183,h=0.6704)

print T.RHSquared_a(1.0)


Om_=Om_par
Om_.setValue(0.3)

Ode_=Ode_par
Ode_.setValue(0.1)
T2.updateParams([Ode_,Om_])

Ode_.setValue(0.2)
T3.updateParams([Ode_,Om_])

Ode_.setValue(0.3)
T4.updateParams([Ode_,Om_])

h_=h_par
Ode_.setValue(0.1)
Om_.setValue(0.45)
h_.setValue(0.55)
T5.updateParams([Ode_,Om_,h_])

zLOWZ  = 0.32
zCMASS = 0.57
zLyaA  = 2.34-0.04
zLyaC  = 2.36+0.04
z_CMB = 1090.43

zl= arange(0,4,0.2)
fig = plt.figure(figsize=(9,6))

#def plot_errorbar(z,val, yerr=0, color=0, fmt=0, markersize=0,label=None, empty=True):
#    if empty:
#        mfc='white'
#    else:
#        mfc=color
#    pylab.errorbar(z,val/fixer(z), yerr=yerr/fixer(z), color=color, fmt=fmt, markersize=markersize, lw=2,markerfacecolor=mfc)
#    if label>0:
#        if (mfc=='white'):
#            pylab.plot ([],[],fmt,color='black',label=label,markersize=markersize,markerfacecolor=mfc)
#        else:
#            pylab.plot ([],[],fmt,color='black',label=label,markersize=markersize)

ax = fig.add_subplot(1,1,1)

y1=[T.DaOverrd(z) for z in zl]
y2=[T2.DaOverrd(z)/(((1-0.0*T2.Omega_de(1/(1+z)))**(0.5))) for z in zl]
y3=[T3.DaOverrd(z)/(((1-0.0*T3.Omega_de(1/(1+z)))**(0.5))) for z in zl]
y4=[T4.DaOverrd(z)/(((1-0.0*T4.Omega_de(1/(1+z)))**(0.5))) for z in zl]

#y4=[T4.DaOverrd(z)/(((1-0.0*T4.Omega_de(1/(1+z)))**(0.5))*fixer(z)) for z in zl]
#y3=[T3.DaOverrd(z)/(((1-0.0*T3.Omega_de(1/(1+z)))**(0.5))*fixer(z)) for z in zl]
#y2=[T2.DaOverrd(z)/(((1-0.0*T2.Omega_de(1/(1+z)))**(0.5))*fixer(z)) for z in zl]
#y1=[T.DaOverrd(z)/(fixer(z)) for z in zl]
#y5=[T5.DaOverrd(z)/(((1-0.0*T5.Omega_de(1/(1+z)))**(0.5))*fixer(z)) for z in zl]
pylab.plot(zl,y1,'k-')
#pylab.plot(zl,y4,'b-',label="Ode=0.1")
#pylab.plot(zl,y3,'g-',label="Ode=0.3")
#pylab.plot(zl,y2,'r-',label="Ode=0.5") 
#pylab.plot(zl,y5,'m-',label="Ode=0.2")
#ax.set_xscale('log')
#pylab.ylim(0,50.0)
#pylab.xlim(0.09,4.0)
pylab.ylabel("$D_a(z)/r_d$") #\\sqrt{z}$")

print 'z=2.34 ',     '  Jose AV    ', 'Hee-Jong'
#print 'Ode=0.1',  T4.DaOverrd(2.34), '39.98210'
print 'Ode=0.2',  T5.DaOverrd(2.34), '40.92519'
#print 'Ode=0.3',  T3.DaOverrd(2.34), '42.19994'

X=[2.34,2.34,2.34,2.34]
X2=[0.57,0.57,0.57,0.57]
Y=[39.0357,39.89219,40.92519, 42.19994]
Y2=[14.66733,15.26124,15.97122,16.83849]
pylab.plot(X,Y,'go')
pylab.plot(X2,Y2,'bo')


plot_errorbar(zCMASS, 9.519*(1+zCMASS), yerr=0.134*(1+zCMASS), color ='red', fmt='d', markersize=6, empty=False)
plot_errorbar(zLyaA,  11.28*(1+zLyaA),  yerr=0.65*(1+ zLyaA),  color ='red', fmt='o', markersize=8, empty=False)
plot_errorbar(zLyaC,  10.8*(1+zLyaC),   yerr=0.4*(1+zLyaC),    color ='red', fmt='*', markersize=8, empty=False)
#plot_errorbar(z_CMB,  94.22,            yerr=0.6,              color ='red', fmt='-p', markersize=6, empty=False)

pylab.legend(loc="upper left")

#ax2 = fig.add_subplot(3,1,2)
#y5=[T5.HIOverrd(z)*z/(((1-T5.Omega_de(1/(1+z)))**(0.5))*fixer(z)) for z in zl]
#y4=[T4.HIOverrd(z)*z/(((1-T4.Omega_de(1/(1+z)))**(0.5))*fixer(z)) for z in zl]
#y3=[T3.HIOverrd(z)*z/(((1-T3.Omega_de(1/(1+z)))**(0.5))*fixer(z)) for z in zl]
#y2=[T2.HIOverrd(z)*z/(((1-T2.Omega_de(1/(1+z)))**(0.5))*fixer(z)) for z in zl]
#y1=[T.HIOverrd(z)*z/fixer(z) for z in zl]
#pylab.plot(zl,y1,'k-',zl,y2,'r-',zl,y3,'g-',zl,y4,'b-',zl,y5,'m-')
#pylab.ylabel("$D_H(z)/(r_d*(1-\\Omega_{de}))$") #\\sqrt{z}$")
#ax2.set_xscale('log')
##pylab.ylim(0,30.0)
#pylab.xlim(0.09,4.0)


#plot_errorbar(zCMASS,20.75*zCMASS, yerr=0.73*zCMASS,color ='blue', fmt='-d', markersize=6,empty=False)
#plot_errorbar(zLyaA, 9.18*zLyaA,   yerr=0.28*zLyaA, color ='blue', fmt='-o', markersize=8,empty=False)
#plot_errorbar(zLyaC, 9.0*zLyaC,    yerr=0.3*zLyaC,  color ='blue', fmt='-*', markersize=8, empty=False)


#ax3 = fig.add_subplot(2,1,2)

#y5=[T5.Omega_de(1/(1+z)) for z in zl]
#y4=[T4.Omega_de(1/(1+z)) for z in zl]
#y3=[T3.Omega_de(1/(1+z)) for z in zl]
#y2=[T2.Omega_de(1/(1+z)) for z in zl]
#pylab.plot(zl,y2,'r-',zl,y3,'g-',zl,y4,'b-',zl,y5,'m-')
#pylab.ylabel("$\\Omega_d(z)$")
##ax3.set_xscale('log')
##pylab.xlim(0.09,4.0)
pylab.xlabel("z")

pylab.savefig("Early.pdf")
pylab.show()
