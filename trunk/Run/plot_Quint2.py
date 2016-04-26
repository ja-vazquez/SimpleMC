#!/usr/bin/env python
#based on http://arxiv.org/pdf/0906.0396.pdf
# or http://arxiv.org/pdf/gr-qc/9711068v2.pdf
#Scalar Field potential V=V0*exp[-lam*\phi]

from RunBase import *
from scipy import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import pylab

#T=QuintCosmology()

if False:

 #T=QuintCosmology()
 T2=LCDMCosmology()

 zCMASS = 0.57
 zLyaA  = 2.34
 zLyaC  = 2.36

 pylab.figure(figsize=(8,8))
 pylab.errorbar(zCMASS,20.75*sqrt(zCMASS), yerr=0.73,color ='red', fmt='-o')
 pylab.errorbar(zLyaA, 9.18*sqrt(zLyaA),  yerr=0.28*zLyaA,color ='blue', fmt='-o')
 pylab.errorbar(zLyaC, 9.0*sqrt(zLyaC),   yerr=0.3*zLyaC,color ='magenta', fmt='-o')


 zl= arange(0.01, 3.0 ,0.05)

 lam_=lam_par
 B_=B_par
 
 lam_.setValue(6)
 B_.setValue(1) 
 T.updateParams([lam_])
 y1=[T.HIOverrd(z)*sqrt(z) for z in zl]

 lam_.setValue(6)
 B_.setValue(5) 
 T.updateParams([lam_])
 y2=[T.HIOverrd(z)*sqrt(z) for z in zl]

 lam_.setValue(6)
 B_.setValue(10)
 T.updateParams([lam_])
 y3=[T.HIOverrd(z)*sqrt(z) for z in zl]

 y4=[T2.HIOverrd(z)*sqrt(z) for z in zl]
 
 pylab.plot(zl,y4,'g-', label="LCDM")
 pylab.plot(zl,y1,'r-', label="$\lambda = 0.1$")
 pylab.plot(zl,y2,'k-',label="$\lambda  = 0.5$")
 pylab.plot(zl,y3,'b-',label="$\lambda  = 1.0$")
 pylab.legend(loc="lower left")

 plt.ylabel("D_H")
 plt.xlabel("$z$")

 pylab.show()


####
pts=[]
file='/astro/u/jvazquez/work/SimpleMC/trunk/chains/Quint/Quint_all2_phy_BBAO_7'
pnames=open(file+'.txt').readlines()
T=QuintCosmology()

number =3.0
number2=10.0
number3=0.00625
number4=1.001

#lpnames = len(pnames)
#for l in range(850,890): #lpnames):
#  vals =pnames[l].split()[1:]
#  number = 8.0 #float(vals[4])
#  number2 =float(vals[4])
#  number3 =0.00625 #float(vals[6])
#  number4 =float(vals[6])
#  print number	


if True: 
  lam_=lam_par
  V0_ =V0_par
  A_  =A_par
  B_  =B_par
 
  lam_.setValue(number)#8.0)
  V0_.setValue(number2) #10.0) 
  A_.setValue(number3) #0.00625)
  B_.setValue(number4) #1.001)

  #T=QuintCosmology()
  T.updateParams([lam_, V0_, A_, B_])
  
  sol = T.Ini_phi()
  xx, yy = sol.T

  H_0=0.6821
  Om0=0.3038
  Or0= Om0*1E-4
  Cte=sqrt(3.0)*H_0 
  lna = np.linspace(-25., 0, 200)

  pts.append((0.5*yy**2+T.Pot(xx,0)/Cte**2)*(1.0/T.hub((lna), sol.T))**2)

if True:
#   H_0=0.6821
#   Om0=0.3038
#   Or0= Om0*1E-4
#   Cte=sqrt(3.0)*H_0 
   #print 'pts', len(pts)
   fig = plt.figure(figsize=(6,4))
   ax=fig.add_subplot(1,1,1)
#   lna = np.linspace(-25., 0, 200)
   #ax.plot(lna, (0.5*yy**2+T.Pot(xx,0)/Cte**2)*(1.0/T.hub((lna), sol.T))**2, 'k-',label = "$\Omega_{\phi}$")
   for i in range(0,1):
     ax.plot(lna, pts[i], 'k-',label = "$\Omega_{\phi}$")
   #ax.plot(lna, Om0/(exp(lna)**3)*(1.0/T.hub((lna), sol.T))**2, 'b--', label = "$\Omega_{dm}$")
#   ax.plot(lna, Or0/(exp(lna)**4)*(1.0/T.hub((lna), sol.T))**2, 'r--', label = "$\Omega_r$")
   #ax.legend(loc="upper left")

   ylabel("$\\Omega$")
   xlabel("$\ln a$")
#   #ax.set_yscale('log')
   ax.axis([-15, 0, 0, 1.3])
  
#   ax2=fig.add_subplot(1,2,2)
#   ax2.plot(lna, xx, 'k-')
##   ax2.set_yscale('log')
##   ax2.plot(lna, (0.5*yy**2*Cte**2-T.Pot(xx,0))/(0.5*yy**2*Cte**2+T.Pot(xx,0)), 'k-')
#   ylabel("$\phi$")
#   xlabel("$\ln a$")
   plt.show()

