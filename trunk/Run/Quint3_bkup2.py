#!/usr/bin/env python
from scipy import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import pylab

H_0 =0.7
Ol0= 0.7
Om0= 0.299
Or0= 0.001

lna_ini = -9.0

V0 = .5
lam =1.

def Pot(x, i):
    if i==0: 
	return V0*exp(-lam*x)/H_0**2	
    if i==1: 
	return -lam*V0*exp(-lam*x)/H_0**2 
    else:
	print 'wrong choice' 
	stop()


def hub(a, x_vec):
    x, y, z  = x_vec
    return H_0*sqrt(0.5*y**2 + Pot(x,0) + Om0/a**3 + Or0/a**4)	


def RHS(x_vec, lna):
    a=exp(lna)	
    x, y, z = x_vec
    return [y/hub(a,x_vec), -3*y -Pot(x,1)/hub(a,x_vec), -3*z]


def solver(x0):
    y0 = [0, x0, 0.299*(exp(lna_ini)**(-3.0))]
    y_result = odeint(RHS, y0, lna)
    return y_result


lna = np.linspace(lna_ini, 0, 100)

lowr, highr = 0.0, 1e9
tol, tol1 =100, 100
Ttol= 1e-2
count=0


while (abs(tol)>Ttol):
 mid=(lowr+highr)/2.0 

 sol = solver(mid) 
 sol1 = solver(lowr)

 Omegal= (0.5*sol[-1,1]**2+Pot(sol[-1,0],0))*H_0**2/hub(1.0, sol[-1])**2 
 Omegal1= (0.5*sol1[-1,1]**2+Pot(sol1[-1,0],0))*H_0**2/hub(1.0, sol1[-1])**2 
 tol = Ol0 - Omegal
 tol1 = Ol0 - Omegal1
 print 'Omega_l',Omegal,Omegal1, 'mid', lowr, highr

 if(abs(tol) < Ttol):
   print 'nothing', abs(tol), abs(tol1)
   break
 else:
 #   mid = highr - Omegal1*((highr-lowr)/(Omegal1-Omegal)) 
 #   lowr = highr
 #   highr = mid
   if(tol*tol1>0):
     lowr = mid
   else:
     highr = mid

 count+=1
 if (count >50):
    stop()

#-----------------------------

xx, yy, zz= sol.T

#phi =interp1d(lna, xx)
#phidot = interp1d(lna, yy)
#Om = interp1d(lna, zz)

fig = plt.figure(figsize=(9,8))
fig.add_subplot(1,1,1)

pylab.plot(lna, (0.5*yy**2+Pot(xx,0))*H_0**2/hub(exp(lna), sol.T)**2,'k-')
pylab.plot(lna, zz*H_0**2/hub(exp(lna), sol.T)**2,'b-')
pylab.plot(lna, 0.7/(exp(lna)**0)*H_0**2/hub(exp(lna), sol.T)**2,'k--')
pylab.plot(lna, Om0/(exp(lna)**3)*H_0**2/hub(exp(lna), sol.T)**2,'b--')
pylab.plot(lna, Or0/(exp(lna)**4)*H_0**2/hub(exp(lna), sol.T)**2,'r--')
#pylab.xlim(-3,0)
pylab.ylim(0,1)
plt.show()
