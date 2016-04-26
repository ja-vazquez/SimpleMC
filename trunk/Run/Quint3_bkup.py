#!/usr/bin/env python
from scipy import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import pylab

H_0 =0.7
lam= 1.0
Ol0= 0.7
V0 = 0.5

def pot(x, i):
    if i==0: 
	return V0*exp(-lam*x)/H_0**2	
    if i==1: 
	return -lam*V0*exp(-lam*x)/H_0**2 
    else:
	print 'wrong choice' 
	stop()

def hub(x_vec):
    x, y, z  = x_vec
    return H_0*sqrt(0.5*x**2 + pot(x,0) +z)	

def harmonic(x_vec, t):
    x, y, z  = x_vec
    return [y/hub(x_vec), -3*y-pot(x,1)/hub(x_vec), -3*z]


lowr, highr = 0, 10**4
tol, tol1 =100, 100
c=0

t_output = np.linspace(-3, 0, 500)


def solver(x0):
    y0 = [0, x0, 0.3*(exp(-3)**(-3))]
    y_result = odeint(harmonic, y0, t_output)
    return y_result

while (abs(tol)>1e-3 and abs(tol1)>1e-3):
 mid=(lowr+highr)/2.0
 x0 = mid
 x1 = lowr 

 #y0 = [0, x0, 0.3*(exp(-3)**(-3))]
 y_result = solver(x0) #odeint(harmonic, y0, t_output) 

 #y1 = [0, x1, 0.3*(exp(-3)**(-3))]
 y1_result = solver(x1) #odeint(harmonic, y1, t_output)

 Omegal= (0.5*y_result[-1,1]**2+pot(y_result[-1,0],0))/(0.5*y_result[-1,1]**2+pot(y_result[-1,0],0)+y_result[-1,2])
 tol = Omegal - Ol0
 Omegal1= (0.5*y1_result[-1,1]**2+pot(y1_result[-1,0],0))/(0.5*y1_result[-1,1]**2+pot(y1_result[-1,0],0)+y1_result[-1,2])
 tol1 = Omegal1 - Ol0

 print 'Omegal',Omegal, Omegal1

 if(tol*tol1>0):
   lowr = mid
 else:
   highr = mid
 c+=1
 if (c>100):
    stop()

#-----------------------------

xx, yy, zz= y_result.T

#phi =interp1d(t_output, xx)
#phidot = interp1d(t_output, yy)
#Om = interp1d(t_output, zz)

fig = plt.figure(figsize=(9,8))
fig.add_subplot(1,1,1)

pylab.plot(t_output, (0.5*yy**2+pot(xx,0))/(0.5*yy**2+pot(xx,0)+zz),'k-')
pylab.plot(t_output, zz/(0.5*yy**2+pot(xx,0)+zz),'k--')
pylab.xlim(-3,0)
pylab.ylim(0,1)
plt.show()
