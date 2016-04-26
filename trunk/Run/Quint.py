#!/usr/bin/env python
#This only LCDM model
from scipy import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import pylab

H_0 =0.7
Om0= 0.3

#def hub(x_vec):
#    return H_0*sqrt(dot(x_vec,x_vec))

def harmonic(x_vec, t):
    x, y  = x_vec
    return [0, -3*y]

lowr = 0
highr = 1
tol =100
c=0

t_output = np.linspace(-3, 0, 500)

while (abs(tol)>1e-3):
 mid=(lowr+highr)/2.0
 Om = mid
 
 y0 = [1.0-Om, Om]
 sol = odeint(harmonic, y0, t_output) 

 Omegam= (sol[-1,1])/(sol[-1,0]+sol[-1,1])
 tol = Omegam - Om0
 print Omegam
 if(tol >0): 
   highr = mid
 else:
   lowr = mid  
 c+=1
 if (c>100):
    stop()

#-----------------------------

xx, yy = sol.T

Om = interp1d(t_output, yy)
Ol = interp1d(t_output, xx)

fig = plt.figure(figsize=(9,8))
fig.add_subplot(1,1,1)

print 'Omega_m today',(Om(0.0))/(Om(0.0)+Ol(0.0))

pylab.plot(t_output, (xx)/(xx+yy),'k-')
pylab.plot(t_output, yy/(xx+yy),'k--')
pylab.xlim(-3,0)
pylab.ylim(0,1)
plt.show()
