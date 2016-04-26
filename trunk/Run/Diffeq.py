#!/usr/bin/env python

from scipy.integrate import odeint
from pylab import *
import matplotlib.pyplot  as pyplot


H0  = 0.7
Ox0 = 0.19
Or0 = 0.1
Om0 = 0.3
OL0 = 1-Ox0-Or0-Om0

Lambda = 0.01


def Hub2(y,t):
   return H0*H0*(OL0+y[:,0]+y[:,1]+y[:,2]) 

def deriv(y,t):
   factor = Lambda*y[0]/(t*H0*(OL0+y[0]+y[1]+y[2])**(0.5))
   return array([ -3*y[0]/t -factor, -4*y[1]/t +factor, -3*y[2]/t ])    



time= linspace(1.0,0.0001,100000)
yinit= array([Ox0, Or0, Om0])
y =odeint(deriv,yinit,time)

fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')

ax.plot(time,y[:,0]*H0*H0/Hub2(y,time),label='$O_x$')
ax.plot(time,y[:,1]*H0*H0/Hub2(y,time),label='$O_r$')
ax.plot(time,y[:,2]*H0*H0/Hub2(y,time),label='$O_{dm}$')
ax.plot(time,OL0*H0*H0/Hub2(y,time),label='$O_L$')
text(0.01, .9, r'$\lambda=-1$')
legend(loc='upper right')
xlabel('t')
ylabel('y')
show()




