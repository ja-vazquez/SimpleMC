#!/usr/bin/env python
# Scalar Field, tracking solutions
from pylab import *
from scipy import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import pylab

H_0 =0.7
Om0= 0.2999
Or0= 0.0001

lna_ini = -30.0
lna_fin = 0.0
 
lam = 3.8
V0 =1E30


def Pot(x, i):
    fact = 3*H_0**2
    if i==0: 
	return V0*exp(-lam*x)/fact
    if i==1: 
	return -lam*V0*exp(-lam*x)/fact
    else:
	print 'wrong choice' 
	stop()


def hub(a, x_vec):
    x, y  = x_vec
    return H_0*sqrt(0.5*y**2 + Pot(x,0) + Om0/a**3  + Or0/a**4)	


def RHS(x_vec, lna):
    a=exp(lna)	
    x, y = x_vec
    return [3.0*y*H_0/hub(a,x_vec), -3*y -Pot(x,1)*H_0/hub(a,x_vec)]

def solver(x0):
    phi0=0
    phidot0= 0
    print phidot0
    y0= [phi0, phidot0]
    y_result = odeint(RHS, y0, lna)
    return y_result


lna = np.linspace(lna_ini, lna_fin, 2000)

#lowr, highr = -1.e1, 1e2
#tol, tol1 =100, 100
#Ttol= 1e-3
#count=0

		#search initial conditions 

if True:
#while (abs(tol)>Ttol):
 mid= 1. #(lowr+highr)/2.0 

 sol = solver(mid) 
# sol1 = solver(lowr)

# Omegal= (0.5*sol[-1,1]**2+Pot(sol[-1,0],0))*H_0**2/hub(1.0, sol[-1])**2 
# Omegal1= (0.5*sol1[-1,1]**2+Pot(sol1[-1,0],0))*H_0**2/hub(1.0, sol1[-1])**2 
# tol = Ol0 - Omegal
# tol1 = Ol0 - Omegal1
# print 'Omega_l',Omegal, '\phi_ini= ', mid

# if(abs(tol) < Ttol):
#   print 'reach tolerance', abs(tol)
#   break
# else:
#   if(tol*tol1>0):
#     lowr = mid
#   else:
#     highr = mid
 
# count+=1
# if (count >50):
#    break

#-----------------------------

xx, yy = sol.T

IHub = interp1d(lna, H_0**2/hub(exp(lna), sol.T)**2)
#phi =interp1d(lna, xx)

fig = plt.figure(figsize=(10,5))

#ax = fig.add_subplot(1,3,1)
#ax.plot(lna, (0.5*yy**2+Pot(xx,0)), 'k-', label = "$\\rho_{\phi}$")
#ax.plot(lna, Om0/exp(lna)**3, 'k--', label="$\\rho_{DM}$")
#ax.legend(loc="lower left")
#ylabel("$\ln \\rho$")
#xlabel("$\ln a$")
#ax.set_yscale('log')

ax=fig.add_subplot(1,2,1)
ax.plot(lna, (0.5*yy**2+Pot(xx,0))*(H_0/hub(exp(lna), sol.T))**2, 'k-')
ax.plot(lna, Om0/(exp(lna)**3)*(H_0/hub(exp(lna), sol.T))**2, 'b--', label = "$\\rho_{dm}$") 
ax.plot(lna, Or0/(exp(lna)**4)*(H_0/hub(exp(lna), sol.T))**2, 'r--', label = "$\\rho_r$") 
ylabel("$\\Omega_{\phi}$")
xlabel("$\ln a$")
pylab.xlim([-13,0])

ax2=fig.add_subplot(1,2,2)
ax2.plot(lna, (H_0**2*0.5*yy**2-Pot(xx,0)/(3*H_0**2))/(H_0**2*0.5*yy**2+Pot(xx,0)/(3*H_0**2)), 'k-')
ylabel("$\phi/H_0$")
xlabel("$\ln a$")
pylab.xlim([-13,0])

#ax.plot(lna, Om0/(exp(lna)**3), 'b--', label = "$\\rho_{dm}$") 
#ax.plot(lna, Or0/(exp(lna)**4), 'r--', label = "$\\rho_r$") 
#plt.plot(lna, (IHub(lna)), 'r-', label="Hubble")

pylab.savefig('Quint_2.pdf')
plt.show()




