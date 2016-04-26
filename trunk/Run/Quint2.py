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


if True:
 lam=1E-2

 Ol= 0.72
 Om = 0.278
 h0 =0.72

 ini = [0, sqrt(Ol), sqrt(Om), h0]

 def harmonic(x_vec, t):
    x, y, z, h  = x_vec
    func = 1.5*(2*x**2 + z**2+ 1.33*(1 -x**2 -y**2 -z**2)) 
    return [-3*x+ lam*sqrt(1.5)*y**2+ x*func, -lam*sqrt(1.5)*x*y + y*func, func*z-1.5*z, -h*func] 


 t_output = np.linspace(0, -6.0, 500)
 y_result = odeint(harmonic, ini, t_output)

 xx, yy, zz, hh= y_result.T 

 plt.subplot(311)
 oo = 1- xx**2 - yy**2 - zz**2 
 plt.plot(t_output, xx**2+yy**2, 'k-',label="$\Omega_{\phi}$")
 plt.plot(t_output, zz**2,'b-', label = "$\Omega_{dm}$")
 plt.plot(t_output, oo,'r-', label = "$\Omega_r$")
 plt.ylabel("$\Omega_i$")
 pylab.legend(loc="upper left")

 plt.subplot(312)
 pylab.plot(t_output, (xx**2-yy**2)/(xx**2 + yy**2),'k-',label="$w_{\phi}$")
 plt.ylabel("$w_{\phi}$")
 pylab.legend(loc="upper left")

 #interpol
# logar= linspace(0.0,-1.5,500)
# ilogar= logar[::-1]
# hubble = interp1d(ilogar, hh[::-1])


# plt.subplot(313)
# pylab.plot(ilogar, (hubble(ilogar)), 'r-', label="Hubble")
# pylab.legend(loc="upper right")
# plt.ylabel("Hubble")

 plt.xlabel("$\ln a$")
 plt.show()
