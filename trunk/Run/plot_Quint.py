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
from matplotlib import gridspec
import pylab, sys
import scipy.optimize as optimize
import Useful as Usf

params1 = Usf.setting_plot()
pylab.rcParams.update(params1)


name = 'Quint_vl_phi'


pts, pts2, pts3, pts4 =[], [], [], []
Oma, Ora = [], []
lam = []

H_0, Om0 = 0.6821, 0.3038
Cte=sqrt(3.0)*H_0
Or0= sqrt(3.05)*Om0*1E-4


T=QuintCosmology_try()

lna = np.linspace(-20., 5, 300)
z = np.exp(-lna)

lam_, V0_ = lam_par, V0_par
A_, lB_   = A_par, lB_par
V0_f = 1.0


def ini_phi(xx):
     return (T.Pot(xx,1)**2/T.Pot(xx,0)/4. - T.Pot(xx,0))*exp(-140)*4.E4*(3/2.36) - 1.


A_f, B_f   = 0.5, 10.

for rr in range(4):
 lam_f = 3 + rr/2.
 lam.append(lam_f) 

 if True:
   lam_.setValue(lam_f)
   A_.setValue(A_f)
   lB_.setValue(B_f) 
   V0_.setValue(V0_f)

   T.updateParams([lam_, V0_, A_, lB_])
  
   sol = T.Ini_phi()
   xx, yy = sol.T
 
   pts.append(T.O_phi(lna))
   pts3.append(T.w_ede(lna))
   Oma.append(Om0/(exp(lna)**3)*(1.0/T.hub((lna), sol.T))**2)
   Ora.append(Or0/(exp(lna)**4)*(1.0/T.hub((lna), sol.T))**2)
   pts4.append(xx) 

if True:
   fig = plt.figure(figsize=(16,5))
   ax=fig.add_subplot(1,3,1)
   for i in range(4):
       ax.plot(z, pts[i], color = Usf.colour(i+1))
       ax.plot(z, Oma[i], 'g--')
       ax.plot(z, Ora[i], 'y--')
   plt.ylim([-0.01,1.01])
   ax.set_xscale('log')
   plt.xlim([1.e-1, 1e8])
   plt.xlabel("$z$")
   plt.ylabel("$\\Omega_{ede}$")
   ax.grid(True)
   
   ax2=fig.add_subplot(1,3,2)
   for i in range(4):
       ax2.plot(z, pts3[i], color = Usf.colour(i+1))
   plt.ylim(ymin=-1.01)
   ax2.set_xscale('log')
   plt.xlim([1.e-1, 1e8])
   plt.xlabel("$z$")
   plt.ylabel("$w$")
   ax2.grid(True)

   ax3=fig.add_subplot(1,3,3)
   for i in range(4):
       ax3.plot(z, pts4[i], color = Usf.colour(i+1), 
	label = "$\\lambda$ =%1.1f"%(lam[i]))
   plt.legend(loc="lower left")
   plt.xlim([1.e-1, 1e8])
   ax3.set_xscale('log')
   plt.xlabel("$z$")
   plt.ylabel("$\phi$")
   ax3.grid(True)

plt.savefig(name+".pdf")
plt.show()




if False: 
   ax4=fig.add_subplot(1,3,4)
   for i in range(4):
       ax4.plot(lna, pts4[i])
   plt.xlabel("$\ln a$")
   plt.ylabel("$\dot \phi^2$")
   pylab.yscale('log')
   ax4.grid(True)

   ax5=fig.add_subplot(2,3,5)
   ax5.plot(lna, T.Pot(xx,0))
   plt.xlabel("$\ln a$")
   plt.ylabel("$V(\phi)$")
   pylab.yscale('log') 
   ax5.grid(True)

   #ax6=fig.add_subplot(2,3,6)
   #ax6.plot(xx, (T.Pot(xx,1)**2/T.Pot(xx,0)/4. - T.Pot(xx,0))*exp(-140)*4.E4*(3/2.36))
   #plt.xlim([-35,-30])   
   #plt.xlabel("$\ln a$")
   #plt.ylabel("$V,{\phi}(\phi)$")
   #pylab.yscale('log')
   #ax6.grid(True)
