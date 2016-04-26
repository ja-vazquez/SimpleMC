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

name = 'Quintb_vl2'

params1 = {'backend': 'pdf',
               'axes.labelsize': 20,
               'text.fontsize': 18,
               'xtick.labelsize': 20,
               'ytick.labelsize': 20,
               'legend.draw_frame': False,
               'legend.fontsize': 12,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}#
pylab.rcParams.update(params1)

def colour(x):
    if x==1: return 'red'
    if x==2: return 'blue'
    if x==3: return 'black'
    if x==4: return 'magenta'
    if x==5: return 'cyan'
    if x==6: return 'orange'
    if x==7: return 'green'
    if x==8: return 'yellow'
    if x==9: return 'purple'
    if x>9: return 'black' # print("Increased colouring")

rr =5

# if vary = 1, otherwise =0
vl, vA, vB = 0, 1, 0

if vl == 1: lam_f = 3.0 
else: lam_f = 6.0

if vA == 1: A_f = 0.5
else: A_f   = 0.5 #0.005
        
if vB == 1: B_f = 5. #0.5
else: B_f   = 7. 

pts, pts2, pts2b, pts3, pts4, pts5 =[], [], [], [], [], []

T=QuintCosmology_try()
TL =  LCDMCosmology()

#Ini values
H_0, Om0 = 0.6821, 0.3038
Cte=sqrt(3.0)*H_0
Or0= sqrt(3.05)*Om0*1E-4
lna = np.linspace(-35., 10, 300)

zz = np.linspace(0, 1./exp(-25)-1, 200)
zll = np.linspace(-4, 0, 200)
zl = np.logspace(-1,1,50)

lam_, V0_ = lam_par, V0_par
A_, lB_   = A_par, lB_par


vall, valA, valB = [], [], []

#lam, bb = ['4','5','6','7','10'], ['17','13','11','10','6']
for r in range(0,rr):
   lam_f, A_f, B_f =lam_f + vl*r, A_f + vA*0.2*r, B_f + vB*1.*r
#   lam_f, A_f, B_f = float(lam[r]), 0.5, float(bb[r])
   vall.append(lam_f)
   valA.append(A_f) 
   valB.append(B_f)

   lam_.setValue(lam_f)
   A_.setValue(A_f)
   lB_.setValue(B_f) 

   T.updateParams([lam_, V0_, A_, lB_])
#   y4=[T.HIOverrd(z) for z in zl] 
   sol = T.Ini_phi()
   xx, yy = sol.T
  
   Omegal= (0.5*sol[-1,1]**2+T.Pot(sol[-1,0],0)/Cte**2)/T.hub(0.0, sol[-1])**2
  
   if True: #abs(Omegal -0.7) < 0.1:
      pts.append(T.O_phi(lna)) 
      pts2.append(Om0/(exp(lna)**3)*(1.0/T.hub((lna), sol.T))**2)
      pts2b.append(Or0/(exp(lna)**4)*(1.0/T.hub((lna), sol.T))**2)
      pts3.append((0.5*yy**2*Cte**2-T.Pot(xx,0))/(0.5*yy**2*Cte**2+T.Pot(xx,0)))     
      pts4.append(xx)
      pts5.append(T.Pot(xx,0))
#      y1=[T.DaOverrd(z) for z in zl]
#      pts3.append(y1)
#      y2=[T.HIOverrd(z) for z in zl]
#      pts4.append(y2)

if True:
   fig = plt.figure(figsize=(20,6))
   #gs = gridspec.GridSpec(1,4, width_ratios=[1, 1, 1, 91]) 
   #ax=fig.add_subplot(gs[0])
   ax = plt.subplot(1, 3, 1)

   for i in range(0,len(pts)):
     ax.plot(lna, pts[i], color = colour(i+1))
	#label = "$\\lambda$ =%1.1f,  A=%1.3f,  B=%1.1f"%(vall[i],valA[i],valB[i]))
     ax.plot(lna, pts2[i], 'g--')
     ax.plot(lna, pts2b[i], 'y--')	
   ax.grid(True)
   #plt.legend(loc="upper left")
   ylabel("$\\Omega$")
   xlabel("$\ln a$")
   #ax.axis([-15, 0, 0, 1.3])
   plt.ylim([-0.01,1.01])
   plt.xlim([-30,5])
    
   #ax2=fig.add_subplot(gs[1]) 
   ax2= plt.subplot(1, 3, 2)
   for i in range(0,len(pts)):   
     ax2.plot(lna, np.array(pts3[i]), color = colour(i+1))
#	label = "$\\lambda$ =%1.1f,  A=%1.3f,  B=%1.1f"%(vall[i],valA[i],valB[i]))
   #plt.legend(loc="upper left")
   ax2.grid(True)
   plt.ylim([-1.01,1.01]) 
   ylabel("$w$")
   xlabel("$\ln a$")
   plt.xlim([-30,5])

   #ax3=fig.add_subplot(gs[2])
   ax3 = plt.subplot(1, 3, 3)
   for i in range(0,len(pts)):
      ax3.plot(lna, np.array(pts4[i]), color = colour(i+1),
	label = "$\\lambda$ =%1.0f,  $A\\lambda^2$=%1.1f,  $B$=%1.0f"%(vall[i],valA[i],valB[i]))
   plt.legend(loc="upper left")
   ax3.grid(True)
   ylabel("$\phi$")  
   xlabel("$\ln a$")
#   plt.xlim([-30,5])

   #ax4 = plt.subplot(2, 2, 4)
   #for i in range(0,len(pts)):
   #   ax4.plot(lna, np.array(pts5[i]), color = colour(i+1)) 	
   #ax4.set_yscale('log') 

   plt.tight_layout()
   plt.savefig(name+".pdf")
   plt.show()



