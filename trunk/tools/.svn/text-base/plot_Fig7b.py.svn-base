#!/usr/bin/env python

#Error bars mainly got them from 
#table 3 at
#http://arxiv.org/pdf/1108.2635v1.pdf

from RunBase import *
import matplotlib.pyplot as plt
import matplotlib.ticker
import pylab 
import math as N

pylab.figure(figsize=(6,6))
params1 = {'backend': 'pdf',
               'axes.labelsize': 30,
               'text.fontsize': 20,
               'xtick.labelsize': 22,
               'ytick.labelsize': 22,
               'legend.draw_frame': False,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}#
pylab.rcParams.update(params1)


T=LCDMCosmology(Obh2=0.0223983,Om=0.302422,h=0.681841)
T2=WeirdCDMCosmology()

mu_=  mu_par
amp_= Amp_par
sig_= sig_par

mu_.setValue(2.04296)
amp_.setValue(-0.0638545)
sig_.setValue(0.266033)

T2.updateParams([amp_, sig_,mu_])

zz=10**(arange(-1,2,0.01))
de=[]
print zz
Tx=T2
Ol=1-Tx.Om
for z in zz:
    a=1./(1+z)
    mattercont = LCDMCosmology.RHSquared_a(Tx,a)-Ol
    lambdacont = Tx.RHSquared_a(a)-mattercont
    de.append(lambdacont)

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(1,1,1)


pylab.plot (zz,de,'r-',lw=3, label='Tuned Oscillation')
pylab.plot (zz,-1.*array(de),'r--',lw=3)

pylab.plot (zz,[1-T.Om]*size(zz),'b-',lw=3, label='$\\Lambda$CDM')
pylab.xlim(0.1,20.)
pylab.ylim(0.1,10.)
pylab.loglog()
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))

pylab.xlabel('$z$')
pylab.ylabel('$\\rho_{\\rm DE}(z)/\\rho_{\\rm crit}(z=0)$')
pylab.tight_layout()
pylab.legend(loc='upper left')
pylab.tight_layout()

pylab.savefig('W2.pdf')

pylab.show()

