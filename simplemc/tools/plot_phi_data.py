#!/usr/bin/env python

#Error bars mainly got from  table 3,
#http://arxiv.org/pdf/1108.2635v1.pdf
#and MGS paper


from simplemc.models.PhiCDMCosmology import PhiCosmology
from simplemc.models.LCDMCosmology import LCDMCosmology
#from owa0CDMCosmology import owa0CDMCosmology
#from ParamDefs import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab
#import sys

model = 'pow'


params1 = {'backend': 'pdf',
               'axes.labelsize':  18,
               'xtick.labelsize': 18,
               'ytick.labelsize': 15,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 18,
               'text.usetex': True}
pylab.rcParams.update(params1)
zl=np.arange(0, 3.0, 0.05)
lna   = np.linspace(-6, 0, 500)

plaw=0.5
eboss=False

def fixer(z):
    if plaw>0:
        return z**plaw
    else:
        return log(1.+z)

def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
                txt.set_color(line.get_color())

#-------------------

#Planck best fit cosmology
T2=LCDMCosmology()
x1=[0.68/(T2.RHSquared_a(np.exp(la)))      for la in lna]
#z1= [0.7/np.sqrt(T2.RHSquared_a(np.exp(la)) for la in np.linspace(-6, 0, 500)]

hh = []
ww = []
hh2 = []
ww2 = []
zz = []
PP = []

if model== 'pow_exp': min, max = (0.1, 1.1)
elif model== 'cosh':  min, max = (0.5, 1.6)
else:                 min, max = (-0.1, 0.1)

steps    = 4.
step     = (max-min)/steps

mymap    = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['red','blue'])
Z        = [[0,0],[0,0]]
levels   = np.arange(min, max+step, step)
CS3      = plt.contourf(Z, levels, cmap=mymap)

beta, mu, alpha, ilam = 3.14, 3.14, 3.14, 3.14
#Quintess
for i in np.arange(min, max, step):

    beta, mu, ilam = 0, 2.0, -0.5
    T= PhiCosmology(beta=beta, mu=mu, alpha=alpha, ilam= ilam, eps=1, curv=i)


    y1= [T.Omegaphi(la) for la in lna]
    y2= [T.Omegak(la) for la in lna]
    y3= [1- T.Omegaphi(la)-T.Omegak(la) for la in lna]
    #y1=[T.HIOverrd(z)*z/fixer(z) for z in zl]
    #y2=[T.DaOverrd(z)/fixer(z)   for z in zl]
    #y1=[T.w_de(1./(1+z))     for z in zl]
    #y2=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
    ww.append(y1)
    hh.append(y2)
    ww2.append(y3)
    PP.append(i)
    zz.append(lna)



if False: #not(model=='exp'):
    for i in np.arange(min, max, step):

        beta, mu, ilam = 0, 1, -0.1
        T= PhiCosmology(beta=beta, mu=mu, alpha=alpha, ilam= ilam, eps=1, curv=i)

        y1=[T.HIOverrd(z)*z/fixer(z) for z in zl]
        y2=[T.DaOverrd(z)/fixer(z)   for z in zl]
        #y1=[T.w_de(1./(1+z))     for z in zl]
        #y1=[T.distance_modulus(z) for z in zl ]
        #y2=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
        ww.append(y1)
        hh.append(y2)
        PP.append(i)
        zz.append(zl)
"""
#Phantom
for i in np.arange(min, max, step):
    if   model== 'exp':       T= PhiCosmology(beta=i, mu=0, alpha=1., eps=1)
    elif model== 'pow':       T= PhiCosmology(beta=0, mu=i, ilam=-0.75, eps=-1)
    elif model== 'exp_pow2':  T= PhiCosmology(beta=-i, mu=0, alpha=2, ilam=1.1, eps=-1)
    elif model== 'pow_exp':   T= PhiCosmology(beta=i, mu=1, alpha=1, ilam=0.5)
    elif model== 'exp_pow_a': T= PhiCosmology(beta=1.5, mu=0, alpha=i, ilam=-1.0)
    elif model== 'cosh': T= PhiCosmology(mu='st', beta=-1, alpha=i, ilam=-0.1)

    y12=[T.w_de(1./(1+z))     for z in zl]
    y22=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
    ww2.append(y12)
    hh2.append(y22)

if not(model=='exp'):
    for i in np.arange(min, max, step):
        if model== 'pow':         T= PhiCosmology(beta=0, mu=-i, ilam=0.75, eps=-1)
        elif model== 'exp_pow2':  T= PhiCosmology(beta=i, mu=0, alpha=2, ilam=-0.25, eps=-1)
        elif model== 'pow_exp':   T= PhiCosmology(beta=i, mu=-1, alpha=1, ilam=-0.75)
        elif model== 'exp_pow_a': T= PhiCosmology(beta=1., mu=0, alpha=i, ilam=-1.0)
        elif model== 'cosh':      T= PhiCosmology(mu='st', beta=1, alpha=i, ilam=1.0)

        y12=[T.w_de(1./(1+z))     for z in zl]
        y22=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
        ww2.append(y12)
        hh2.append(y22)
"""

#T2 = LCDMCosmology()



fig, ax1  = plt.subplots(1, figsize=(7,6))
fig.suptitle('$\phi^2$', fontsize=17,  y=0.95)

for i, (x,w,w2,w3, z) in enumerate(zip(zz, ww, ww2, hh, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls = '-' if i <=steps-1 else '--'

        l1,  = ax1.plot(x, w , color=(r,g,b), linestyle='-',
                        label='$\Omega_\phi$' if i==steps-1 else None, lw=1.)
        l1,  = ax1.plot(x, w2 , color=(r,g,b), linestyle='--',
                        label='$\Omega_m$' if i==steps-1 else None, lw=1.)
        l1,  = ax1.plot(x, w3 , color=(r,g,b), linestyle='-.',
                        label='$\Omega_k$' if i==steps-1 else None, lw=1.)
ax1.plot(lna, x1, 'k-')
ax1.plot(lna, np.ones(len(x1))-x1, 'k-')
cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.78])
cbar = pylab.colorbar(CS3, cax=cbaxes)  #fraction=0.046, pad=0.04)
cbar.set_label('$\Omega_k$', rotation=0, fontsize=20, labelpad=-20)
cbar.ax.tick_params(labelsize=15)
ax1.legend(loc = 'best', frameon=False, fontsize=25)
ax1.set_ylabel('$\Omega(z)$', fontsize=20, labelpad = -5)
ax1.set_xlabel("$\ln a$", fontsize=20)
"""
for i, (x,w,z) in enumerate(zip(zz, ww2, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls = '-' if i <=steps-1 else '--'

        s = steps-1 if model =='exp' else 2*steps-1
        if model =='exp': ls = '--'
        lb ={ 'exp': '$V=e^{\\beta \phi}$', 'exp_pow2':'$V=e^{\\beta \phi^2}$',
            'pow':'$V=\phi^{-\\mu}$', 'pow_exp': '$V=\phi^{1} e^{-\\beta \phi}$',
            'exp_pow_a':'$V=e^{\\beta \phi^\\alpha}, \\beta=1.0$', 'pow2_exp_pow2':'test',
              'cosh':'$V=\cos(\\alpha \phi)+1$'}
        l1,  = ax1.plot(x, w , color=(r,g,b), linestyle=ls, lw=1.5,
                        label=lb[model] if i==s else None)
#ax1.set_ylim(-2, 0)
#ax1.plot(zl, -np.ones(len(zl)) , color='k', linestyle=':', lw=2)
ax1.legend(loc = 'best', frameon=False, fontsize=25)
ax1.set_ylabel('$w(z)$', fontsize=20, labelpad = -5)
"""
ax1.grid()
"""
for i, (x,w,z) in enumerate(zip(zz, hh, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls = '-' if i <=steps-1 else '--'
        l2, = ax2.plot(x, w , color=(r,g,b) , linestyle=ls, lw=1.5)


for i, (x,w,z) in enumerate(zip(zz, hh2, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls = '-' if i <=steps-1 else '--'
        if model =='exp': ls = '--'
        l2,  = ax2.plot(x, w , color=(r,g,b), linestyle=ls, lw=1.5)

#ax2.plot(zl, x1 , color='k', linestyle=':', lw=2)
ax2.set_ylabel('$H(z)/(1+z)$', fontsize=20)
#ax2.set_ylim(55, 75)
ax2.grid()

cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.78])
cbar = pylab.colorbar(CS3, cax=cbaxes)  #fraction=0.046, pad=0.04)
#cbar.set_label('$\\lambda_i$', rotation=0, fontsize=18, labelpad=-10)

lab = {'exp': '$\\beta$','exp_pow2':'$\\beta$','pow':'$\\mu$','pow_exp':'$\\beta$',
       'exp_pow_a':'$\\alpha$', 'pow2_exp_pow2':'$\\beta$', 'cosh': '$\\alpha$'}
cbar.set_label(lab[model], rotation=0, fontsize=20, labelpad=-10)
cbar.ax.tick_params(labelsize=15)
ax2.set_xlabel("$z$", fontsize=20)
"""
plt.savefig('Omegas.pdf', bbox_inches='tight')
plt.show()
