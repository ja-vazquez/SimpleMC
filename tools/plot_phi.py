#!/usr/bin/env python

#Error bars mainly got from  table 3,
#http://arxiv.org/pdf/1108.2635v1.pdf
#and MGS paper

import sys
sys.path = ["Models", "Cosmo"] + sys.path

from PhiCDMCosmology import PhiCosmology
from LCDMCosmology import LCDMCosmology
from ParamDefs import *
import matplotlib.pyplot as plt
import pylab
import numpy as np
import matplotlib as mpl


model = 'exp' #pow, exp, exp_beta, pow2
T = PhiCosmology()

if model=='pow':
    T.poten  = 'pow'
    T.phi0  = 1.5
    name    = 'phi_pow'
elif model=='pow2':
    #change to varylam
    T.poten ='pow2'
    name   ='phi_pow2'
elif model=='exp':
    T.poten  = 'exp'
    name    = "phi_exp"
elif model=='exp_beta':
    T.poten  = 'exp'
    name    = "phi_exp_beta"

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


def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
                txt.set_color(line.get_color())

#-------------------

#Planck best fit cosmology
T2=LCDMCosmology()
x1=[67.4*np.sqrt(T2.RHSquared_a(1./(1+z)))/(1+z)      for z in zl]

hh = []
ww = []
hh2 = []
ww2 = []
zz = []
PP = []

min, max = (0.1, 2.0)
steps    = 5
step     = (max-min)/steps

mymap    = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['red','blue'])
Z        = [[0,0],[0,0]]
levels   = np.arange(min, max+step, step)
CS3      = plt.contourf(Z, levels, cmap=mymap)

ilam_  = plam_par
ialp_  = palp_par
ibeta_ = pbeta_par
if model=='exp':      T.alpha=1.
if model=='exp_beta': T.beta =1.5
if model=='pow2':     T.alpha=2.
for i in np.arange(min, max, step):
    if model=='pow' or model=='exp_beta':
        ialp_.setValue(i)
        T.updateParams([ialp_])
    elif model=='exp':
        ibeta_.setValue(i)
        T.updateParams([ibeta_])
    elif model=='pow2':
        ilam_.setValue(i)
        T.updateParams([ilam_])
    y1=[T.w_de(1./(1+z))     for z in zl]
    y2=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
    ww.append(y1)
    hh.append(y2)
    PP.append(i)
    zz.append(zl)
if model=='pow2':     T.alpha=-2.
if model =='pow' or model=='exp_beta' or model=='pow2':
    if model=='exp_beta': T.beta =1.
    for i in np.arange(min, max, step):
        if model=='pow':
            ialp_.setValue(-i)
            T.updateParams([ialp_])
        elif model=='exp_beta':
            ialp_.setValue(i)
            T.updateParams([ialp_])
        elif model=='pow2':
            ilam_.setValue(i)
            T.updateParams([ilam_])

        y1=[T.w_de(1./(1+z))     for z in zl]
        y2=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
        ww.append(y1)
        hh.append(y2)
        PP.append(i)
        zz.append(zl)

T.qp=-1
if model=='exp':      T.alpha=1.
if model=='exp_beta': T.beta =1.5
if model=='pow2':     T.alpha=2.
for i in np.arange(min, max, step):
    if model=='pow' or model=='exp_beta':
        ialp_.setValue(i)
        T.updateParams([ialp_])
    elif model=='exp':
        ibeta_.setValue(i)
        T.updateParams([ibeta_])
    elif model=='pow2':
        ilam_.setValue(i)
        T.updateParams([ilam_])
    y12=[T.w_de(1./(1+z))     for z in zl]
    y22=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
    ww2.append(y12)
    hh2.append(y22)
if model=='pow2':     T.alpha=-2.
if model =='pow' or model=='exp_beta' or model=='pow2':
    if model=='exp_beta': T.beta =1.
    for i in np.arange(min, max, step):
        if model=='pow':
            ialp_.setValue(-i)
            T.updateParams([ialp_])
        elif model=='exp_beta':
            ialp_.setValue(i)
            T.updateParams([ialp_])
        elif model=='pow2':
            ilam_.setValue(i)
            T.updateParams([ilam_])

        y12=[T.w_de(1./(1+z))     for z in zl]
        y22=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
        ww2.append(y12)
        hh2.append(y22)


fig, (ax1, ax2)  = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0}, figsize=(7,10))
#fig.suptitle(fname, fontsize=17,  y=0.95)

for i, (x,w,z) in enumerate(zip(zz, ww, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls = '-' if i <=5 else '--'
        if model=='pow':        lb= '$V=\phi^\\alpha$'
        elif model=='pow2':     lb= '$V=\phi^2$'
        elif model=='exp':      lb= '$V=e^{\\beta \phi}$'
        elif model=='exp_beta': lb= '$V=e^{\\beta \phi^\\alpha}, \\beta=1.5$'
        l1,  = ax1.plot(x, w , color=(r,g,b), linestyle=ls ,label=lb if i==steps-1 else None)

for i, (x,w,z) in enumerate(zip(zz, ww2, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls = '-' if i <=5 else '--'
        if model=='pow':        lb= '$V=\phi^{-\\alpha}$'
        elif model=='pow2':     lb= '$V=\phi^{-2}$'
        elif model=='exp':      lb=  None
        elif model=='exp_beta': lb='$V=e^{\\beta \phi^\\alpha}, \\beta=1.0$'
        l1,  = ax1.plot(x, w , color=(r,g,b), linestyle=ls, lw=1, label=lb if i==2*steps-1 else None)

ax1.legend(loc = 'best', frameon=False, fontsize=20)
ax1.set_ylabel('$w(z)$', fontsize=20, labelpad = -10)
ax1.grid()

for i, (x,w,z) in enumerate(zip(zz, hh, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls = '-' if i <=5 else '--'
        l2, = ax2.plot(x, w , color=(r,g,b) , linestyle=ls)
for i, (x,w,z) in enumerate(zip(zz, hh2, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls = '-' if i <=5 else '--'
        l2,  = ax2.plot(x, w , color=(r,g,b), linestyle=ls)

plt.plot(zl, x1 , color='k', linestyle=':', lw=3)
ax2.set_ylabel('$H(z)/(1+z)$', fontsize=20)
ax2.set_ylim(55, 75)
ax2.grid()

cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.78])
cbar = pylab.colorbar(CS3, cax=cbaxes)  #fraction=0.046, pad=0.04)
#cbar.set_label('$\\lambda_i$', rotation=0, fontsize=18, labelpad=-10)
if model=='pow':        lab= '\\alpha'
elif model=='pow2':     lab= '\lambda_i'
elif model=='exp':      lab= '$\\beta$'
elif model=='exp_beta': lab= '\\alpha'
cbar.set_label(lab, rotation=0, fontsize=18, labelpad=-10)
cbar.ax.tick_params(labelsize=12)
ax2.set_xlabel("$z$", fontsize=18)

#plt.savefig(name + '.pdf', bbox_inches='tight')
plt.show()
