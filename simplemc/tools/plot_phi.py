#!/usr/bin/env python

#Error bars mainly got from  table 3,
#http://arxiv.org/pdf/1108.2635v1.pdf
#and MGS paper



from simplemc.models.PhiCDMCosmology import PhiCosmology
from simplemc.models.LCDMCosmology import LCDMCosmology
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab
import sys


model = sys.argv[1]
models= ['exp', 'pow' , 'exp_pow2', 'pow_exp', 'exp_pow_a',
         'pow2_exp2', 'cosh']
if model not in models:
    sys.exit('models allowed: %s'%models)


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

if model== 'pow_exp': min, max = (0.1, 1.1)
elif model== 'cosh':  min, max = (0.5, 1.6)
else:                 min, max = (0.5, 2.1)

steps    = 4.
step     = (max-min)/steps

mymap    = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['red','blue'])
Z        = [[0,0],[0,0]]
levels   = np.arange(min, max+step, step)
CS3      = plt.contourf(Z, levels, cmap=mymap)

beta, mu, alpha, ilam = 3.14, 3.14, 3.14, 3.14
#Quintess
for i in np.arange(min, max, step):

    if   model=='exp':      beta, mu, alpha= -i, 0, 1
    elif model=='pow':      beta, mu, ilam = 0, i, -0.75
    elif model=='exp_pow2': beta, mu, alpha, ilam = i, 0, 2, -0.50
    elif model== 'pow_exp': beta, mu, alpha, ilam = -i, -1, 1, 0.5
    elif model== 'exp_pow_a':beta, mu, alpha, ilam = -1.2, 0, i, 0.75
    elif model== 'cosh':    beta, mu, alpha, ilam = -1, 'st', i, -0.1
    #elif model== 'pow2_exp2':beta, mu, alpha, ilam = -i, 2, 2, -0.25
    T= PhiCosmology(beta=beta, mu=mu, alpha=alpha, ilam= ilam, curv=0)

    """
    elif model== 'cosh':      T= PhiCosmology(mu='st', beta=-1, alpha=i, ilam=-0.1)
    """
    y1=[T.w_de(1./(1+z))     for z in zl]
    y2=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
    ww.append(y1)
    hh.append(y2)
    PP.append(i)
    zz.append(zl)

if not(model=='exp'):
    for i in np.arange(min, max, step):
        if model=='pow':         beta, mu, ilam = 0, -i, 0.75
        elif model=='exp_pow2':  beta, mu, alpha, ilam = -i, 0, 2, 0.15
        elif model== 'pow_exp':  beta, mu, alpha, ilam = -i, 1, 1, -0.5
        elif model== 'exp_pow_a':beta, mu, alpha, ilam = -0.8, 0, i, 0.75
        elif model== 'cosh':     beta, mu, alpha, ilam = -1, 'st', i, 0.1
        #elif model== 'pow2_exp2':beta, mu, alpha, ilam = i, 2, 2, -0.75
        T= PhiCosmology(beta=beta, mu=mu, alpha=alpha, ilam= ilam, curv=0)


        #elif model== 'cosh':      T= PhiCosmology(mu='st', beta=1, alpha=i, ilam=1)

        y1=[T.w_de(1./(1+z))     for z in zl]
        y2=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
        ww.append(y1)
        hh.append(y2)
        PP.append(i)
        zz.append(zl)


#Phantom
for i in np.arange(min, max, step):
    if   model=='exp':      beta, mu, alpha= i, 0, 1
    elif model=='pow':      beta, mu, ilam = 0, i, 0.75
    elif model=='exp_pow2': beta, mu, alpha, ilam = i, 0, 2, 0.15
    elif model== 'pow_exp': beta, mu, alpha, ilam = -i, -1, 1, -0.75
    elif model== 'exp_pow_a':beta, mu, alpha, ilam = -1.2, 0, i, -0.45
    elif model== 'cosh':    beta, mu, alpha, ilam = 1, 'st', i, -1.5
    #elif model== 'pow2_exp2':beta, mu, alpha, ilam = -i, 2, 2, 0.75

    T= PhiCosmology(beta=beta, mu=mu, alpha=alpha, ilam= ilam, curv=0)


    #elif model== 'cosh': T= PhiCosmology(mu='st', beta=-1, alpha=i, ilam=-0.1)

    y12=[T.w_de(1./(1+z))     for z in zl]
    y22=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
    ww2.append(y12)
    hh2.append(y22)

if not(model=='exp'):
    for i in np.arange(min, max, step):
        if   model=='exp':      beta, mu, alpha= i, 0, 1
        elif model=='pow':      beta, mu, ilam = 0, -i, -0.75
        elif model=='exp_pow2': beta, mu, alpha, ilam = -i, 0, 2, -0.50
        elif model== 'pow_exp': beta, mu, alpha, ilam = -i, 1, 1, 0.75
        elif model== 'exp_pow_a':beta, mu, alpha, ilam = -0.8, 0, i, -0.75
        elif model== 'cosh':     beta, mu, alpha, ilam = 1, 'st', i, 1.5
        T= PhiCosmology(beta=beta, mu=mu, alpha=alpha, ilam= ilam, curv=0)

        #elif model== 'exp_pow_a': T= PhiCosmology(beta=1., mu=0, alpha=i, ilam=-1.0)
        #elif model== 'cosh':      T= PhiCosmology(mu='st', beta=1, alpha=i, ilam=1.0)

        y12=[T.w_de(1./(1+z))     for z in zl]
        y22=[67.4*np.sqrt(T.RHSquared_a(1./(1+z)))/(1+z)     for z in zl]
        ww2.append(y12)
        hh2.append(y22)


fig, (ax1, ax2)  = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0, "height_ratios":[1,0.5]},
                                figsize=(7,9))
if model== 'exp_pow_a':
    fig.suptitle('$V=e^{\\beta \phi^\\alpha}$', fontsize=25,  y=0.88)

for i, (x,w,z) in enumerate(zip(zz, ww, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls  = '-' if i <=steps-1 else '--'
        lb2 = None
        if model== 'exp_pow_a':
            ls = '-' if i <=steps-1 else '--'
            if i== 2*steps-1:
              lb2 = '$\\beta=-1.0$'

        lb ={ 'exp':'$V=e^{-\\beta \phi}$', 'exp_pow2':'$V=e^{\\beta \phi^2}$',
            'pow':'$V=\phi^\\mu$', 'pow_exp': '$V=\phi^{-1} e^{-\\beta \phi}$',
            'exp_pow_a':'$\\beta=-1.5$', 'pow2_exp2':'test',
              'cosh':'$V=\cosh(\\alpha \phi)-1$'}

        l1,  = ax1.plot(x, w , color=(r,g,b), linestyle=ls,
                        label=lb[model] if i==steps-1 else lb2, lw=1.5)


for i, (x,w,z) in enumerate(zip(zz, ww2, PP)):
        b    = (float(z)-min)/(max-min)
        g, r = 0, 1-b
        ls = '-' if i <=steps-1 else '--'

        s = steps-1 if model =='exp' else 2*steps-1
        lb2 = None
        if model== 'exp_pow_a':
            ls = '-.' if i <=steps-1 else ':'
            s= steps-1
            if i== 2*steps-1:
              lb2 = '$\\beta=1.0$'


        if model =='exp': ls = '--'
        lb ={ 'exp': '$V=e^{\\beta \phi}$', 'exp_pow2':'$V=e^{-\\beta \phi^2}$',
            'pow':'$V=\phi^{-\\mu}$', 'pow_exp': '$V=\phi^{1} e^{-\\beta \phi}$',
            'exp_pow_a':'$\\beta=1.5$', 'pow2_exp2':'test',
              'cosh':'$V=\cos(\\alpha \phi)+1$'}
        l1,  = ax1.plot(x, w , color=(r,g,b), linestyle=ls, lw=1.5,
                        label=lb[model] if i==s else lb2)
#ax1.set_ylim(-2, 0)
ax1.plot(zl, -np.ones(len(zl)) , color='k', linestyle=':', lw=2)
ax1.legend(loc = 'best', frameon=False, fontsize=20)
ax1.set_ylabel('$w(z)$', fontsize=20, labelpad = -5)
ax1.grid()

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

ax2.plot(zl, x1 , color='k', linestyle=':', lw=2)
ax2.set_ylabel('$H(z)/(1+z)$', fontsize=20)
ax2.set_ylim(55, 75)


rd_fid_DR12 = 147.78

zLyaA = 2.37
zLyaC = 2.35
zCombBAO1 = 0.38
zCombBAO2 = 0.51
zCombBAO3 = 0.61
fact = (300000./rd_fid_DR12)

def ersys(x, y):
    return np.sqrt(x**2 + y**2)

#https://arxiv.org/abs/1904.03430
ax2.errorbar(zLyaC,  fact/9.20/(1+zLyaC),        yerr=0.36*fact/(1+zLyaC)/9.20**2,
              color='green', fmt='s', markersize=6, elinewidth=1.5)
#https://arxiv.org/abs/1904.03430
ax2.errorbar(zLyaA,  fact/8.86/(1+zLyaA),        yerr=0.29*fact/(1+zLyaA)/8.86**2,
              color='green', fmt='s', markersize=6, elinewidth=1.5)

ax2.errorbar(zCombBAO1,  81.21/(1+zCombBAO1),       yerr=ersys(2.17, 0.97)/(1+zCombBAO1),
              color='green', fmt='d', markersize=6, elinewidth=1.5)
ax2.errorbar(zCombBAO2,  90.90/(1+zCombBAO2),       yerr=ersys(2.07, 1.08)/(1+zCombBAO1),
              color='green', fmt='d', markersize=6, elinewidth=1.5)
ax2.errorbar(zCombBAO3,  98.96/(1+zCombBAO3),       yerr=ersys(2.21, 1.18)/(1+zCombBAO1),
              color='green', fmt='d', markersize=6, elinewidth=1.5)


ax2.grid()

cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.78])
cbar = pylab.colorbar(CS3, cax=cbaxes)  #fraction=0.046, pad=0.04)
#cbar.set_label('$\\lambda_i$', rotation=0, fontsize=18, labelpad=-10)

lab = {'exp': '$\\beta$','exp_pow2':'$\\beta$','pow':'$\\mu$','pow_exp':'$\\beta$',
       'exp_pow_a':'$\\alpha$', 'pow2_exp2':'$\\beta$', 'cosh': '$\\alpha$'}
cbar.set_label(lab[model], rotation=0, fontsize=20, labelpad=-10, y=0.65)
cbar.ax.tick_params(labelsize=15)
ax2.set_xlabel("$z$", fontsize=20)

plt.savefig(model + '.pdf', bbox_inches='tight')
plt.show()
