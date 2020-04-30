#!/usr/bin/env python

#Error bars mainly got from  table 3,
#http://arxiv.org/pdf/1108.2635v1.pdf
#and MGS paper

import sys
sys.path = ["Models"] + sys.path

from QuintomCosmology import QuintomCosmology
from LCDMCosmology import LCDMCosmology
from ParamDefs import *
import matplotlib.pyplot as plt
import matplotlib.ticker
import pylab
import math as N
import numpy as np
import matplotlib as mpl
from matplotlib import cm


beta  = 0

#fname = 'Quintessence'
#fname = 'Phantom'
#fname = 'Quintmphi'
#mphan = 1.2
#fname = 'Quintmphan'
#mphi  = 1.2
#fname = 'Quintomcopphi'
#mphan = 1.2
#beta  = 4.0
#fname = 'Quintomcopphan'
#mphi  = 1.2
#beta  = 6.0
fname = 'Quintomcopbeta'
mphi  = 1.2
mphan = 0.1





addlabel  = False
addlegend = False

plaw=0.5
eboss=False


params1 = {'backend': 'pdf',
               'axes.labelsize':  18,
                #'text.fontsize': 18,
               'xtick.labelsize': 18,
               'ytick.labelsize': 18,
                #'legend.draw_frame': False,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}
pylab.rcParams.update(params1)


zLOWZ  = 0.32 
zCMASS = 0.57
zLyaA  = 2.34-0.04
zLyaC  = 2.36+0.04

z6dFGS   = 0.106
zMGS     = 0.15
zSDSS1   = 0.2
zSDSS2   = 0.35
zWiggleZ1=0.44
zWiggleZ2= 0.6
zWiggleZ3= 0.73

z_CMB = 1090.43

zCombBAO1 = 0.38
zCombBAO2 = 0.51
zCombBAO3 = 0.61

rd_EHtoCAMB =153.19/149.28
rd_fid_DR12 = 147.78
rd_fid_DR7  = 151.84  

zl=np.arange(0, 3, 0.05)

fmt1 = 'o'
fmt2 = 'o'
empty1= False
empty2= False
alpha= 1.0


def fixer(z):
    if plaw>0:
        return z**plaw
    else:
        return log(1.+z)

#============Functions -------


### Plotting -  Error bars
def plot_errorbar(z,val, yerr=0, color=0, fmt=0, markersize=0,label=None, empty=True, alpha=1, ax=None):
    if empty:
        mfc='white'
        lw=1
    else:
        mfc=color
        lw=1
    color = 'black'

    ax.errorbar(z,val/fixer(z), yerr=yerr/fixer(z),
                color='blue', marker='o', ls='None',
                elinewidth =2, capsize=3, capthick = 1, alpha=1, markersize=4)

    if addlabel:
        if label>0:
            if (mfc=='white'):
                ax.plot ([],[],fmt,color='black',label=label,markersize=markersize,markerfacecolor=mfc)
            else:
                ax.plot ([],[],fmt,color='black',label=label,markersize=markersize)

def ersys(x, y):
    return np.sqrt(x**2 + y**2)


def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
                txt.set_color(line.get_color())

#-------------------



#Planck best fit cosmology
#T2=LCDMCosmology()
#x1=[67.4*np.sqrt(T2.RHSquared_a(1./(1+z)))     for z in zl]
#x2=[T2.HIOverrd(z)*z/fixer(z) for z in zl]
#x3=[T2.DaOverrd(z)/fixer(z)   for z in zl]
#PLK-15
#T=LCDMCosmology(Obh2=0.02225,Om=0.3156,h=0.6727)


if fname == 'Quintessence': T = QuintomCosmology(varymquin=True)
if fname == 'Phantom':      T = QuintomCosmology(varymphan=True)
if fname == 'Quintmphi':  T = QuintomCosmology(varymquin=True, varymphan=True)
if fname == 'Quintmphan': T = QuintomCosmology(varymquin=False, varymphan=True)
if fname == 'Quintomcopphi': T = QuintomCosmology(varymquin=True, varybeta=True)
if fname == 'Quintomcopphan': T = QuintomCosmology(varymphan=True, varybeta=True)
if fname == 'Quintomcopbeta': T = QuintomCosmology(varymquin=True, varymphan=True, varybeta=True)

hh = []
ww = []
da = []
dh = []
zz = []
PP = []

if fname == 'Quintomcopbeta':
    min, max = (4., 8.)
else:
    min, max = (0.1, 2.5)
if beta < 0:  min, max = (-4, -1.)

steps    = 2
step     = (max-min)/steps

mymap    = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['red','green'])
Z        = [[0,0],[0,0]]
levels   = np.arange(min, max+step, step)
CS3      = plt.contourf(Z, levels, cmap=mymap)

mquin_ = mquin_par
mphan_ = mphan_par
beta_  = beta_par

for i in np.arange(min, max, step):
    if fname == 'Quintessence':
        mquin_.setValue(i)
        T.updateParams([mquin_])
        name = fname
        mlabel = '$m_\phi$'
    if fname == 'Phantom':
        mphan_.setValue(i)
        T.updateParams([mphan_])
        name = fname
        mlabel = '$m_\psi$'
    if fname == 'Quintmphi':
        mquin_.setValue(i)
        mphan_.setValue(mphan)
        T.updateParams([mquin_, mphan_])
        name = 'Quintom, $m_{\psi}$=%0.1f'%mphan
        mlabel = '$m_\phi$'
    if fname == 'Quintmphan':
        mphan_.setValue(i)
        mquin_.setValue(mphi)
        T.updateParams([mquin_, mphan_])
        name = 'Quintom, $m_{\phi}$=%0.1f'%mphi
        mlabel = '$m_\psi$'
    if fname == 'Quintomcopphi':
        mquin_.setValue(i)
        mphan_.setValue(mphan)
        beta_.setValue(beta)
        T.updateParams([mquin_, mphan_, beta_])
        name = 'Quintom, $m_{\psi}$=%0.1f, $\\beta=%0.1f$'%(mphan, beta)
        mlabel = '$m_\phi$'
    if fname == 'Quintomcopphan':
        mquin_.setValue(mphi)
        mphan_.setValue(i)
        beta_.setValue(beta)
        T.updateParams([mquin_, mphan_, beta_])
        name = 'Quintom, $m_{\phi}$=%0.1f, $\\beta=%0.1f$'%(mphi, beta)
        mlabel = '$m_\psi$'
    if fname == 'Quintomcopbeta':
        mquin_.setValue(mphi)
        mphan_.setValue(mphan)
        beta_.setValue(i)
        T.updateParams([mquin_, mphan_, beta_])
        name = 'Quintom, $m_{\phi}$=%0.1f, $m_{\psi}$=%0.1f'%(mphi, mphan)
        mlabel = '$\\beta$'

    y1=[T.w_de(1./(1+z))         for z in zl]
    y2=[T.Hubble_a(1./(1+z))     for z in zl]
    y3=[T.HIOverrd(z)*z/fixer(z) for z in zl]
    y4=[T.DaOverrd(z)/fixer(z)   for z in zl]


    ww.append(y1)
    hh.append(y2)
    dh.append(y3)
    da.append(y4)
    PP.append(i)
    zz.append(zl)


fig, (ax1, ax2, ax3, ax4)  = plt.subplots(4, sharex=True, gridspec_kw={'hspace': 0}, figsize=(7,10))



fig.suptitle(name, fontsize=17,  y=0.95)

for x,w,z in zip(zz, ww, PP):
        g    = (float(z)-min)/(max-min)
        b, r = 0, 1-g
        l1, = ax1.plot(x, w , color=(r,g,b)) #, linestyle='-.')
if (fname == 'Quintessence') or (fname =='Quintomcopphi'): ax1.set_ylabel('$w(z)$', fontsize=20)
ax1.axhline(y=-1, color='k', linestyle='--')
if beta <0 : ax1.set_ylim(-3, 0.)

#ax1.set_xscale('log')
for x,w,z in zip(zz, hh, PP):
        g    = (float(z)-min)/(max-min)
        b, r = 0, 1-g
        l2, = ax2.plot(x, w , color=(r,g,b)) #, linestyle='-.')
#ax2.plot(zl, x1 , color='k', linestyle='--')
dataHz = np.loadtxt('data/Hz_all.dat')
redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
ax2.errorbar(redshifts, obs, errors, xerr=None,
                 color='blue', marker='o', ls='None',
                elinewidth =2, capsize=3, capthick = 1, alpha=1, markersize=4)
#ax2.legend(loc='lower right', frameon=False)
if (fname == 'Quintessence') or (fname =='Quintomcopphi'): ax2.set_ylabel('$H(z)$', fontsize=20)


cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.78])
cbar = pylab.colorbar(CS3, cax=cbaxes) #,  pad=0.04) #fraction=0.046, pad=0.04)
cbar.set_label(mlabel, rotation=0, fontsize=18, labelpad=-10)
cbar.ax.tick_params(labelsize=12)

#ax2.set_xscale('log')
for x,w,z in zip(zz, dh, PP):
        g    = (float(z)-min)/(max-min)
        b, r = 0, 1-g
        l3, = ax3.plot(x, w , color=(r,g,b)) #, linestyle='-.')
if (fname == 'Quintessence') or (fname == 'Quintomcopphi'): ax3.set_ylabel("${\\rm zD_H(z)}/r_d\\sqrt{z}$")
#ax3.plot(zl, x2 , color='k', linestyle='--')


#ax4.set_xscale('log')
for x,w,z in zip(zz, da, PP):
        g    = (float(z)-min)/(max-min)
        b, r = 0, 1-g
        l4, = ax4.plot(x, w , color=(r,g,b)) #, linestyle='-.')
ax4.set_xlim(0.05, 3)
if (fname == 'Quintessence') or (fname == 'Quintomcopphi'): ax4.set_ylabel("${\\rm D_M(z)}/r_d\\sqrt{z}$")
#ax4.plot(zl, x3, color='k', linestyle='--')

if addlegend:
    if plaw>0:
        if eboss:
            legend1=pylab.legend([l1,l2,l3],["$D_%s(z)/r_d\\sqrt{z}$"%st for st in ['M','V','H']],loc="lower center", bbox_to_anchor = (0.7, 0.0))
        else:
            legend1=pylab.legend([l1,l2,l3],["$%sD_%s(z)/r_d\\sqrt{z}$"%st for st in [('','M'),('','V'),('z','H')]],loc="lower center")
    else:
        legend1=pylab.legend([l1,l2,l3],["$D_%s(z)/r_d\log(1+z)$"%st for st in ['A','V','H']],loc="lower center")

    pylab.gca().add_artist(legend1)
    legend1.draw_frame(False)
    color_legend(legend1)




pylab.legend(loc="lower right")

fact = (300000./rd_fid_DR12)
"""
plot_errorbar(z6dFGS,    2.97*rd_EHtoCAMB,   yerr=rd_EHtoCAMB*0.015/0.336**2,  color ='blue', fmt='o',
              markersize=4, empty=empty2, label="$\\rm{6dFGS}$", alpha=alpha, ax=ax3)
plot_errorbar(zMGS,      4.464,    yerr=0.168,               color ='blue', fmt='o',
              markersize=4, label="$\\rm{SDSS\ MGS}$", empty=empty2, alpha=alpha, ax=ax3)

plot_errorbar(zSDSS1,    5.2493*rd_EHtoCAMB, yerr=rd_EHtoCAMB*0.0061/0.1905**2,color ='blue', fmt=fmt1,
              markersize=4, empty=empty1, alpha=alpha, ax=ax3)
plot_errorbar(zSDSS2,    1348./rd_fid_DR12, yerr=26./rd_fid_DR12 ,color ='blue', fmt=fmt1,
              markersize=4, label="$\\rm{SDSS\ DR7}$", empty=empty1, alpha=alpha, ax=ax3)

plot_errorbar(zWiggleZ1, 1695./rd_fid_DR12 ,yerr=82./rd_fid_DR12 ,color ='blue', fmt=fmt2,
              markersize=4, label="$\\rm{WiggleZ}$", empty=empty1, alpha=alpha, ax=ax3)
plot_errorbar(zWiggleZ2, 2194./rd_fid_DR12 ,yerr=100./rd_fid_DR12 ,color ='blue', fmt=fmt2,
              markersize=4, empty=empty1, alpha=alpha, ax=ax3)
plot_errorbar(zWiggleZ3, 2486./rd_fid_DR12 ,yerr=85./rd_fid_DR12 ,color ='blue', fmt=fmt2,
              markersize=4, empty=empty1, alpha=alpha, ax=ax3)
"""

if eboss:
    plot_errorbar(zLOWZ,     8.467,              yerr=0.17,        color ='blue', fmt='o',
                  markersize=4, label="$\\mathrm{LOWZ}$",empty=False, ax=ax3)
else:
    pass	


plot_errorbar(zCombBAO1,  1512.4/rd_fid_DR12,     yerr=ersys(22.5, 11.0)/rd_fid_DR12,       
		color ='red', fmt='o', markersize=4, empty=empty2, label="$\\rm{BOSS\ Galaxy\ DR12}$", ax=ax4)
plot_errorbar(zCombBAO2,  1975.2/rd_fid_DR12,     yerr=ersys(26.6, 14.1)/rd_fid_DR12,
		color ='red', fmt='o', markersize=4, empty=empty2, ax=ax4)
plot_errorbar(zCombBAO3,  2306.7/rd_fid_DR12,  	 yerr=ersys(33.2, 16.7)/rd_fid_DR12, 
		color ='red', fmt='o', markersize=4, empty=empty2, ax=ax4)


plot_errorbar(zCombBAO1,  fact*zCombBAO1/81.21,       yerr=fact*zCombBAO1*ersys(2.17, 0.97)/(81.21)**2,  
		color ='green', fmt='o', markersize=4, empty=empty2, ax=ax3)
plot_errorbar(zCombBAO2,  fact*zCombBAO2/90.90,       yerr=fact*zCombBAO2*ersys(2.07, 1.08)/(90.90)**2,    
		color ='green', fmt='o', markersize=4, empty=empty2, ax=ax3)
plot_errorbar(zCombBAO3,  fact*zCombBAO3/98.96,       yerr=fact*zCombBAO3*ersys(2.21, 1.18)/(98.96)**2, 
		color ='green', fmt='o', markersize=4, empty=empty2, ax=ax3)


plot_errorbar(zLyaA,  11.28*(1+zLyaA),  yerr=0.65*(1+ zLyaA),  color ='red', fmt='o', markersize=4,
              label="$\\rm{BOSS}\ \\mathrm{Ly}\\alpha\\mbox{-}\\rm{auto}\ \\rm{DR11}$", empty=empty2, ax=ax4)
plot_errorbar(zLyaA,  9.18*zLyaA,       yerr=0.28*zLyaA,       color ='green', fmt='o', markersize=4,empty=empty2, ax=ax3)

plot_errorbar(zLyaC,  10.8*(1+zLyaC),   yerr=0.4*(1+zLyaC),    color ='red', fmt='o', markersize=4,
              label="$\\rm{BOSS}\ \\mathrm{Ly}\\alpha\\mbox{-}\\rm{cross}\ \\rm{DR11}$", empty=empty2, ax=ax4)
plot_errorbar(zLyaC,  9.0*zLyaC,        yerr=0.3*zLyaC,        color ='green', fmt='o', markersize=4, empty=empty2, ax=ax3)


#Axis
ax4.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax4.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

#ax4.xaxis.set_minor_formatter(plt.ScalarFormatter())
#ax4.xaxis.set_major_locator(plt.FixedLocator([0.1,1.0]))
#ax4.xaxis.set_minor_locator(plt.FixedLocator([0.2,0.5,2]))

ax4.set_xlabel("$z$")

"""
if plaw>0:
    ax2.set_ylabel("${\\rm distance}/r_d\\sqrt{z}$")
    pylab.xlabel("$z$")
    #pylab.tight_layout()
    if (eboss):
        pylab.legend(loc='lower center', numpoints=1,bbox_to_anchor = (0.4, 0.0))
        pylab.ylim(3,30)
        pylab.xlim(0.0,2.6)
        pylab.savefig("Fig1_"+fname+"_DR12_eboss.pdf")
    else:
        pylab.legend(loc='upper left', numpoints=1, frameon=True)
        pylab.ylim(6,30)
        pylab.xlim(0.08,3.0)
        pylab.savefig("Fig1_"+fname+"_DR12.pdf")

else:
    pylab.legend(loc='lower left', numpoints=1)
    pylab.ylabel("$D(z)/r_d \log(1+z)$")
    if eboss:
        pylab.ylim(12,35)
        pylab.xlim(0.0,2.5)
    else:
        pylab.ylim(15,35)
        pylab.xlim(0.08,3.0)
    pylab.xlabel("$z$")
    pylab.tight_layout()
    if (eboss):
        pylab.savefig("Fig1_"+fname+"_DR12_v2_eboss.pdf")
    else:
        pylab.savefig("Fig1_"+fname+"_DR12_v2.pdf")
"""


#pylab.savefig("Fig1_"+fname+".pdf", bbox_inches='tight')
pylab.show()

