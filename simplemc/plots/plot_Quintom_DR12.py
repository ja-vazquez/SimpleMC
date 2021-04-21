#!/usr/bin/env python

#TODO Add Omega_de to the plots

from simplemc.plots.plot_Quintom_variables import *
from simplemc.models.QuintomCosmology import QuintomCosmology
from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker
import numpy as np
import pylab


steps = 9
coupling = 0
zl = np.arange(0, 3, 0.05)

fname = 'Quintessence'

# ---
if fname == 'Quintessence':
    T = QuintomCosmology(vary_mquin=True)
    name = fname
    mlabel = '$m_\phi$'

if fname == 'Phantom':
    T = QuintomCosmology(vary_mphan=True)
    name = fname
    mlabel = '$m_\psi$'

if fname == 'Quintom_mquin':
    T = QuintomCosmology(vary_mquin=True, vary_mphan=True)
    mphan = 1.2
    name = 'Quintom, $m_{\psi}$=%0.1f'%mphan
    mlabel = '$m_\phi$'

if fname == 'Quintom_mphan':
    T = QuintomCosmology(vary_mquin=False, vary_mphan=True)
    mphi = 1.2
    name = 'Quintom, $m_{\phi}$=%0.1f'%mphi
    mlabel = '$m_\psi$'

if fname == 'Quintom_coupling_mquin':
    T = QuintomCosmology(vary_mquin=True, vary_coupling=True)
    mphan = 1.2
    coupling = 4.0
    name = 'Quintom, $m_{\psi}$=%0.1f, $\\beta=%0.1f$'%(mphan, coupling)
    mlabel = '$m_\phi$'

if fname == 'Quintom_coupling_mphan':
    T = QuintomCosmology(vary_mphan=True, vary_coupling=True)
    mphi = 1.0 #1.2
    coupling = 10 #6.0
    name = 'Quintom, $m_{\phi}$=%0.1f, $\\beta=%0.1f$'%(mphi, coupling)
    mlabel = '$m_\psi$'

if fname == 'Quintom_coupling_both':
    T = QuintomCosmology(vary_mquin=True, vary_mphan=True, vary_coupling=True)
    mphi = 2.0
    mphan = 1.0
    coupling = -1
    name = 'Quintom, $m_{\phi}$=%0.1f, $m_{\psi}$=%0.1f'%(mphi, mphan)
    mlabel = '$\\beta$'



if fname == 'Quintom_coupling_both':
    min, max = (4., 8.)
else:
    min, max = (0.1, 2.5)
if coupling < 0:
    min, max = (-10, -1.)


step = (max-min)/steps

mquin_ = mquin_par
mphan_ = mphan_par
coupling_ = coupling_par

hh = []
ww = []
da = []
dh = []
zz = []
PP = []


for i in np.arange(min, max, step):
    if fname == 'Quintessence':
        mquin_.setValue(i)
        T.updateParams([mquin_])

    if fname == 'Phantom':
        mphan_.setValue(i)
        T.updateParams([mphan_])


    if fname == 'Quintom_mquin':
        mquin_.setValue(i)
        mphan_.setValue(mphan)
        T.updateParams([mquin_, mphan_])


    if fname == 'Quintom_mphan':
        mphan_.setValue(i)
        mquin_.setValue(mphi)
        T.updateParams([mquin_, mphan_])

    if fname == 'Quintom_coupling_mquin':
        mquin_.setValue(i)
        mphan_.setValue(mphan)
        coupling_.setValue(coupling)
        T.updateParams([mquin_, mphan_, coupling_])

    if fname == 'Quintom_coupling_mphan':
        mquin_.setValue(mphi)
        mphan_.setValue(i)
        coupling_.setValue(coupling)
        T.updateParams([mquin_, mphan_, coupling_])

    if fname == 'Quintom_coupling_both':
        mquin_.setValue(mphi)
        mphan_.setValue(mphan)
        coupling_.setValue(i)
        T.updateParams([mquin_, mphan_, coupling_])


    T.call_functions()
    ww.append([T.w_de(1./(1+z)) for z in zl])
    hh.append([T.Hubble_a(1./(1+z)) for z in zl])
    dh.append([T.HIOverrd(z)*z/fixer(z) for z in zl])
    da.append([T.DaOverrd(z)/fixer(z) for z in zl])
    PP.append(i)
    zz.append(zl)

#Planck best fit cosmology
T2 = LCDMCosmology()
x1 = [67.4*np.sqrt(T2.RHSquared_a(1./(1+z))) for z in zl]
x2 = [T2.HIOverrd(z)*z/fixer(z) for z in zl]
x3 = [T2.DaOverrd(z)/fixer(z) for z in zl]
#PLK-15
#T=LCDMCosmology(Obh2=0.02225,Om=0.3156,h=0.6727)




params1 = {'backend': 'pdf',
               'axes.labelsize': 18,
               'xtick.labelsize': 18,
               'ytick.labelsize': 18,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}
pylab.rcParams.update(params1)


## --- Plotting --
fig, (ax1, ax2, ax3, ax4)= plt.subplots(4, sharex=True, gridspec_kw={'hspace': 0}, figsize=(7,10))
fig.suptitle(name, fontsize=17,  y=0.95)

## -- Plot 1
for x, w, z in zip(zz, ww, PP):
    g = (float(z) - min)/(max - min)
    b, r = 0, 1 - g
    ax1.plot(x, w, color=(r, g, b))
if (fname == 'Quintessence') or (fname == 'Quintomcopphi'):
    ax1.set_ylabel('$w(z)$', fontsize=20)
ax1.axhline(y=-1.0, color='k', linestyle='--')
if coupling < 0:
    ax1.set_ylim(-3, 0.)


## -- Plot 2
for x, w, z in zip(zz, hh, PP):
    g = (float(z)-min)/(max-min)
    b, r = 0, 1-g
    (l2,) = ax2.plot(x, w, color=(r, g, b))

dataHz = np.loadtxt('simplemc/data/Hz_all.dat')
redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
ax2.errorbar(redshifts, obs, errors, xerr=None, color='blue', marker='o', ls='None',
                elinewidth =2, capsize=3, capthick = 1, alpha=1, markersize=4)
ax2.plot(zl, x1 , color='k', linestyle='--')
if (fname == 'Quintessence') or (fname =='Quintomcopphi'):
    ax2.set_ylabel('$H(z)$', fontsize=20)

Z = [[0,0] , [0,0]]
mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['red','green'])
levels = np.arange(min, max+step, step)
CS3 = plt.contourf(Z, levels, cmap=mymap)

cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.78])
cbar = pylab.colorbar(CS3, cax=cbaxes)
cbar.set_label(mlabel, rotation=0, fontsize=18, labelpad=-10)
cbar.ax.tick_params(labelsize=12)


## -- Plot 3
for x, w, z in zip(zz, dh, PP):
    g = (float(z)-min)/(max-min)
    b, r = 0, 1-g
    (l3,) = ax3.plot(x, w, color=(r ,g, b))
ax3.plot(zl, x2 , color='k', linestyle='--')
if (fname == 'Quintessence') or (fname == 'Quintomcopphi'):
    ax3.set_ylabel("${\\rm zD_H(z)}/r_d\\sqrt{z}$")


## -- Plot 4
for x, w, z in zip(zz, da, PP):
    g = (float(z)-min)/(max-min)
    b, r = 0, 1-g
    (l4,) = ax4.plot(x, w, color=(r, g, b))
ax4.plot(zl, x3, color='k', linestyle='--')
ax4.set_xlim(0.05, 3)
if (fname == 'Quintessence') or (fname == 'Quintomcopphi'):
    ax4.set_ylabel("${\\rm D_M(z)}/r_d\\sqrt{z}$")



plot_errorbar(zCombBAO1, 1512.4/rd_fid_DR12, yerr=ersys(22.5, 11.0)/rd_fid_DR12,
              color ='red', fmt='o', markersize=4, empty=empty2,
              label="$\\rm{BOSS\ Galaxy\ DR12}$", ax=ax4)
plot_errorbar(zCombBAO2, 1975.2/rd_fid_DR12, yerr=ersys(26.6, 14.1)/rd_fid_DR12,
              color ='red', fmt='o', markersize=4, empty=empty2, ax=ax4)
plot_errorbar(zCombBAO3, 2306.7/rd_fid_DR12, yerr=ersys(33.2, 16.7)/rd_fid_DR12,
              color ='red', fmt='o', markersize=4, empty=empty2, ax=ax4)


plot_errorbar(zCombBAO1,  fact*zCombBAO1/81.21, yerr=fact*zCombBAO1*ersys(2.17, 0.97)/(81.21)**2,
              color ='green', fmt='o', markersize=4, empty=empty2, ax=ax3)
plot_errorbar(zCombBAO2,  fact*zCombBAO2/90.90, yerr=fact*zCombBAO2*ersys(2.07, 1.08)/(90.90)**2,
              color ='green', fmt='o', markersize=4, empty=empty2, ax=ax3)
plot_errorbar(zCombBAO3,  fact*zCombBAO3/98.96, yerr=fact*zCombBAO3*ersys(2.21, 1.18)/(98.96)**2,
              color ='green', fmt='o', markersize=4, empty=empty2, ax=ax3)


plot_errorbar(zLyaA, 11.28*(1+zLyaA), yerr=0.65*(1+ zLyaA), color ='red', fmt='o',
              markersize=4, label="$\\rm{BOSS}\ \\mathrm{Ly}\\alpha\\mbox{-}\\rm{auto}\ \\rm{DR11}$",
              empty=empty2, ax=ax4)
plot_errorbar(zLyaA, 9.18*zLyaA, yerr=0.28*zLyaA, color ='green', fmt='o',
              markersize=4,empty=empty2, ax=ax3)
plot_errorbar(zLyaC, 10.8*(1+zLyaC), yerr=0.4*(1+zLyaC), color ='red', fmt='o',
              markersize=4, label="$\\rm{BOSS}\ \\mathrm{Ly}\\alpha\\mbox{-}\\rm{cross}\ \\rm{DR11}$",
              empty=empty2, ax=ax4)
plot_errorbar(zLyaC, 9.0*zLyaC, yerr=0.3*zLyaC, color ='green', fmt='o',
              markersize=4, empty=empty2, ax=ax3)


#Axis
ax4.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax4.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#ax4.xaxis.set_minor_formatter(plt.ScalarFormatter())
#ax4.xaxis.set_major_locator(plt.FixedLocator([0.1,1.0]))
#ax4.xaxis.set_minor_locator(plt.FixedLocator([0.2,0.5,2]))
ax4.set_xlabel("$z$")


#pylab.savefig("Fig1_"+fname+".pdf", bbox_inches='tight')
pylab.show()

