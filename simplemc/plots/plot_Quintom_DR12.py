#!/usr/bin/env python

from simplemc.plots.plot_Quintom_variables import *
from simplemc.models.QuintomCosmology import QuintomCosmology
from simplemc.cosmo.paramDefs import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker
import numpy as np


beta  = 0
zl=np.arange(0, 3, 0.05)

fname = 'Quintessence'

if fname == 'Quintessence':
    T = QuintomCosmology(vary_mquin=True)
if fname == 'Phantom':
    T = QuintomCosmology(vary_mphan=True)
if fname == 'Quintom_mphi':
    T = QuintomCosmology(vary_mquin=True, vary_mphan=True)
    mphan = 1.2
if fname == 'Quintmphan':
    T = QuintomCosmology(varymquin=False, varymphan=True)
    mphi  = 1.2
if fname == 'Quintomcopphi':
    T = QuintomCosmology(varymquin=True, varybeta=True)
    mphan = 1.2
    beta  = 4.0
if fname == 'Quintomcopphan':
    T = QuintomCosmology(varymphan=True, varybeta=True)
    mphi  = 1.2
    beta  = 6.0
if fname == 'Quintomcopbeta':
    T = QuintomCosmology(varymquin=True, varymphan=True, varybeta=True)
    mphi  = 1.2
    mphan = 0.1


hh = []
ww = []
da = []
dh = []
zz = []
PP = []

if fname == 'Quintomcopbeta':
    min, max = (4., 8.)
else:
    min, max = (0.01, 1.5)
if beta < 0:  min, max = (-4, -1.)

steps    = 5
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

    T.call_functions()
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
ax1.axhline(y=-1.05, color='k', linestyle='--')
if beta <0 : ax1.set_ylim(-3, 0.)

#ax1.set_xscale('log')
for x,w,z in zip(zz, hh, PP):
        g    = (float(z)-min)/(max-min)
        b, r = 0, 1-g
        l2, = ax2.plot(x, w , color=(r,g,b)) #, linestyle='-.')
#ax2.plot(zl, x1 , color='k', linestyle='--')
dataHz = np.loadtxt('simplemc/data/Hz_all.dat')
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

