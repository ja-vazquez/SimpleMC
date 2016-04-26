#!/usr/bin/env python

#Error bars mainly got them from 
#table 3 at
#http://arxiv.org/pdf/1108.2635v1.pdf

from RunBase import *
import matplotlib.pyplot as plt
import matplotlib.ticker 
import pylab 
import math as N
plaw=0.5


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

zLOWZ  = 0.32 
zCMASS = 0.57
zLyaA  = 2.34
zLyaC  = 2.36

z6dFGS   = 0.106
zMGS     = 0.15
zSDSS1   = 0.2
zSDSS2   = 0.35
zWiggleZ1=0.44
zWiggleZ2= 0.6
zWiggleZ3= 0.73
z_CMB = 1090.43

rd_EHtoCAMB =153.19/149.28

zl=arange(0,8,0.03)
zle=np.exp(zl)-1


y1=[T.DaOverrd(z)/z**plaw   for z in zle]
y1W=[T2.DaOverrd(z)/z**plaw for z in zle]

y2=[T.HIOverrd(z)*z/z**plaw for z in zle]
y2W=[T2.HIOverrd(z)*z/z**plaw for z in zle]

y3=[T.DVOverrd(z)/z**plaw   for z in zle]
y3W=[T2.DVOverrd(z)/z**plaw   for z in zle]

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(1,1,1)

l1,=plt.plot(zle,y1,'r-',lw=3) 
pylab.plot(zle,y1W,'r--',lw=2)

l2,=plt.plot(zle,y3,'b-',lw=3) 
pylab.plot(zle,y3W,'b--',lw=2)

l3,=plt.plot(zle,y2,'g-',lw=3) 
pylab.plot(zle,y2W,'g--',lw=2)

legend1=pylab.legend([l1,l2,l3],["$D_%s(z)/r_d\\sqrt{z}$"%st for st in ['M','V','H']],loc="upper left")
pylab.gca().add_artist(legend1)


def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
                txt.set_color(line.get_color())

legend1.draw_frame(False)
color_legend(legend1)


### Plotting -  Error bars 
def plot_errorbar(z,val, yerr=0, color=0, fmt=0, markersize=0,label=None):
    pylab.errorbar(z,val/z**plaw, yerr=yerr/z**plaw, color=color, lw=2, capsize=4, fmt=fmt, markersize=markersize)
    if label>0:
        pylab.plot ([],[],fmt,color='black',label=label,markersize=markersize)

pylab.legend(loc="lower right")

plot_errorbar(z6dFGS, 2.97*rd_EHtoCAMB, yerr=rd_EHtoCAMB*0.015/0.336**2,  color ='green', fmt='+', markersize=8, label="$\\rm{6dFGS}$")
plot_errorbar(zMGS,   666/148.651,      yerr=25/148.651,       color ='blue', fmt='^', markersize=8, label="$\\rm{MGS}$")
plot_errorbar(zLOWZ,  8.467,            yerr=0.17,             color ='blue', fmt='v', markersize=8, label="$\\mathrm{LOWZ}$")
plot_errorbar(zCMASS, 9.519*(1+zCMASS), yerr=0.134*(1+zCMASS), color ='red', fmt='d', markersize=8,label="$\\mathrm{CMASS}$")
plot_errorbar(zLyaA,  11.28*(1+zLyaA),  yerr=0.65*(1+ zLyaA),  color ='red', fmt='o', markersize=8,label="$\\mathrm{Ly}\\alpha\ \\rm{auto}$")
plot_errorbar(zLyaC,  10.8*(1+zLyaC),   yerr=0.4*(1+zLyaC),    color ='red', fmt='*', markersize=8,label="$\\mathrm{Ly}\\alpha\ \\rm{cross}$")
plot_errorbar(z_CMB,  94.22,            yerr=0.6,              color ='red', fmt='p', markersize=8,label="$\\mathrm{CMB}\ D_M/r_D$")

plot_errorbar(zCMASS,20.75*zCMASS, yerr=0.73*zCMASS,color ='green', fmt='d', markersize=8)
plot_errorbar(zLyaA, 9.18*zLyaA,   yerr=0.28*zLyaA, color ='green', fmt='o', markersize=8)
plot_errorbar(zLyaC, 9.0*zLyaC,    yerr=0.3*zLyaC,  color ='green', fmt='*', markersize=8)


legend2 = pylab.legend(loc='upper right', numpoints=1)
legend2.draw_frame(False)

ax.set_xscale('log')
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
ax.set_yticks([0, 10 ,20 ,30])

pylab.xlim(0.05,1200)
pylab.ylim(0,31)
pylab.xlabel("$z$")
pylab.ylabel("$D(z)/r_d\\sqrt{z}$")
pylab.tight_layout()
pylab.savefig("FigW.pdf")
pylab.show()
