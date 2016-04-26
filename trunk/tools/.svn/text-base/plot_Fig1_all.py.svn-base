#!/usr/bin/env python

#Error bars mainly got from  table 3,
#http://arxiv.org/pdf/1108.2635v1.pdf
#and MGS paper

from RunBase import *
import matplotlib.pyplot as plt
import matplotlib.ticker
import pylab 
import math as N

plaw=0.5
#for division by log(1+z) use this
#plaw=-1

eboss=False 
#True


params1 = {'backend': 'pdf',
               'axes.labelsize': 20,
               'text.fontsize': 18,
               'xtick.labelsize': 20,
               'ytick.labelsize': 20,
               'legend.draw_frame': False,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}#
pylab.rcParams.update(params1)


#Planck best fit cosmology
T=LCDMCosmology(Obh2=0.022032,Om=0.3183,h=0.6704)

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

rd_EHtoCAMB =153.19/149.28

zl=arange(0,4,0.01)

def fixer(z):
    if plaw>0:
        return z**plaw
    else:
        return log(1.+z)


y1=[T.DaOverrd(z)/fixer(z)   for z in zl]
y2=[T.HIOverrd(z)*z/fixer(z) for z in zl]
y3=[T.DVOverrd(z)/fixer(z)   for z in zl]

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(1,1,1)
if (not eboss):
    ax.set_xscale('log')


l1,=plt.plot(zl,y1,'r-',lw=2)
l2,=plt.plot(zl,y3,'b-',lw=2)
l3,=plt.plot(zl,y2,'g-',lw=2)

if plaw>0:
    if eboss:
        legend1=pylab.legend([l1,l2,l3],["$D_%s(z)/r_d\\sqrt{z}$"%st for st in ['M','V','H']],loc="lower center", bbox_to_anchor = (0.7, 0.0))
    else:
        legend1=pylab.legend([l1,l2,l3],["$%sD_%s(z)/r_d\\sqrt{z}$"%st for st in [('','M'),('','V'),('z','H')]],loc="lower center")
else:
    legend1=pylab.legend([l1,l2,l3],["$D_%s(z)/r_d\log(1+z)$"%st for st in ['A','V','H']],loc="lower center")


pylab.gca().add_artist(legend1)
def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
                txt.set_color(line.get_color())

legend1.draw_frame(False)
color_legend(legend1)


### Plotting -  Error bars 
def plot_errorbar(z,val, yerr=0, color=0, fmt=0, markersize=0,label=None, empty=True):
    if empty:
        mfc='white'
        lw=1
    else:
        mfc=color
        lw=2
    pylab.errorbar(z,val/fixer(z), yerr=yerr/fixer(z), color=color, fmt=fmt, markersize=markersize, lw=lw, capthick=lw,capsize=2+2*lw,markerfacecolor=mfc)
    if label>0:
        if (mfc=='white'):
            pylab.plot ([],[],fmt,color='black',label=label,markersize=markersize,markerfacecolor=mfc)
        else:
            pylab.plot ([],[],fmt,color='black',label=label,markersize=markersize)

pylab.legend(loc="lower right")


plot_errorbar(z6dFGS,    2.97*rd_EHtoCAMB,   yerr=rd_EHtoCAMB*0.015/0.336**2,  color ='blue', fmt='+', markersize=8, label="$\\rm{6dFGS}$")
plot_errorbar(zMGS,      666/148.651,    yerr=25/148.651,               color ='blue', fmt='^', markersize=8, label="$\\rm{MGS}$")
plot_errorbar(zSDSS1,    5.2493*rd_EHtoCAMB, yerr=rd_EHtoCAMB*0.0061/0.1905**2,color ='blue', fmt='p', markersize=8, label="$\\rm{SDSS-II}$")
plot_errorbar(zSDSS2,    9.1157*rd_EHtoCAMB, yerr=rd_EHtoCAMB*0.0036/0.1097**2,color ='blue', fmt='p', markersize=8)
plot_errorbar(zWiggleZ1, 10.917*rd_EHtoCAMB ,yerr=rd_EHtoCAMB*0.0071/0.0916**2,color ='blue', fmt='s', markersize=8, label="$\\rm{WiggleZ}$")
plot_errorbar(zWiggleZ2, 13.774*rd_EHtoCAMB ,yerr=rd_EHtoCAMB*0.0034/0.0726**2,color ='blue', fmt='s', markersize=8)
plot_errorbar(zWiggleZ3, 16.891*rd_EHtoCAMB ,yerr=rd_EHtoCAMB*0.0031/0.0592**2,color ='blue', fmt='s', markersize=8)

if eboss:
    plot_errorbar(zLOWZ,     8.467,              yerr=0.17,        color ='blue', fmt='v', markersize=8, label="$\\mathrm{LOWZ}$",empty=False)
else:
    plot_errorbar(zLOWZ,     8.467,              yerr=0.17,        color ='blue', fmt='v', markersize=8, label="$\\mathrm{LOWZ}$",empty=False)


plot_errorbar(zCMASS, 9.519*(1+zCMASS), yerr=0.134*(1+zCMASS), color ='red', fmt='d', markersize=6,label="$\\mathrm{CMASS}$", empty=False)
plot_errorbar(zCMASS, 20.75*zCMASS,     yerr=0.73*zCMASS,      color ='green', fmt='-d', markersize=6,empty=False)

plot_errorbar(zLyaA,  11.28*(1+zLyaA),  yerr=0.65*(1+ zLyaA),  color ='red', fmt='o', markersize=8,label="$\\mathrm{Ly}\\alpha\ \\rm{auto}$", empty=False)
plot_errorbar(zLyaA,  9.18*zLyaA,       yerr=0.28*zLyaA,       color ='green', fmt='-o', markersize=8,empty=False)

plot_errorbar(zLyaC,  10.8*(1+zLyaC),   yerr=0.4*(1+zLyaC),    color ='red', fmt='*', markersize=8,label="$\\mathrm{Ly}\\alpha\ \\rm{cross}$",empty=False)
plot_errorbar(zLyaC,  9.0*zLyaC,        yerr=0.3*zLyaC,        color ='green', fmt='-*', markersize=8, empty=False)


#Axis
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.xaxis.set_minor_formatter(plt.ScalarFormatter())
ax.xaxis.set_major_locator(plt.FixedLocator([0.1,1.0]))
ax.xaxis.set_minor_locator(plt.FixedLocator([0.2,0.5,2]))

plt.yticks(range(0, 50, 10))



if plaw>0:
    pylab.ylabel("${\\rm distance}/r_d\\sqrt{z}$")
    pylab.xlabel("$z$")
    pylab.tight_layout()
    if (eboss):
        pylab.legend(loc='lower center', numpoints=1,bbox_to_anchor = (0.4, 0.0))
        pylab.ylim(3,30)
        pylab.xlim(0.0,2.6)
        pylab.savefig("Fig1_all_eboss.pdf")
    else:
        pylab.legend(loc='upper left', numpoints=1)
        pylab.ylim(6,30)
        pylab.xlim(0.08,3.0)
        pylab.savefig("Fig1_all.pdf")

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
        pylab.savefig("Fig1_all_v2_eboss.pdf")
    else:
        pylab.savefig("Fig1_all_v2.pdf")


pylab.show()

