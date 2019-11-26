#!/usr/bin/env python

# Error bars mainly got from  table 3,
# http://arxiv.org/pdf/1108.2635v1.pdf
# and MGS paper

from RunBase import *
import matplotlib.pyplot as plt
import matplotlib.ticker
import pylab
import math as N


plaw = 0.5
# for division by log(1+z) use this
# plaw=-1


params1 = {'backend': 'pdf',
           'axes.labelsize': 20,
           'text.fontsize': 18,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           #               'legend.draw_frame': False,
           'legend.fontsize': 16,
           'lines.markersize': 6,
           'font.size': 20,
           'text.usetex': True}
pylab.rcParams.update(params1)


# Planck best fit cosmology
T = LCDMCosmology(Obh2=0.022, Om=0.31, h=0.676)
# PLK-15
# T=LCDMCosmology(Obh2=0.02225,Om=0.3156,h=0.6727)

zLOWZ = 0.32
zCMASS = 0.57
zLyaA = 2.33
zLyaC = 2.40

z6dFGS = 0.106
zMGS = 0.15
zSDSS1 = 0.2
zSDSS2 = 0.35
zWiggleZ1 = 0.44
zWiggleZ2 = 0.6
zWiggleZ3 = 0.73

z_CMB = 1090.43

zCombBAO1 = 0.38
zCombBAO2 = 0.51
zCombBAO3 = 0.61

zEBQSO = 1.52


rd_EHtoCAMB = 153.19/149.28
rd_fid_DR12 = 147.78
rd_fid_DR7 = 151.84

zl = arange(0, 4, 0.01)


def fixer(z):
    if plaw > 0:
        return z**plaw
    else:
        return log(1.+z)


y1 = [T.DaOverrd(z)/fixer(z) for z in zl]
y2 = [T.HIOverrd(z)*z/fixer(z) for z in zl]
y3 = [T.DVOverrd(z)/fixer(z) for z in zl]

fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(1, 1, 1)
ax.set_xscale('log')

l1, = plt.plot(zl, y1, 'r-', lw=2)
l2, = plt.plot(zl, y3, 'b-', lw=2)
l3, = plt.plot(zl, y2, 'g-', lw=2)

if plaw > 0:
    legend1 = pylab.legend([l1, l2, l3], ["$%sD_%s(z)/r_d\\sqrt{z}$" % st for st in [
                           ('', 'M'), ('', 'V'), ('z', 'H')]], loc="lower center")
else:
    legend1 = pylab.legend([l1, l2, l3], ["$D_%s(z)/r_d\log(1+z)$" %
                                          st for st in ['A', 'V', 'H']], loc="lower center")


pylab.gca().add_artist(legend1)


def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
    for line, txt in zip(leg.get_lines(), leg.get_texts()):
        txt.set_color(line.get_color())


legend1.draw_frame(False)
color_legend(legend1)


def ersys(x, y):
    return sqrt(x**2 + y**2)


# Plotting -  Error bars
def plot_errorbar(z, val, yerr=0, color=0, fmt=0, markersize=0, label=None, empty=True, alpha=1):
    if empty:
        mfc = 'white'
        lw = 1
    else:
        mfc = color
        lw = 2
    pylab.errorbar(z, val/fixer(z), yerr=yerr/fixer(z), color=color, fmt=fmt, markersize=markersize,
                   lw=lw, capthick=lw, capsize=2+2*lw, markerfacecolor=mfc, alpha=alpha)
    if label > 0:
        if (mfc == 'white'):
            pylab.plot([], [], fmt, color='black', label=label,
                       markersize=markersize, markerfacecolor=mfc)
        else:
            pylab.plot([], [], fmt, color='black',
                       label=label, markersize=markersize)


pylab.legend(loc="lower right")

fmt1 = '^'
fmt2 = 's'
empty1 = True
empty2 = False
alpha = 1.0

# Errorbars from DR12 Full-shape
fact = (300000./rd_fid_DR12)

# 666/148.651,    yerr=25/148.651
plot_errorbar(z6dFGS,    2.97*rd_EHtoCAMB,   yerr=rd_EHtoCAMB*0.015/0.336**2,
              color='blue', fmt='>', markersize=6, empty=empty2, label="$\\rm{6dFGS}$", alpha=alpha)
plot_errorbar(zMGS,      4.464,    yerr=0.168,               color='blue',
              fmt='p', markersize=6, label="$\\rm{SDSS\ MGS}$", empty=empty2, alpha=alpha)

plot_errorbar(zSDSS1,    5.2493*rd_EHtoCAMB, yerr=rd_EHtoCAMB*0.0061/0.1905 **
              2, color='blue', fmt=fmt1, markersize=6, empty=empty1, alpha=alpha)
plot_errorbar(zSDSS2,    1348./rd_fid_DR12, yerr=26./rd_fid_DR12, color='blue',
              fmt=fmt1, markersize=6, label="$\\rm{SDSS\ DR7}$", empty=empty1, alpha=alpha)

plot_errorbar(zWiggleZ1, 1695./rd_fid_DR12, yerr=82./rd_fid_DR12, color='blue',
              fmt=fmt2, markersize=6, label="$\\rm{WiggleZ}$", empty=empty1, alpha=alpha)
plot_errorbar(zWiggleZ2, 2194./rd_fid_DR12, yerr=100./rd_fid_DR12,
              color='blue', fmt=fmt2, markersize=6, empty=empty1, alpha=alpha)
plot_errorbar(zWiggleZ3, 2486./rd_fid_DR12, yerr=85./rd_fid_DR12,
              color='blue', fmt=fmt2, markersize=6, empty=empty1, alpha=alpha)


plot_errorbar(zCombBAO1,  1512.4/rd_fid_DR12,     yerr=ersys(22.5, 11.0)/rd_fid_DR12,
              color='red', fmt='d', markersize=8, empty=empty2, label="$\\rm{BOSS\ Galaxy\ DR12}$")
plot_errorbar(zCombBAO2,  1975.2/rd_fid_DR12,     yerr=ersys(26.6, 14.1)/rd_fid_DR12,
              color='red', fmt='d', markersize=8, empty=empty2)
plot_errorbar(zCombBAO3,  2306.7/rd_fid_DR12,  	 yerr=ersys(33.2, 16.7)/rd_fid_DR12,
              color='red', fmt='d', markersize=8, empty=empty2)


plot_errorbar(zCombBAO1,  fact*zCombBAO1/81.21,       yerr=fact*zCombBAO1*ersys(2.17, 0.97)/(81.21)**2,
              color='green', fmt='d', markersize=8, empty=empty2)
plot_errorbar(zCombBAO2,  fact*zCombBAO2/90.90,       yerr=fact*zCombBAO2*ersys(2.07, 1.08)/(90.90)**2,
              color='green', fmt='d', markersize=8, empty=empty2)
plot_errorbar(zCombBAO3,  fact*zCombBAO3/98.96,       yerr=fact*zCombBAO3*ersys(2.21, 1.18)/(98.96)**2,
              color='green', fmt='d', markersize=8, empty=empty2)


plot_errorbar(zEBQSO,  3855/rd_fid_DR12, yerr=170/rd_fid_DR12,  color='blue',
              fmt='v', markersize=8, label="$\\rm{eBOSS\ QSO}$", empty=empty2, alpha=alpha)

plot_errorbar(zLyaA,  37.77,  yerr=2.13,  color='red', fmt='o', markersize=8,
              label="$\\rm{BOSS}\ \\mathrm{Ly}\\alpha\\mbox{-}\\rm{auto}\ \\rm{DR12}$", empty=empty2)
plot_errorbar(zLyaA,  9.07*zLyaA,       yerr=0.31*zLyaA,
              color='green', fmt='-o', markersize=8, empty=empty2)

plot_errorbar(zLyaC,  35.7,   yerr=1.5,    color='red', fmt='*', markersize=8,
              label="$\\rm{BOSS}\ \\mathrm{Ly}\\alpha\\mbox{-}\\rm{cross}\ \\rm{DR12}$", empty=empty2)
plot_errorbar(zLyaC,  9.01*zLyaC,        yerr=0.32*zLyaC,
              color='green', fmt='-*', markersize=8, empty=empty2)


# Axis
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.xaxis.set_minor_formatter(plt.ScalarFormatter())
ax.xaxis.set_major_locator(plt.FixedLocator([0.1, 1.0]))
ax.xaxis.set_minor_locator(plt.FixedLocator([0.2, 0.5, 2]))

plt.yticks(list(range(0, 50, 10)))


if plaw > 0:
    pylab.ylabel("${\\rm distance}/r_d\\sqrt{z}$")
    pylab.xlabel("$z$")
    pylab.tight_layout()
    pylab.legend(loc='upper left', numpoints=1, frameon=False)
    pylab.ylim(6, 32)
    pylab.xlim(0.08, 3.0)
    pylab.savefig("Fig1_DR14.pdf")

else:
    pylab.legend(loc='lower left', numpoints=1)
    pylab.ylabel("$D(z)/r_d \log(1+z)$")
    pylab.ylim(15, 35)
    pylab.xlim(0.08, 3.0)
    pylab.xlabel("$z$")
    pylab.tight_layout()
    pylab.savefig("Fig1_DR14_v2.pdf")


# pylab.show()
