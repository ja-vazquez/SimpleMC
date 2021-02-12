
import numpy as np

addlabel = False
addlegend = False
eboss = False
plaw = 0.5


zLOWZ = 0.32
zCMASS = 0.57
zLyaA = 2.34 - 0.04
zLyaC = 2.36 + 0.04

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

rd_EHtoCAMB = 153.19/149.28
rd_fid_DR12 = 147.78
rd_fid_DR7 = 151.84

fact = (300000./rd_fid_DR12)

fmt1 = 'o'
fmt2 = 'o'
empty1 = False
empty2 = False
alpha = 1.0


#============Functions -------

def fixer(z):
    if plaw > 0:
        return z**plaw
    else:
        return np.log(1.+z)


# Plotting -  Error bars
def plot_errorbar(z,val, yerr=0, color=0, fmt=0, markersize=0,label=None, empty=True, alpha=1, ax=None):
    if empty:
        mfc='white'
        lw=1
    else:
        mfc=color
        lw=1
    color = 'black'

    ax.errorbar(z, val/fixer(z), yerr=yerr/fixer(z),
                color='blue', marker='o', ls='None',
                elinewidth =2, capsize=3, capthick = 1, alpha=1, markersize=4)

    if addlabel:
        if label>0:
            if mfc == 'white':
                ax.plot([], [], fmt, color='black', label=label,
                        markersize=markersize, markerfacecolor=mfc)
            else:
                ax.plot([], [], fmt, color='black', label=label,
                        markersize=markersize)


def ersys(x, y):
    return np.sqrt(x**2 + y**2)


def color_legend(leg):
        """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
                txt.set_color(line.get_color())

#-------------------



