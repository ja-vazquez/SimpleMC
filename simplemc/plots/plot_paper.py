#!/usr/bin/env python

from scipy import *
from cosmich import *
from plot_aux import *
import pylab
import sys

params1 = {'backend': 'pdf',
           'axes.labelsize': 32,
           'xtick.labelsize': 24,
           'ytick.labelsize': 24,
           'legend.fontsize': 20,
           'lines.markersize': 6,
           'font.size': 20,
           'text.usetex': True}
pylab.rcParams.update(params1)

# -----------

dire = 'chains_SimpleMC/'
# for decay!
dire = 'chains/'
# 3dire="decay2/"
#dire = 'chains_140630_157330/'
show = True  # False


if len(sys.argv) > 1:
    lst = sys.argv[1:]
else:
    lst = ['S5', 'S3', 'Early', 'HLadder', 'smoking', 'SlowRDE']

T = plot_aux()

for l in lst:
    if l == 'S5':

        datasetl = ['BBAO+Planck', 'SN+Planck', 'BBAO+SN+Planck']
        extra = 'phy'

        model = ['LCDM']
        params = ['h', 'Om']
        T.Plotting_2d(dire, model, extra, datasetl,
                      params, legcolor=True, Nb=120)
        pylab.title('$\\Lambda$CDM')
        pylab.xlim(0.64, 0.74)
        pylab.ylim(0.23, 0.37)
        pylab.tight_layout()
        pylab.savefig('deplot20.pdf')
        if show:
            pylab.show()

        model = ['oLCDM']
        params = ['Om', 'Ok']
        T.Plotting_2d(dire, model, extra, datasetl,
                      params, Nb=40, legcolor=True)
        pylab.title('$o\\Lambda$CDM')
        pylab.xlim(0.18, 0.45)
        pylab.ylim(-0.05, 0.05)
        pylab.tight_layout()
        pylab.savefig('deplot21.pdf')
        if show:
            pylab.show()

        model = ['wCDM']
        params = ['Om', 'w']
        T.Plotting_2d(dire, model, extra, datasetl, params,
                      Nb=40, loc="upper left", legcolor=True)
        pylab.title('$w$CDM')
        xll, xlh = pylab.xlim()
        pylab.plot([xll, xlh], [-1, -1], 'k:')
        pylab.xlim(0.21, 0.38)
        pylab.ylim(-1.4, -0.6)
        pylab.tight_layout()
        pylab.savefig('deplot22.pdf')
        if show:
            pylab.show()

        model = ['owCDM']
        params = ['w', 'Ok']
        T.Plotting_2d(dire, model, extra, datasetl, params,
                      Nb=40, loc="upper left", legcolor=True)
        pylab.title('$ow$CDM')
        pylab.axis([-1.6, -0.5, -0.05, 0.05])
        pylab.plot([-1, -1], [-0.05, 0.05], 'k:')
        #pylab.plot ([-1.6,-0.5],[0,0],'k:')
        pylab.tight_layout()
        pylab.savefig('deplot23.pdf')
        if show:
            pylab.show()

    if l == 'S5' or l == 'S5g':

        datasetl = ['BBAO+Planck', 'SN+Planck', 'BBAO+SN+Planck']
        extra = 'phy'

        model = ['waCDM']
        params = ['w_z', 'wa']
        T.Plotting_waDE_2d(dire, model, extra, datasetl,
                           params, legcolor=True, Nb=20)
        pylab.title('$w_ow_a$CDM')
        pylab.xlim(-1.5, -0.3)
        pylab.xlabel('$w_{(z=0.266)}$')
        pylab.tight_layout()
        pylab.savefig('deplot24.pdf')
        if show:
            pylab.show()

    if l == 'S5' or l == 'S5h':

        datasetl = ['BBAO+Planck', 'SN+Planck', 'BBAO+SN+Planck']
        extra = 'phy'

        model = ['owaCDM']
        params = ['w_z', 'Ok']
        T.Plotting_waDE_2d(dire, model, extra, datasetl,
                           params, loc='upper left', legcolor=True, Nb=25)
        pylab.title('$ow_ow_a$CDM')
        pylab.xlabel('$w_{(z=0.266)}$')
        pylab.ylim(-0.05, 0.05)
        pylab.xlim(-2, -0.3)
        pylab.tight_layout()
        pylab.savefig('deplot25.pdf')
        if show:
            pylab.show()

#######################
    if l == 'Decay2':
        datasetl = ['BBAO+SN+Planck']
        extra = 'phy'
        model = ['Decay01', 'Decay05', 'Decay']
        params = ['lambdaf', 'Om']
        T.Plotting_2d(dire, model, extra, datasetl, params, loc='upper right',
                      legcolor=True, Nbin=[25, 30, 30], labels=['$f_x=0.1$', '$f_x=0.5$', '$f_x=1$'])
        # pylab.xlabel('$f_d$')
        # pylab.xlabel('$\\lambda$')
        pylab.tight_layout()
        pylab.xlim(0, 1.2)
        pylab.savefig('decay2d.eps')
        if show:
            pylab.show()


# Extra

    if l == 'smoking' or l == 'neff':

        extra = 'phy'

        datasetl = ['NFC:r_PL_LBAO', 'NFC:r_PL_BBAO']
        model = ['Neff']
        params = ['nnu', 'h']
        fcol = ['red', 'green']
        lw = [3, 3]
        T.Plotting_2d(dire, model, extra, datasetl, params,
                      lwidth=lw, fillcolor=fcol, loc='upper left')
        pylab.ylim(0.63, 0.9)
        pylab.xlabel('$N_{\\rm eff}$')
        pylab.tight_layout()
        pylab.savefig('neff-h.pdf')
        if show:
            pylab.show()
        if False:
            datasetl = ['LBAO+Planck', 'BBAO+Planck', 'BBAO+SN+Planck',
                        'NFC:LBAO_Nnu_r', 'NFC:BBAO_Nnu_r']
            model = ['NeffLCDM']
            params = ['Nnu', 'h']
            T.Plotting_2d(dire, model, extra, datasetl,
                          params, loc='upper left')
            pylab.tight_layout()
            pylab.savefig('neff-h.pdf')
            if show:
                pylab.show()

    if l == 'smoking' or l == 'mnu':
        extra = 'phy'
        datasetl = ['BBAO+Planck', 'mnu_AL_PL_BBAO', 'mnu_PL_BBAO']
        model = ['nuLCDM']
        params = ['mnu']
        T.Plotting_1d(dire, model, extra, datasetl,
                      params, Nb=80, minmax=[0.0, 1.0])
        lims = [0.56, 0.44, 0.25]
        y1, y2 = pylab.ylim()
        for l, c, i in zip(lims, ['r', 'b', 'g'], list(range(3))):
            pylab.plot([l, l], [y1, y2], c+'--', lw=2)
            pylab.text(l, 3+i, '$<$%geV' % l, color=c)

        pylab.xlim(0, 1.0)
        pylab.tight_layout()
        pylab.savefig('mnu.pdf')
        if show:
            pylab.show()

    if l == 'smoking' or l == 'mnu2':
        extra = 'phy'
        datasetl = ['BBAO+Planck', 'BBAO+SN+Planck']
        model = ['nuLCDM', 'nuwCDM', 'nuoLCDM']
        params = ['mnu']
        leg = ['$\\Lambda$CDM+$m_\\nu$', '$\\Lambda$CDM+$m_\\nu$ (w SN)',
               'wCDM+$m_\\nu$', 'wCDM+$m_\\nu$ (w SN)',
               'o$\\Lambda$CDM+$m_\\nu$', 'o$\\Lambda$CDM+$m_\\nu$ (w SN)']
        lcolor = ['red', 'red', 'blue', 'blue', 'green', 'green']
        lw = [4, 2, 4, 2, 4, 2]
        T.Plotting_1d(dire, model, extra, datasetl, params, Nb=40,
                      linecolor=lcolor, lwidth=lw, legend=leg, minmax=[0.0, 1.0])
        pylab.xlim(0, 1.0)
        pylab.xlabel('$\\sum m_\\nu$')
        pylab.tight_layout()
        pylab.savefig('mnu2.pdf')
        if show:
            pylab.show()

    if l == 'smoking' or l == 'decay':
        extra = 'phy'

        datasetl = ['BBAO+SN+Planck']
        model = ['Decay', 'Decay05', 'Decay01']
        params = ['lambdaf']
        T.Plotting_1d(dire, model, extra, datasetl, params, Nb=[60, 60, 58]*2,
                      legend=['$f_x=1$', '$f_x=0.5$', '$f_x=0.1$']*2, minmax=[0.0, 1.0])
        pylab.tight_layout()
        pylab.xlim(0., 0.8)
        pylab.ylim(0, pylab.ylim()[1])
        pylab.savefig('decay.eps')
        if show:
            pylab.show()

    if l == 'smoking' or l == 'tlight':
        extra = 'phy'
        datasetl = ['BBAO+SN+PlRd']
        model = ['TLight']
        params = ['beta']
        T.Plotting_1d(dire, model, extra, datasetl, params, Nb=70)
        pylab.xlim(-0.2, 0.4)
        pylab.tight_layout()
        pylab.savefig('tlight.pdf')
        if show:
            pylab.show()

    if l == 'Early':

        extra = 'phy'

        datasetl = ['BBAO+SN+Planck']
        model = ['EarlyDE', 'EarlyDE_rd_DE']
        params = ['Ode', 'Om']
        T.Plotting_EDE_2d(dire, model, extra, datasetl, params, Nb=40)
        pylab.xlabel('$\\Omega^e_{\\rm de}$')
        pylab.xlim(0, 0.6)
        lgd = pylab.legend()
        lgd.set_visible(False)
        pylab.tight_layout()
        pylab.savefig('EarlyDE_rd.pdf')
        if show:
            pylab.show()

    elif l == 'SlowRDE':

        extra = 'phy'

        datasetl = ['BBAO+Planck', 'SN+Planck', 'BBAO+SN+Planck']
        model = ['SlowRDE']
        params = ['dw', 'h']
        T.Plotting_2d(dire, model, extra, datasetl,
                      params, Nb=35, legcolor=True)
        pylab.xlim(-0.5, 0.9)
        pylab.ylim(0.55, 0.8)
        pylab.tight_layout()
        pylab.savefig('SlowRDE.pdf')
        if show:
            pylab.show()


# Second set of Plots

    if l == 'S3':

        datasetl = ['GBAO', 'LBAO', 'BBAO', 'PLA:base:planck_lowl_lowLike:']
        model = ['LCDM']
        params = ['h', 'Om']
        extra = 'phy'
        T.Plotting_2d(dire, model, extra, datasetl, params,
                      loc='upper left', legcolor=True)
        pylab.axis([0.5, 1.0, 0.0, 1.0])
        pylab.tight_layout()
        pylab.savefig('lcdm.pdf')
        if show:
            pylab.show()

        extra = 'pre'
        datasetl = ['GBAO+PlDa', 'LBAO+PlDa', 'BBAO', 'BBAO+PlDa']
        model = ['oLCDM']
        params = ['Ol']
        T.Plotting_1d(dire, model, extra, datasetl, params,
                      Ol='True', loc='upper left', legcolor=True)
        pylab.xlim(0.0, 1.2)
        pylab.tight_layout()
        pylab.savefig('pre0.pdf')
        if show:
            pylab.show()

        params = ['Pr']
        T.Plotting_1d(dire, model, extra, datasetl, params,
                      loc='upper left', legcolor=True)
        pylab.xlim(15, 38)
        pylab.tight_layout()
        pylab.savefig('pre1.pdf')
        if show:
            pylab.show()

    if l == 'S3' or l == 'S3d':
        datasetl = ['BBAO', 'BBAO+PlDa']
        model = ['oLCDM']
        extra = 'pre'
        params = ['Om', 'Ol']
        lw = [4, 4, 4]
        lcolor = ['green', 'black']
        T.Plotting_2d(dire, model, extra, datasetl, params, lwidth=lw,
                      fillcolor=lcolor, loc='upper left', Ol='True', legcolor=True, Nb=40)
        pylab.ylabel('$\\Omega_{\Lambda}$')
        pylab.plot([0, 1], [1, 0], 'k:')
        pylab.axis([0.0, 0.5, 0.0, 1.3])
        pylab.tight_layout()
        pylab.savefig('pre2.pdf')
        if show:
            pylab.show()


# H plot

    if l == 'HLadder':

        T.plotH(1.0, dire, 'owaCDM_phy_GBAO+SN+PlRd', 'b',
                '$ow_0w_a$-CDM GBAO+SN+$\mathbf{r_d}$')
        T.plotH(1.2, dire, 'PolyCDM_phy_GBAO+SN+PlRd', 'b',
                '$\\mbox{PolyCDM GBAO+SN}+\mathbf{r_d}$')
        pylab.text(56, 1.1, r'Inverse', fontsize=16)
        pylab.text(56, 1.0, r'Distance Ladder', fontsize=16)

        T.plotH(1.6, '', '/astro/u/anze/work/Planck/PLA/base/planck_lowl_lowLike/base_planck_lowl_lowLike',
                'g', '$\Lambda$CDM Planck (full)')
        T.plotH(1.8, '', '/astro/u/anze/work/Planck/PLA/base/WMAP/base_WMAP',
                'g', '$\Lambda$CDM WMAP (full)')
        pylab.text(56, 1.7, r'CMB+$\Lambda$CDM', fontsize=16)

        T.plotH(2.2, '', 'N: 73.8 2.4', 'r', 'Riess++')
        T.plotH(2.4, '', 'N: 74.3 2.1', 'r', 'Freedman++')
        T.plotH(2.6, '', 'N: 72.5 2.5', 'r', 'Efstathiou')
        pylab.text(56, 2.5, r'Distance Ladder', fontsize=16)

        pylab.xlim(55, 90)
        pylab.ylim(0.8, 2.8)
        pylab.grid(axis='both')
        pylab.yticks([])
        pylab.xlabel("$H_0 [$kms$^{-1}$Mpc$^{-1}]$")
        pylab.tight_layout()
        pylab.savefig('Ayaaaa-boom.pdf')
        if show:
            pylab.show()

    if l == 'HladderOk':

        datasetl = ['GBAO+SN+PlRd']
        model = ['PolyCDM']
        extra = 'phy'
        params = ['h', 'Ok']
        T.Plotting_2d(dire, model, extra, datasetl, params, Nb=30)
        pylab.xlabel('$h$')
        pylab.ylabel('$\\Omega_k$')
        pylab.tight_layout()
        pylab.savefig('hok.pdf')
        if show:
            pylab.show()

    if l == 'mnuh':
        datasetl = ['Planck']
        model = ['nuLCDM']
        extra = 'phy'
        params = ['mnu', 'h']
        T.Plotting_2d(dire, model, extra, datasetl, params, Nb=30)
        pylab.xlabel('$\\sum m_\\nu$')
        pylab.ylabel('$h$')
        pylab.tight_layout()
        pylab.savefig('mnuh.pdf')
        if show:
            pylab.show()
