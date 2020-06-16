#!/usr/bin/env python

# Error bars mainly got from  table 3,
# http://arxiv.org/pdf/1108.2635v1.pdf
# and MGS paper

from simplemc.cosmo.cosmoApprox import rd_cuesta_approx
import matplotlib.pyplot as plt
import numpy as np

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


rd_EHtoCAMB = 153.19/149.28
rd_fid_DR12 = 147.78

fmt1 = '^'
fmt2 = 's'
empty1 = True
empty2 = False
alpha = 1.0

# Errorbars from DR12 Full-shape
fact = (300000./rd_fid_DR12)
rd_fid_DR14_LRG = rd_cuesta_approx(0.0220, 0.1190+0.0220, 0.0006, 3.04)

class Data_plots:
    """
    Cosmological calculator to plot the basic functions
    """
    def __init__(self ):
        self.plaw = .5
        # for division by log(1+z) use this
        # plaw=-1


    # Plotting -  Error bars
    def plot_errorbar(self, z, val, yerr=0, color=0, fmt=0, markersize=0, \
                      label=None, empty=True, alpha=1):
        if empty:
            mfc = 'white'
            lw = 1
        else:
            mfc = color
            lw = 2
        plt.errorbar(z, val/self.fixer(z), yerr=yerr/self.fixer(z), color=color, fmt=fmt,
                     markersize=markersize, lw=lw, capthick=lw, capsize=2+2*lw,
                     markerfacecolor=mfc, alpha=alpha)
        if label:
            if (mfc == 'white'):
                plt.plot([], [], fmt, color='black', label=label,
                           markersize=markersize, markerfacecolor=mfc)
            else:
                plt.plot([], [], fmt, color='black',
                           label=label, markersize=markersize)
        #elif label < 0:
        else:
            if (mfc == 'white'):
                return plt.plot([], [], fmt, color='black', markersize=markersize, markerfacecolor=mfc)
            else:
                return plt.plot([], [], fmt, color='black', markersize=markersize)



    def fixer(self, z):
        if self.plaw > 0:
            return z**self.plaw
        else:
            return np.log(1.+z)



    def ersys(self, x, y):
        return np.sqrt(x**2 + y**2)




    def DaOverrd_data(self):
        self.plot_errorbar(zCombBAO1,  1512.4/rd_fid_DR12,     yerr=self.ersys(22.5, 11.0)/rd_fid_DR12,
                      color='red', fmt='d', markersize=8, empty=empty2, label="$\\rm{BOSS\ Galaxy\ DR12}$")
        self.plot_errorbar(zCombBAO2,  1975.2/rd_fid_DR12,     yerr=self.ersys(26.6, 14.1)/rd_fid_DR12,
                      color='red', fmt='d', markersize=8, empty=empty2)
        self.plot_errorbar(zCombBAO3,  2306.7/rd_fid_DR12,  	 yerr=self.ersys(33.2, 16.7)/rd_fid_DR12,
                      color='red', fmt='d', markersize=8, empty=empty2)

        self.plot_errorbar(zLyaA,  37.77,  yerr=2.13,  color='red', fmt='o', markersize=8,
              label="$\\rm{BOSS}\ \\mathrm{Ly}\\alpha-\\rm{auto}\ \\rm{DR12}$", empty=empty2)

        #https://arxiv.org/abs/1904.03430
        self.plot_errorbar(zLyaC,  36.3,   yerr=1.8,    color='red', fmt='*', markersize=8,
              label="$\\rm{eBOSS}\ \\mathrm{Ly}\\alpha -\\rm{cross}\ \\rm{DR14}$", empty=empty2)
        self.plot_errorbar(0.81,  1.81*10.75,   yerr=0.43,    color='red',
              fmt='x', markersize=8, label="$\\rm{DES} Y1$", empty=empty2)
        #plt.legend(loc='best', numpoints=1, frameon=False, fontsize=13)
        #plt.show()




    def HIOverrd_data(self):
        self.plot_errorbar(zCombBAO1,  fact*zCombBAO1/81.21,
                    yerr=fact*zCombBAO1*self.ersys(2.17, 0.97)/(81.21)**2,
                    color='green', fmt='d', markersize=8, empty=empty2, label="$\\rm{BOSS\ Galaxy\ DR12}$")
        self.plot_errorbar(zCombBAO2,  fact*zCombBAO2/90.90,
                    yerr=fact*zCombBAO2*self.ersys(2.07, 1.08)/(90.90)**2,
                    color='green', fmt='d', markersize=8, empty=empty2)
        self.plot_errorbar(zCombBAO3,  fact*zCombBAO3/98.96,
                    yerr=fact*zCombBAO3*self.ersys(2.21, 1.18)/(98.96)**2,
                    color='green', fmt='d', markersize=8, empty=empty2)

        #https://arxiv.org/abs/1904.03430
        self.plot_errorbar(zLyaA,  9.07*zLyaA,       yerr=0.31*zLyaA,
              color='green', fmt='-o', markersize=8, empty=empty2,
                label="$\\rm{BOSS}\ \\mathrm{Ly}\\alpha-\\rm{auto}\ \\rm{DR12}$")
        self.plot_errorbar(zLyaC,  9.2*zLyaC,        yerr=0.36*zLyaC,
              color='green', fmt='-*', markersize=8, empty=empty2,
                label="$\\rm{eBOSS}\ \\mathrm{Ly}\\alpha -\\rm{cross}\ \\rm{DR14}$")
        #plt.legend(loc='best', numpoints=1, frameon=False, fontsize=13)
        #plt.show()




    def DVOverrd_data(self):
        self.plot_errorbar(z6dFGS,    2.97*rd_EHtoCAMB,   yerr=rd_EHtoCAMB*0.015/0.336**2,
              color='blue', fmt='>', markersize=6, empty=empty2, label="$\\rm{6dFGS}$", alpha=alpha)
        self.plot_errorbar(zMGS,      4.464,    yerr=0.168,               color='blue',
              fmt='<', markersize=6, label="$\\rm{SDSS\ MGS}$", empty=empty2, alpha=alpha)


        self.plot_errorbar(zSDSS1,    5.2493*rd_EHtoCAMB, yerr=rd_EHtoCAMB*0.0061/0.1905 **
              2, color='blue', fmt=fmt1, markersize=6, empty=empty1, alpha=alpha)
        self.plot_errorbar(zSDSS2,    1348./rd_fid_DR12, yerr=26./rd_fid_DR12, color='blue',
              fmt=fmt1, markersize=6, label="$\\rm{SDSS\ DR7}$", empty=empty1, alpha=alpha)

        self.plot_errorbar(zWiggleZ1, 1695./rd_fid_DR12, yerr=82./rd_fid_DR12, color='blue',
              fmt=fmt2, markersize=6, label="$\\rm{WiggleZ}$", empty=empty1, alpha=alpha)
        self.plot_errorbar(zWiggleZ2, 2194./rd_fid_DR12, yerr=100./rd_fid_DR12,
              color='blue', fmt=fmt2, markersize=6, empty=empty1, alpha=alpha)
        self.plot_errorbar(zWiggleZ3, 2486./rd_fid_DR12, yerr=85./rd_fid_DR12,
              color='blue', fmt=fmt2, markersize=6, empty=empty1, alpha=alpha)

        self.plot_errorbar(0.72,  2353/rd_fid_DR14_LRG, yerr=62/rd_fid_DR12,  color='blue',
              fmt='s', markersize=8, label="$\\rm{eBOSS\ LRG\ DR14}$", empty=False, alpha=alpha)
        self.plot_errorbar(1.52,  3843/rd_fid_DR12, yerr=147/rd_fid_DR12,  color='blue', fmt='d',
              markersize=8, label="$\\rm{eBOSS\ QSO\ DR14}$", empty=False, alpha=alpha)

        #plt.legend(loc='best', numpoints=1, frameon=False, fontsize=13)
        #plt.show()



    def Hubble_data(self):
        dataHz = np.loadtxt('simplemc/data/Hz_all.dat')
        redshifts, obs, errors = [dataHz[:,i] for i in [0,1,2]]
        plt.errorbar(redshifts, obs, errors, xerr=None,
            color='purple', marker='o', ls='None',
            elinewidth =2, capsize=5, capthick = 1, label=None)
        #plt.show()


    def SNIa_data(self):
        dataSN = np.loadtxt('simplemc/data/pantheon_lcparam_full_long_zhel.txt', usecols=[1,2,3, 4, 5])

        redshifts, obs, errors = [dataSN[:,i] for i in [0,3,4]]

        plt.errorbar(redshifts, obs, errors, xerr=None,
            color='purple', marker='o', ls='None',
            elinewidth =2, capsize=5, capthick = 1, label=None)
        plt.ylim(18, 30)
        #names = ['name', 'z', 'mu', 'error']
        #result = pd.read_table('simplemc/data/sn_z_mu_dmu_union2.txt', sep='\s+', names=names, index_col='z')
        #result=result.sort_index()
        #plt.figure(figsize=(14,7))
        #result['mu'].plot(yerr=result['error'], linestyle='None', label = 'SN')

        #hub_CDM    = self.logatoz(self.hubble(self.lna))
        #mu = self.mu(self.zvals, hub_CDM)
        #plt.plot(self.zvals, mu, 'o',  markersize=2, label = '$\Omega_{DM} = 0.24, \Omega_\Lambda=0.76$')

        #mu_SF = self.mu(self.zvals, self.hub_SF_z)
        #plt.plot(self.zvals, mu_SF, label = 'SF', color = 'r')

        #plt.xlabel(r'Corrimiento al rojo - z', fontsize = 20)
        #plt.ylabel(r'Distancia modular - $\mu$', fontsize = 20)
        #plt.title('Supernova Tipo Ia')
        #plt.legend(loc='lower right', frameon=False)
        #plt.savefig('SN_models.pdf')
        #plt.show()

    def fs8_data(self):
        datafs8 = np.loadtxt('simplemc/data/Growth_tableII.txt')
        redshifts, obs, errors = [datafs8[:,i] for i in [0,1,2]]
        plt.errorbar(redshifts, obs, errors, xerr=None,
            color='purple', marker='o', ls='None',
            elinewidth =2, capsize=5, capthick = 1)




if __name__=='__main__':
    D= Data_plots()
    #D.DaOverrd_data()
    #D.HIOverrd_data()
    #D.DVOverrd_data()
    #D.Hubble_data()
    D.SNIa_data()
    #D.fs8_data()