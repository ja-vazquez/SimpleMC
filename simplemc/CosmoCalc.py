#!/usr/bin/env python


from simplemc.plots.Data_plots import Data_plots
from simplemc.runbase import ParseModel
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import sys


class CosmoCalc:
    """
    Cosmological calculator to plot the basic functions
    """
    def __init__(self, model, funct=None, param=None, minparam=None,
                    maxparam=None, nsteps=3, savepdf=False, zmax=3,
                    plot_data=False):
        self.T     = ParseModel(model)
        self.model = model
        self.funct = funct
        self.param = param
        self.Opar  = None
        if param:
            assert(maxparam > minparam)
            self.minparam = minparam
            self.maxparam = maxparam
            for i, par in enumerate(self.T.freeParameters()):
                if self.param == par.name: self.Opar = par
            if self.Opar is None:
                print('No param within model')
                sys.exit(1)

        self.nsteps    = nsteps
        self.savepdf   = savepdf
        self.plot_data = plot_data

        self.plaw = 0.5
        self.zmax = zmax
        self.zl = np.linspace(0.05, self.zmax, 50)

        self.func_dic={'Hubble': self.Hubble, 'DaOverrd': self.DaOverrd,
                       'HIOverrd': self.HIOverrd, 'DVOverrd': self.DVOverrd,
                       'HubInvOverz': self.HubInvOverz, 'SNIa': self.SNIa,
                       'fs8':self.fs8, 'Age': self.Age}
        try:
            self.function  = self.func_dic[self.funct]
        except:
            print('function not in list:', list(self.func_dic.keys()))
            sys.exit(1)
        print (model)



    def run_plot(self, **kwargs):

        fig = plt.figure(figsize=(9, 6))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xscale('log')

        if self.plot_data: self.selec_data(self.funct)
        if self.param is None:
            y= [self.function(z) for z in self.zl]
            plt.plot(self.zl, y, **kwargs)

        else:
            for j in np.linspace(self.minparam, self.maxparam, self.nsteps):
                self.Opar.setValue(j)
                self.T.updateParams([self.Opar])
                y= [self.function(z) for z in self.zl]
                label = '{}={:.2f}'.format(self.param, j)
                plt.plot(self.zl, y, label=label, **kwargs)

        plt.grid()
        plt.title(self.model)
        plt.xlim(0.05, self.zmax+0.1)
        plt.xlabel("$z$", fontsize=25)
        plt.ylabel(self.ylabel, fontsize=20)
        plt.legend(loc='best', numpoints=1, frameon=False, fontsize=13)

        # Axis
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

        ax.xaxis.set_minor_formatter(plt.ScalarFormatter())
        ax.xaxis.set_major_locator(plt.FixedLocator([0.1, 0.5]))
        ax.xaxis.set_minor_locator(plt.FixedLocator([0.2, 1.0, 2]))
        #plt.yticks(list(range(0, 50, 10)))

        if self.savepdf: plt.savefig('sm_{}_{}.pdf'.format(self.model, self.funct))
        plt.show()


    def selec_data(self, name):
        D = Data_plots()
        if name == 'Hubble': D.Hubble_data()
        if name == 'DaOverrd' : D.DaOverrd_data()
        if name == 'DVOverrd' : D.DVOverrd_data()
        if name == 'HIOverrd': D.HIOverrd_data()
        if name == 'SNIa': D.SNIa_data()
        if name == 'fs8': D.fs8_data()
        else: pass

    #--- CC --
    def Hubble(self, z):
        function = 100*self.T.h*np.sqrt(self.T.RHSquared_a(1./(1+z)))
        self.ylabel   = r'$H(z) [km/s Mpc^{-1}]$'

        return function


    #--- BAO --
    def DaOverrd(self, z):
        function = 1*self.T.DaOverrd(z)/self.fixer(z)
        self.ylabel = "$D_M(z)/r_d \sqrt{z}$"

        return function


    #--- BAO --
    def HIOverrd(self, z):
        function = z*self.T.HIOverrd(z)/self.fixer(z)
        self.ylabel="$zD_H(z)/r_d \sqrt{z}$"

        return function


    #--- BAO --
    def DVOverrd(self, z):
        function = 1.*self.T.DVOverrd(z)/self.fixer(z)
        self.ylabel ="$D_v(z)/r_d \sqrt{z}$"

        return function


    def HubInvOverz(self, z):
        function = 100*self.T.h*np.sqrt(self.T.RHSquared_a(1./(1+z)))/(1+z)
        self.ylabel="$H(z)/(1+z)$"
        return function


    #--- SN --
    def SNIa(self, z):
        function = self.T.distance_modulus(z) + 23.8
        self.ylabel="$d_L(z)$"

        return function



    #--- fs8 --
    def fs8(self, z):
        function = 1.*self.T.fs8(z)
        self.ylabel="$f\sigma_8(z)$"

        return function


    #--- Age --
    def Age(self):
        AgeofUniverse= self.T.Age()
        print('{:.2f}Gys'.format(AgeofUniverse))
        return AgeofUniverse



    def plot_vars(self, y, ylabel=None, **kwargs):
        plt.plot(self.zl, y, **kwargs)
        plt.ylabel(ylabel, fontsize=20)
        plt.xlabel("$z$", fontsize=20)
        plt.title(self.model)
        plt.legend(loc='best')
        plt.grid()


    def fixer(self, z):
        if self.plaw > 0:
            return z**self.plaw
        else:
            return np.log(1.+z)


if __name__ == "__main__":
    #C = CosmoCalc('LCDM', 'fs8')
    #C.run_plot(lw='2')

    #C = CosmoCalc('LCDM', 'DaOverrd', 'h', 0.4, 0.9)
    #C.run_plot(lw='2')

    C = CosmoCalc('owaCDM', 'fs8', 'wa', -0.5, 0.5, 5, zmax=0.8)
    C.run_plot(lw='1')

    #C = CosmoCalc('LCDM', 'Age')
    #C.Age()

