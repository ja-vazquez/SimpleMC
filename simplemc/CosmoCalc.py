#!/usr/bin/env python



from simplemc.runbase import ParseModel
import matplotlib.pyplot as plt
import numpy as np
import sys


class CosmoCalc:
    """
    Cosmological calculator to plot the basic functions
    """
    def __init__(self, model, funct=None, param=None, minparam=None,
                    maxparam=None, nsteps=4, savepdf=False, zmax=3):
        self.T     = ParseModel(model)
        self.model = model
        self.funct = funct
        self.param = param
        if param:
            assert(maxparam > minparam)
            self.minparam = minparam
            self.maxparam = maxparam
            for i, par in enumerate(self.T.freeParameters()):
                if self.param == par.name: self.Opar = par
                else: continue

        self.nsteps  = nsteps
        self.savepdf = savepdf

        self.zl = np.linspace(0, zmax, 50)

        self.func_dic={'Hubble': self.Hubble, 'DaOverrd': self.DaOverrd,
                       'HIOverrd': self.HIOverrd, 'DVOverrd': self.DVOverrd,
                       'HubInvOverz': self.HubInvOverz, 'SNIa': self.SNIa, 'fs8':self.fs8,
                       'Age': self.Age}
        try:
            self.function  = self.func_dic[self.funct]
        except:
            print('function not in list:', list(self.func_dic.keys()))
            sys.exit(1)
        print (model)



    def run_plot(self, **kwargs):

        if self.param is None:
            y=[self.function(z) for z in self.zl]
            self.plot_vars(y, ylabel=self.ylabel, **kwargs)

        else:
            for j in np.linspace(self.minparam, self.maxparam, self.nsteps):
                self.Opar.setValue(j)
                self.T.updateParams([self.Opar])
                y= [self.function(z) for z in self.zl]
                self.plot_vars(y, ylabel=self.ylabel,
                               label='{}={:.2f}'.format(self.param, j), **kwargs)

        if self.savepdf: plt.savefig('sm_{}_{}.pdf'.format(self.model, self.funct))
        plt.show()


    #--- CC --
    def Hubble(self, z):
        function = 100*self.T.h*np.sqrt(self.T.RHSquared_a(1./(1+z)))
        self.ylabel   = "$H(z)$"
        return function


    #--- BAO --
    def DaOverrd(self, z):
        function = 1*self.T.DaOverrd(z)
        self.ylabel = "$D_M(z)/r_d$"
        return function


    #--- BAO --
    def HIOverrd(self, z):
        function = z*self.T.HIOverrd(z)
        self.ylabel="$zD_H(z)$"
        return function


    #--- BAO --
    def DVOverrd(self, z):
        function = 1.*self.T.DVOverrd(z)
        self.ylabel ="$D_v(z)/r_d$"
        return function


    def HubInvOverz(self, z):
        function = 100*self.T.h*np.sqrt(self.T.RHSquared_a(1./(1+z)))/(1+z)
        self.ylabel="$H(z)/(1+z)$"
        return function


    #--- SN --
    def SNIa(self, z):
        function = self.T.distance_modulus(z) + 43
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




if __name__ == "__main__":
    #C = CosmoCalc('LCDM', 'fs8')
    #C.run_plot(lw='2')

    #C = CosmoCalc('LCDM', 'DaOverrd', 'h', 0.4, 0.9)
    #C.run_plot(lw='2')

    C = CosmoCalc('owaCDM', 'fs8', 'wa', -0.5, 0.5, 5, zmax=0.8)
    C.run_plot(lw='1')

    #C = CosmoCalc('LCDM', 'Age')
    #C.Age()

