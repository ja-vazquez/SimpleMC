#!/usr/bin/env python



from simplemc.runbase import ParseModel
import matplotlib.pyplot as plt
import numpy as np
import sys


class CosmoCalc:
    """
    Cosmological calculator to plot the basic functions
    """
    def __init__(self, model, funct, param=None, minparam=None, maxparam=None,
                    nsteps=4, savepdf=False):
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

        self.zl = np.arange(0, 3, 0.1)

        print (model)



    def run_plot(self, **kwargs):
        if 'Hubble' in self.funct:
            function = self.Hubble2
            ylabel   = "$H(z)$"

        if self.param is None:
            y=[function(z) for z in self.zl]
            self.plot_vars(y, ylabel=ylabel, **kwargs)

        else:
            for j in np.linspace(self.minparam, self.maxparam, self.nsteps):
                self.Opar.setValue(j)
                self.T.updateParams([self.Opar])
                y= [function(z) for z in self.zl]
                self.plot_vars(y, ylabel=ylabel,
                               label='{}={:.2f}'.format(self.param, j), **kwargs)


        if self.savepdf: plt.savefig('sm_{}_{}.pdf'.format(self.model, self.funct))
        plt.show()



    #--- BAO --
    def DaOverrd(self, **kwargs):
        function = self.T.DaOverrd
        self.plot_vars(function, ylabel="$D_M(z)/r_d$", **kwargs)
        return True


    #--- BAO --
    def HIOverrd(self, **kwargs):
        def zHI(z):
            return z*self.T.HIOverrd(z)
        function = zHI
        self.plot_vars(function, ylabel="$zD_H(z)$", **kwargs)
        return True


    #--- BAO --
    def DVOverrd(self, **kwargs):
        function = self.T.DVOverrd
        self.plot_vars(function, ylabel="$D_v(z)/r_d$", **kwargs)


    def HubInvOverz(self, **kwargs):
        def hub(z):
            return np.sqrt(self.T.RHSquared_a(1./(1+z)))/(1+z)
        function = hub
        self.plot_vars(function, ylabel="$H(z)/(1+z)$", **kwargs)
        return True

    #--- CC --
    def Hubble2(self, z):
        #def hub(z):
        return self.T.h*np.sqrt(self.T.RHSquared_a(1./(1+z)))
        #function = hub
        #self.plot_vars(function, ylabel="$H(z)$", **kwargs)
        #return function


    #--- CC --
    def Hubble(self, **kwargs):
        def hub(z):
            return np.sqrt(self.T.RHSquared_a(1./(1+z)))
        function = hub
        #self.plot_vars(function, ylabel="$H(z)$", **kwargs)
        return function


    #--- SN --
    def SNIa(self, **kwargs):
        def SN(z):
            return self.T.distance_modulus(z) + 43
        function = SN
        self.plot_vars(function, ylabel="$d_L(z)$", **kwargs)
        return True


    #--- fs8 --
    def fs8(self, **kwargs):
        function = self.T.fs8
        self.plot_vars(function, ylabel="$f\sigma_8$", **kwargs)
        return True


    #--- Age --
    def Age(self):
        AgeofUniverse= self.T.Age()
        print('{:.2f}Gys'.format(AgeofUniverse))
        return AgeofUniverse



    def plot_vars(self, y, ylabel=None, **kwargs):
        #tmp = [funct(z) for z in self.zl]
        plt.plot(self.zl, y, **kwargs)
        plt.ylabel(ylabel, fontsize=20)
        plt.xlabel("$z$", fontsize=20)
        plt.title(self.model)
        plt.legend(loc='best')
        plt.grid()
        #if self.savepdf: plt.savefig('sm_{}.pdf'.format(str(funct)))
        #plt.show()




if __name__ == "__main__":
    #C = CosmoCalc('LCDM', 'DaOverrd')
    #C.run_plot(lw='2')

    C = CosmoCalc('LCDM', 'Hubble') #, 'h', 0.4, 0.9)
    C.run_plot(lw='2')


"""
T = LCDMCosmology(mnu=0)

zl = arange(0, 3, 0.1)
pylab.figure(figsize=(8, 10))


pylab.subplot(3, 1, 1)
y1 = [T.DaOverrd(z) for z in zl]
pylab.plot(zl, y1, 'k-')
pylab.ylabel("$D_a(z)/r_d$")

pylab.subplot(3, 1, 2)
y1 = [T.HIOverrd(z) for z in zl]
pylab.errorbar(zCMASS, 20.75, yerr=0.73, color='red', fmt='-o')
pylab.errorbar(zLyaA, 9.18,  yerr=0.28, color='blue', fmt='-o')
pylab.errorbar(zLyaC, 9.0,   yerr=0.3, color='magenta', fmt='-o')
pylab.plot(zl, y1, 'k-')
pylab.ylabel("$H^{-1}(z)/r_d$")

pylab.subplot(3, 1, 3)
y1 = [T.DVOverrd(z) for z in zl]
pylab.errorbar(zLOWZ, 8.467, yerr=0.167, color='green', fmt='-o')
pylab.plot(zl, y1, 'k-')
pylab.ylabel("$D_v(z)/r_d$")

pylab.plot([], [], 'g-', label='LOWZ')
pylab.plot([], [], 'r-', label='CMASS')
pylab.plot([], [], 'b-', label='Lyman-$\\alpha$ auto')
pylab.plot([], [], 'magenta', label='Lyman-$\\alpha$ cross')
pylab.legend(loc='lower right')

pylab.xlabel("z")
pylab.savefig("Fig1.pdf")
pylab.show()
"""
