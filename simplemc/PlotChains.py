#!/usr/bin/env python

from simplemc.plots.Simple_Plots import Simple_plots
from simplemc.plots.Plot_elipses import plot_elipses
import matplotlib.pyplot as plt
import numpy as np
import sys

dir_name   = 'simplemc/chains/'
roots      = ['LCDM_phy_HD_mcmc']
params_1D  = ['h', 'Om']
params_2D  = [['h', 'Om']]
labels     = ['HD+SN']



#roots = ['GPantheon_phy_CPantheon_mcmc']
#params_1D = ['zbin%d'%i for i in range(20)]
#labels     = ['Compress_Pantheon']

S = Simple_plots(dir_name, roots, labels)
mean, cov = S.Covariance(params_1D)
print(mean, cov)
plotter = 'getdist'
#Simple_plots, getdist, corner, fgivenx




#1D, 2D and triangular posterior distributions
if plotter == 'Simple_plots':
    S = Simple_plots(dir_name, roots, labels)
    #S.Show_limits(params_1D)
    #S.Covariance(params_1D)
    #S.Plots1D(params_1D)
    S.Plots2D(params_2D)
    #S.triangle(params_1D)


elif plotter == 'corner':
    S = Simple_plots(dir_name, roots)
    S.cornerPlotter(params_1D)


elif plotter == 'getdist':
    sys.path = ['getdist', 'corner'] + sys.path
    from getdist import plots
    fig, ax = plt.subplots()
    g = plots.getSubplotPlotter(chain_dir= dir_name, width_inch=10,
                                analysis_settings={'ignore_rows': 0.2})
    #g.plots_1d(roots, params=params_1D)
    g.plots_2d(roots, param_pairs=params_2D, nx=1, filled=True)
    fig = plt.gcf()
    ax = fig.gca()
    plot_elipses(mean, cov, 0, 1, ax=ax)
    #g.add_legend(labels,  legend_loc='best', ax=ax)
    #g.triangle_plot(roots, params_1D, filled=True) #, plot_3d_with_param='h')
    import matplotlib.patches as mpatches
    red_patch = mpatches.Patch(color='green', label='Fisher')
    blue_patch = mpatches.Patch(color='blue', label='MCMC')

    plt.legend(handles=[red_patch, blue_patch])
    g.export('Plot_Fisher.pdf')
    plt.show()



elif plotter == 'fgivenx':
    S = Simple_plots(dir_name, roots[0])

    z = np.linspace(0,4,100)
    def func(z,theta1):
        Omega_m, h = theta1
        Hz=100*h*(Omega_m*(1+z)**3 + (1-Omega_m))**0.5
        return Hz

    S.fgivenx(['Om', 'h'], z, func, labels=['z','H(z)'])



