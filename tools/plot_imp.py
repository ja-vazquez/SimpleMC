#!/usr/bin/env python
import sys
from Simple_Plots import Simple_Plots
import matplotlib.pyplot as plt


dir_name   = 'chains/'
roots      = ['wCDM_phy_BBAO+JLA'] #, 'wCDM_phy_BBAO+Pantheon+Planck_15']
params_1D  = ['h', 'w', 'Ol', 'Age']
params_2D  = [['h', 'w'], ['Om', 'h']]
labels     = ['BBAO+JLA'] #, 'BBAO+Pantheon+PLK15']


plotter = 'corner'
#[Simple_plots, getdist, corner]


if plotter == 'Simple_plots':
    S = Simple_Plots(dir_name, roots, labels)
    #S.Show_limits(params_1D)
    #S.Covariance(params_1D)
    #S.Plots1D(params_1D)
    #S.Plots2D(params_2D)
    S.triangle(params_1D)


elif plotter == 'getdist':
    sys.path = ['getdist', 'corner'] + sys.path
    from getdist import plots
    g = plots.getSubplotPlotter(chain_dir= dir_name, width_inch=10,
                                analysis_settings={'ignore_rows': 0.2})
    #g.plots_1d(roots, params=params_1D)
    #g.plots_2d(roots, param_pairs=params_2D, nx=2, filled=True)
    g.triangle_plot(roots, params_1D, filled=True) #, plot_3d_with_param='h')
    g.add_legend(labels,  legend_loc='best')
    g.export('Plot_getdist.pdf')
    plt.show()


elif plotter == 'corner':
    S = Simple_Plots(dir_name, roots)
    S.cornerPlotter(params_1D)


