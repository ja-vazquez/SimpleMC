#!/usr/bin/env python


from simplemc.plots.Simple_Plots import Simple_plots
from simplemc.plots.Plot_elipses import plot_elipses
import matplotlib.pyplot as plt
import numpy as np
import sys

#dir_name   = 'simplemc/chains/'
#roots      = ['Anisotropic_phy_Planck_15+CBAO+SN+HD_nested_dynesty_multi',
#        'NumuLCDM_phy_Planck_15+CBAO+SN+HD_nested_dynesty_multi']
params_1D  = ['h', 'Om']
params_2D  = [['h', 'Om']]
labels     = ['HD+SN']
dir_name = 'simplemc/chain_not_github/'
roots = ['LCDM_phy_HD+SN_mcmc']

#roots = ['GPantheon_phy_CPantheon_mcmc']
#params_1D = ['zbin%d'%i for i in range(20)]
#labels     = ['Compress_Pantheon']

#S = Simple_plots(dir_name, roots, labels)
#mean, cov = S.Covariance(params_1D)
#print(mean, cov)
#plotter = 'getdist' #'Simple_plots'
plotter = 'fgivenx4'
#Simple_plots, getdist, corner, fgivenx




#1D, 2D and triangular posterior distributions
if plotter == 'Simple_plots':
    S = Simple_plots(dir_name, roots, labels)
    #S.Show_limits(params_1D)
    #S.Covariance(params_1D)
    S.Plots1D(params_1D)
    #S.Plots2D(params_2D)
    #S.triangle(params_1D)


elif plotter == 'corner':
    S = Simple_plots(dir_name, roots)
    S.cornerPlotter(params_1D)


elif plotter == 'getdist':
    sys.path = ['getdist', 'corner'] + sys.path
    from getdist import plots
    fig, ax = plt.subplots()
    g = plots.getSubplotPlotter(chain_dir= dir_name, width_inch=10,
            analysis_settings={'ignore_rows': 0.2, 'smooth_scale_2D': -1})
    #g.plots_1d(roots, params=params_1D)
    g.plots_2d(roots, param_pairs=params_2D, nx=1, filled=True)
    #fig = plt.gcf()
    #ax = fig.gca()
    #plot_elipses(mean, cov, 0, 1, ax=ax)
    ##g.add_legend(labels,  legend_loc='best', ax=ax)
    ##g.triangle_plot(roots, params_1D, filled=True) #, plot_3d_with_param='h')
    #import matplotlib.patches as mpatches
    #red_patch = mpatches.Patch(color='green', label='Fisher')
    #blue_patch = mpatches.Patch(color='blue', label='MCMC')

    #plt.legend(handles=[red_patch, blue_patch])
    g.export('Plot_LCDM.pdf')
    plt.show()



elif plotter == 'fgivenx':
    S = Simple_plots(dir_name, roots[0], nchains=[1])

    z = np.linspace(0,4,100)
    def func(z,theta1):
        Omega_m, h = theta1
        Hz=100*h*(Omega_m*(1+z)**3 + (1-Omega_m))**0.5
        return Hz

    S.fgivenx(['Om', 'h'], z, func, labels=['z','H(z)'])




elif plotter == 'fgivenx2':
    # e-field Dark energy
    from simplemc.models.PhiCDMCosmology import PhiCosmology
    from simplemc.cosmo.paramDefs import *

    P = PhiCosmology(alpha=1)


    from fgivenx import plot_contours, samples_from_getdist_chains

    params = ['h', 'Om','Ok', 'phimu', 'philam', 'phibeta']
    #dir ='/Users/josevazquezgonzalez/Desktop/Codigos_SimpleMC/new_SimpleMC/SimpleMC/simplemc/chains/'
    dir='/Users/josevazquezgonzalez/Desktop/Codigos_SimpleMC/SimpleMC/chains_epsilon/'
    #file_root = dir+'LCDM_phy_HD+SN_nested_dynesty_multi'
    #file_root=dir + 'Phi_exp_pow2_curv_phy_HD+SN+CBAO+Planck_15_nested_dynesty_multi'
    file_root=dir + 'Phi_pow_exp_curvf_phy_HD+SN+CBAO+Planck_15_nested_dynesty_multi'
    try:
        samples, weights = samples_from_getdist_chains(params, file_root)
    except:
        samples, weights = samples_from_getdist_chains(params, file_root+'_')


    def func(z,theta1):
        a = 1./(1+z)
        h, Om, Ok, phimu, philam, phibeta = theta1
        #print (theta1)
        h_  = h_par
        Om_  = Om_par
        Ok_  = Ok_par
        phimu_ = phimu_par
        philam_ = philam_par
        phibeta_ = phibeta_par

        h_.setValue(h)
        Om_.setValue(Om)
        Ok_.setValue(Ok)
        phimu_.setValue(phimu)
        philam_.setValue(philam)
        phibeta_.setValue(phibeta)

        P.updateParams([h_, Om_, Ok_, phimu_, philam_, phibeta_])
        w = P.w_de(a)
        #print(w)
        #Omega_m, h = theta1
        #Hz=100*h*(Omega_m*(1+z)**3 + (1-Omega_m))**0.5
        return w

    z = np.linspace(0,3,100)

    cbar = plot_contours(func, z, samples, weights=weights,
                         contour_line_levels=[1,2], linewidths = 0.8,
                         colors=plt.get_cmap('winter'), fineness=0.08)
    cbar = plt.colorbar(cbar, ticks=[0, 1, 2])
    cbar.set_ticklabels(['','$1\sigma$','$2\sigma$'])
    plt.grid()
    plt.xlabel('z', fontsize=20)
    plt.ylabel('w(z)', fontsize=20)
    plt.ylim(-1.15, -0.76)
    plt.tight_layout()
    plt.savefig('pow_exp_curv.pdf')
    plt.show()
#    S.fgivenx(['Om', 'h'], z, func, labels=['z','H(z)'])


elif plotter == 'fgivenx3':
    #Simple-Graduated DE
    from simplemc.models.PhiCDMCosmology import PhiCosmology
    from simplemc.cosmo.paramDefs import *
    from fgivenx import plot_contours, samples_from_getdist_chains


    dir='/Users/josevazquezgonzalez/Desktop/Codigos_SimpleMC/new_SimpleMC/SimpleMC/simplemc/chains/'
    #file_root=dir + 'Grad_phy_CBAO+HD+SN_nested_dynesty_multi'
    file_root=dir + 'Grad_Ok_phy_CBAO+HD+SN+Planck_15_nested_dynesty_multi'
    #file_root=dir + 'Grad_Ok_phy_CBAO+HD+SN_nested_dynesty_multi'
    params = ['h', 'Om', 'Ok', 'ggama']
    #params = ['h', 'Om', 'ggama']

    try:
        samples, weights = samples_from_getdist_chains(params, file_root)
    except:
        samples, weights = samples_from_getdist_chains(params, file_root+'_')


    z = np.linspace(0, 3, 100)
    def func(z,theta1):
        h, Om, Ok, ggama = theta1
        #h, Om, ggama = theta1
        #Ok=0
        Hz=100*h*(Om*(1+z)**3 + (1-Om-Ok)*(1+3*ggama*np.log(1+z)) + Ok*(1+z)**2)**0.5
        return Hz/(1+z)



    cbar = plot_contours(func, z, samples, weights=weights,
                         contour_line_levels=[1,2], linewidths = 0.8,
                         colors=plt.get_cmap('GnBu_r'), fineness=0.08)
    cbar = plt.colorbar(cbar, ticks=[0, 1, 2])
    cbar.set_ticklabels(['','$1\sigma$','$2\sigma$'])
    plt.grid()
    plt.title('oDE - BAO+SN+H+PLK')
    plt.xlabel('z', fontsize=20)
    plt.ylabel('H(z)/(1+z)', fontsize=20)


    rd_fid_DR12 = 147.78

    zLyaA = 2.37
    zLyaC = 2.35
    zCombBAO1 = 0.38
    zCombBAO2 = 0.51
    zCombBAO3 = 0.61
    fact = (300000./rd_fid_DR12)

    def ersys(x, y):
        return np.sqrt(x**2 + y**2)

    #https://arxiv.org/abs/1904.03430
    plt.errorbar(zLyaC,  fact/9.20/(1+zLyaC),        yerr=0.36*fact/(1+zLyaC)/9.20**2,
                  color='red', fmt='s', markersize=6, elinewidth=1.5)
    #https://arxiv.org/abs/1904.03430
    plt.errorbar(zLyaA,  fact/8.86/(1+zLyaA),        yerr=0.29*fact/(1+zLyaA)/8.86**2,
                  color='red', fmt='s', markersize=6, elinewidth=1.5)

    plt.errorbar(zCombBAO1,  81.21/(1+zCombBAO1),       yerr=ersys(2.17, 0.97)/(1+zCombBAO1),
                  color='red', fmt='d', markersize=6, elinewidth=1.5)
    plt.errorbar(zCombBAO2,  90.90/(1+zCombBAO2),       yerr=ersys(2.07, 1.08)/(1+zCombBAO1),
                  color='red', fmt='d', markersize=6, elinewidth=1.5)
    plt.errorbar(zCombBAO3,  98.96/(1+zCombBAO3),       yerr=ersys(2.21, 1.18)/(1+zCombBAO1),
                  color='red', fmt='d', markersize=6, elinewidth=1.5)
    plt.errorbar(0.01,  69.8,       yerr=0.8,
                  color='red', fmt='d', markersize=6, elinewidth=1.5)

    plt.plot(z, [100*0.682*(0.30*(1+zi)**3 + (1-0.30))**0.5/(1+zi) for zi in z], 'red', lw=1.5)
    #plt.ylim(55, 80)
    plt.ylim(50, 85)
    plt.tight_layout()
    plt.savefig('GDE_Hz_PLK.pdf')
    plt.show()
#    S.fgivenx(['Om', 'h'], z, func, labels=['z','H(z)'])



elif plotter == 'fgivenx4':
    #new- LCDM model
    from simplemc.models.PhiCDMCosmology import PhiCosmology
    from simplemc.cosmo.paramDefs import *
    from fgivenx import plot_contours, samples_from_getdist_chains
    import matplotlib as mpl

    dir_name ='/Users/josevazquezgonzalez/Desktop/Desktop_Jose/work/Papers/Graduated/new_gDE/'
    file_root=dir_name + 'z2.3_LCDM_PLK18_BAO'
    #params = ['H0', 'omegam', 'w']
    params = ['H0', 'omegam']

    samples, weights = samples_from_getdist_chains(params, file_root)


    z = np.linspace(0, 3.5, 100)
    def func(z,theta1):
        #h, Om, zs= theta1
        h, Om = theta1

        Ode=[(1-Om) if b else -(1-Om) for b in z<2.32]
        #Ode=[(1-Om) if b else -(1-Om) for b in z<zs]

        Hz=h*(Om*(1+z)**3 + Ode)**0.5
        return Hz/(1+z)

    default_color_levels = np.arange(0,2.001, 0.01)
    cbar = plot_contours(func, z, samples, weights=weights,
                         contour_line_levels=[1,2], linewidths = 0.8,
                         colors=plt.get_cmap('GnBu_r'), fineness=0.08, contour_color_levels=default_color_levels)
    cbar = plt.colorbar(cbar, ticks=[0, 1, 2])
    cbar.set_ticklabels(['$0\sigma$','$1\sigma$','$2\sigma$'])
    plt.grid()
    plt.title('new-$\Lambda$CDM (z=2.32) - PLK18+BAO')
    #plt.title('new-$\Lambda$CDM - PLK18')
    plt.xlabel('z', fontsize=20)
    plt.ylabel('H(z)/(1+z)', fontsize=20)

    rd_fid_DR12 = 147.78

    zLyaA = 2.37
    zLyaC = 2.35
    zCombBAO1 = 0.38
    zCombBAO2 = 0.51
    zCombBAO3 = 0.61
    fact = (300000./rd_fid_DR12)

    def ersys(x, y):
        return np.sqrt(x**2 + y**2)

    #https://arxiv.org/abs/1904.03430
    plt.errorbar(zLyaC,  fact/9.20/(1+zLyaC),        yerr=0.36*fact/(1+zLyaC)/9.20**2,
                  color='red', fmt='s', markersize=6, elinewidth=1.5)
    #https://arxiv.org/abs/1904.03430
    plt.errorbar(zLyaA,  fact/8.86/(1+zLyaA),        yerr=0.29*fact/(1+zLyaA)/8.86**2,
                  color='red', fmt='s', markersize=6, elinewidth=1.5)

    plt.errorbar(zCombBAO1,  81.21/(1+zCombBAO1),       yerr=ersys(2.17, 0.97)/(1+zCombBAO1),
                  color='red', fmt='d', markersize=6, elinewidth=1.5)
    plt.errorbar(zCombBAO2,  90.90/(1+zCombBAO2),       yerr=ersys(2.07, 1.08)/(1+zCombBAO1),
                  color='red', fmt='d', markersize=6, elinewidth=1.5)
    plt.errorbar(zCombBAO3,  98.96/(1+zCombBAO3),       yerr=ersys(2.21, 1.18)/(1+zCombBAO1),
                  color='red', fmt='d', markersize=6, elinewidth=1.5)
    plt.errorbar(0.01,  69.8,       yerr=0.8,
                  color='red', fmt='d', markersize=6, elinewidth=1.5)

    plt.plot(z, [100*0.673*(0.316*(1+zi)**3 + (1-0.316))**0.5/(1+zi) for zi in z], 'red', lw=1.5)

    plt.tight_layout()
    plt.savefig('new_23LCDM_PLK_BAO.pdf')
    plt.show()
