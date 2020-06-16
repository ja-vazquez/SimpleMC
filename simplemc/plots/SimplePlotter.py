#!/usr/bin/env python
import sys
from simplemc.plots.Simple_Plots import Simple_plots
import matplotlib.pyplot as plt
import numpy as np
import webbrowser

class SimplePlotter:
    def __init__(self, chainsdir, listpars, path=None, root=None, show=False, weights=None):
        self.chainsdir = chainsdir
        self.filename = path
        self.listpars = listpars
        if root == None:
            self.root = path.replace("{}/".format(chainsdir), "")
        else:
            self.root = root
        self.ndim = len(listpars)
        self.show = show
        self.weights = weights



    def simpleGetdist(self, **kwargs):
        """
        Lewis (2019)
        arXiv:1910.13970v1 [astro-ph.IM]
        """
        from getdist import plots, MCSamples, chains
        smooth2d = kwargs.pop("smooth2d", 0.3)
        smooth1d = kwargs.pop("smooth1d", 0.3)
        burnin = kwargs.pop("burnin", 0.2)
        colors = kwargs.pop("colors", ['red', 'blue', 'black', 'green', 'yellow', 'purple', 'gray'])
        legend_labels = kwargs.pop("legend_labels", [])
        filled = kwargs.pop("filled", False)
        normalized = kwargs.pop("normalized", False)
        shaded = kwargs.pop("shaded", False)
        label = kwargs.pop("label", None)
        roots = kwargs.pop('roots', [self.root])

        g = plots.getSinglePlotter(chain_dir=self.chainsdir, width_inch=10,
                                   ratio=0.9,
                                   analysis_settings={'smooth_scale_2D': smooth2d,
                                                      'smooth_scale_1D': smooth1d,
                                                      'ignore_rows': burnin})
        g.settings.lab_fontsize = 14
        g.settings.legend_fontsize = 12
        g.settings.axes_fontsize = 12
        g.triangle_plot(roots, self.listpars,
                        diag1d_kwargs={'colors':colors},
                        colors=colors,
                        legend_labels=legend_labels,
                        filled=filled,
                        normalized=normalized, shaded=shaded)

        self.image = "{}_getdist.png".format(self.filename)
        self.saveFig(label)


    def simpleCorner(self, **kwargs):
        """
            Corner method.
            Daniel Foreman-Mackey (2016)
            doi = {10.21105/joss.00024},
            url = {https://doi.org/10.21105/joss.00024},
        """
        import corner
        color = kwargs.pop("color", 'g')
        show_titles = kwargs.pop("show_titles", True)
        fill_contours = kwargs.pop("fill_contours", True)
        bins = kwargs.pop("bins", 20)
        smooth1d = kwargs.pop("smooth1d" ,True)
        smooth = kwargs.pop("smooth" ,True)
        plot_contours = kwargs.pop("plot_contours",True)
        plot_density =kwargs.pop("plot_density", True)
        truths = kwargs.pop("truths", None)
        label = kwargs.pop("label", None)

        print("Plotting with Corner!")
        self.readFile()
        figure = corner.corner(self.samples, labels=self.latexnames[0:self.ndim],
                               bins=bins, weights=self.weights, color=color,
                               quantiles=[0.5], show_titles=show_titles, title_fmt='.4f',
                               smooth1d=smooth1d, smooth=smooth, fill_contours=fill_contours,
                               plot_contours=plot_contours, plot_density=plot_density,
                               truths=truths, truth_color='#4682b4', title_kwargs={"fontsize": 12})
        self.image = "{}_corner.png".format(self.filename)
        self.saveFig(label)

    def simpleFgivenx(self, **kwargs):
        """
        Handley, (2018). https://doi.org/10.21105/joss.00849

        """
        from fgivenx import plot_contours, samples_from_getdist_chains

        params = kwargs.pop("params", ['Om', 'h'])
        z = kwargs.pop("z", np.linspace(0, 4, 100))
        func = kwargs.pop("func", self.Hzfunc)
        labels = kwargs.pop("labels", ['z', 'H(z)'])
        interval = kwargs.pop("interval", [0, 4])
        colors = kwargs.pop("colors", plt.get_cmap('Greens'))
        linewidths = kwargs.pop("line_widths", 0.8)
        file_root = self.filename
        try:
            samples, weights = samples_from_getdist_chains(params, file_root)
        except:
            samples, weights = samples_from_getdist_chains(params, file_root+'_')
        # z = np.logspace(-4, 1, 100)
        z = np.linspace(interval[0], interval[1], 100)
        cbar = plot_contours(func, z, samples, weights=weights,
                             contour_line_levels=[1,2], linewidths=linewidths,
                             colors=colors)
        cbar = plt.colorbar(cbar, ticks=[0, 1, 2, 3])
        #cbar.set_ticklabels(['', r'$1\sigma$', r'$2\sigma$', r'$3\sigma$'])
        cbar.set_ticklabels(['', '$1\sigma$', '$2\sigma$'])
        # plt.xscale('log')
        plt.tight_layout()
        plt.grid()
        if labels:
            plt.ylabel(r'$%s$' % labels[1])
            plt.xlabel(r'$%s$' % labels[0])

        self.image = "{}_fgivenx.png".format(self.filename)
        self.saveFig()

    def simplePlot(self, **kwargs):
        """
        Native simplemc plotter
        """
        from simplemc.tools.Simple_Plots import Simple_plots
        type = kwargs.pop('type', 'triangle')
        roots = kwargs.pop('roots', [self.root])
        nchains = kwargs.pop('nchains', None)

        # parlist is a par for 2d
        # list of paramnames for triangle
        # one or more for 1d
        label = kwargs.pop("label", [""])
        colors = kwargs.pop("colors", ['red', 'blue', 'green', 'orange'])
        #1d
        pars1d = kwargs.pop('pars1d', ['Om', 'h'])
        smooth1d = kwargs.pop("smooth1d", 2)
        normpeak1d = kwargs.pop("normpeak1d", True)
        #2d
        pars2d = kwargs.pop('pars2d', [['Om', 'h']])
        pbest2d = kwargs.pop("pbest2d", True)
        solid2d = kwargs.pop("solid2d", True)
        #triangle
        parstriangle = kwargs.pop('parstriangle', ['Om', 'Obh2', 'h'])
        colortriangle = kwargs.pop("colortriangle", "blue")
        fig = Simple_plots(self.chainsdir+"/", roots, label=label,
                           colors=colors, nchains=nchains)
        if type == "triangle" or type == "tri":
            fig.triangle(parstriangle, color=colortriangle)
        elif type == "1D" or type == "1d":
            fig.Plots1D(pars1d, smooth=smooth1d, normpeak=normpeak1d)
        elif type == "2D" or type == "2d":
            fig.Plots2D(pars2d, pbest=pbest2d, solid=solid2d)
        else:
            sys.exit("Invalid option")
        self.image = "{}_{}_simple.png".format(self.filename, type)
        self.saveFig()

    def simplex_vs_y(self, **kwargs):
        """
        Plot 2 columns given a text file
        """
        file = self.filename + ".txt"
        usecols = kwargs.pop("usecols", (0, 1))
        xlabel = kwargs.pop("xlabel", "x")
        ylabel = kwargs.pop("ylabel", "y")
        color = kwargs.pop("color", "g")
        linestyle = kwargs.pop("linestyle", "-")
        linewidth = kwargs.pop("linewidth", 4)

        data = np.loadtxt(file, usecols=usecols)
        x = data[:, 0]
        y = data[:, 1]
        plt.plot(x, y, color=color, linestyle=linestyle,
                 linewidth=linewidth)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        self.image = "{}_{}_vs_{}.png".format(self.filename, xlabel, ylabel)
        self.saveFig()



    def readFile(self):
        """
         This method reads the samples and the .param file.
        """
        labelsfile = open(self.filename + ".paramnames", 'r')
        self.latexnames = []
        self.paramnames = []

        for item in labelsfile:
            self.latexnames.append('$' + item.split('\t\t\t')[1].strip('\n') + '$')
            self.paramnames.append(item.split('\t\t\t')[0].strip('\n'))

        labelsfile.close()
        try:
            npchain = np.loadtxt(self.filename+'.txt')
        except:
            npchain = np.loadtxt(self.filename + '_1.txt')

        self.samples = npchain[:, 2:self.ndim + 2]



    def saveFig(self, label=None):
        if label is not None:
            plt.text(0.6, 2.0, label, transform=plt.gca().transAxes)
        plt.savefig(self.image, bbox_inches='tight')
        if self.show:
            webbrowser.open(self.image)

    def Hzfunc(self, z, theta1):
        Omega_m, h = theta1
        Hz = 100 * h * (Omega_m * (1 + z) ** 3 + (1 - Omega_m)) ** 0.5
        return Hz

