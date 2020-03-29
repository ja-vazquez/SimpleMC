
from scipy.ndimage import gaussian_filter1d
from cosmich import cosmochain
import matplotlib.pyplot as plt
from matplotlib import rcParams


rcParams.update({'backend': 'pdf',
               'axes.labelsize': 15,
               'text.fontsize': 15,
               'xtick.labelsize': 15,
               'ytick.labelsize': 15,
               'legend.fontsize': 10,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True})

class Simple_Plots(cosmochain):
    name = 'jav'

    def __init__(self, dir_name, roots):
        self.dir_name   = dir_name
        self.roots      = roots

        self.label    = None
        self.colors   = ['red', 'blue']

    def color_legend(self, leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
            txt.set_color(line.get_color())


    def Plots1D(self, params, **kwargs):
        plt.figure(figsize=(5*len(params), 4))
        for i, root in enumerate(self.roots):
            for j, param in enumerate(params):
                plt.subplot(1, len(params), j+1)
                C = cosmochain(self.dir_name + root, 'auto')

                xx, yy = C.GetHisto(param, smooth=2, NormPeak=True)
                plt.plot(xx, yy, label=self.label[i], color=self.colors[i])

                plt.xlabel(C.latexname(str(param)))
                plt.ylabel('prob.')

            leg = plt.legend(loc='upper right')
            leg.draw_frame(False)
            self.color_legend(leg)
        plt.tight_layout()
        plt.savefig('Plot_1D.pdf')
        plt.show()



    def Plots2D(self, params_pair, **kwargs):
        plt.figure(figsize=(5*len(params_pair), 4))
        for i, root in enumerate(self.roots):
            for j, params in enumerate(params_pair):
                plt.subplot(1, len(params_pair), j+1)
                C = cosmochain(self.dir_name + root)
                C.Plot2D(params[0], params[1], label=self.label[i],
                         filled=self.colors[i], solid=True)

                plt.xlabel(C.latexname(params[0]))
                plt.ylabel(C.latexname(params[1]))

            leg = plt.legend(loc='upper right')
            leg.draw_frame(False)
            self.color_legend(leg)
        plt.tight_layout()
        plt.savefig('Plot_2D.pdf')
        plt.show()





