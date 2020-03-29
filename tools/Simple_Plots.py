
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

    def __init__(self, dir_name, roots, params):
        self.dir_name = dir_name
        self.roots    = roots
        self.params   = params

        self.label    = None


    def color_legend(self, leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
            txt.set_color(line.get_color())


    def Plots1D(self, **kwargs):
        colors = kwargs.get('color', ['r', 'b'])
        plt.figure(figsize=(5*len(self.params), 4))
        for i, root in enumerate(self.roots):
            for j, param in enumerate(self.params):
                plt.subplot(1, len(self.params), j+1)
                C = cosmochain(self.dir_name + root, 'auto')

                xx, yy = C.GetHisto(param, smooth=2, NormPeak=True)
                y2 = gaussian_filter1d(yy, 2)
                plt.plot(xx, y2, label=self.label[i], color=colors[i])

                plt.xlabel(C.latexname(str(param)))
                plt.ylabel('prob.')

            leg = plt.legend(loc='upper right')
            leg.draw_frame(False)
            self.color_legend(leg)
        plt.tight_layout()
        #plt.savefig('Plot_1D.pdf')
        plt.show()