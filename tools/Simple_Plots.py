
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

        self.Clist = [cosmochain(dir_name + r) for r in roots]



    def Plots1D(self, params, **kwargs):
        plt.figure(figsize=(5*len(params), 4))
        for i, C in enumerate(self.Clist):
            for j, param in enumerate(params):
                plt.subplot(1, len(params), j+1)

                xx, yy = C.GetHisto(param, smooth=2, NormPeak=True)
                plt.plot(xx, yy, label=self.label[i], color=self.colors[i])

                plt.xlabel(C.latexname(str(param)))
                plt.ylabel('prob.')

        self.draw_frame()
        plt.tight_layout()
        plt.savefig('Plot_1D.pdf')
        plt.show()



    def Plots2D(self, params_pair, **kwargs):
        plt.figure(figsize=(5*len(params_pair), 4))
        for i, C in enumerate(self.Clist):
            for j, params in enumerate(params_pair):
                plt.subplot(1, len(params_pair), j+1)
                C.Plot2D(params[0], params[1], label=self.label[i], pbest=True,
                         filled=self.colors[i], solid=True, conts=[0.68, 0.95])

                plt.xlabel(C.latexname(params[0]))
                plt.ylabel(C.latexname(params[1]))

        self.draw_frame()
        plt.tight_layout()
        plt.savefig('Plot_2D.pdf')
        plt.show()




    def plotAlls(self, parlist, new_style=True):
        rcParams.update({'xtick.labelsize': 12,
                         'ytick.labelsize': 12,})
        for i, C in enumerate(self.Clist):
            C.plotAll(color='blue', parlist=parlist, new_style=new_style)
            plt.savefig('Plot_triangle.pdf')
            plt.show()




    def Show_limits(self, params):
        for i, C in enumerate(self.Clist):
            for j, param in enumerate(params):
                x = C.GetLimits(param, returnlims=True, ML=True)
                print ('--'*10)
                print (param, '[3,2,1] sigma, best-fit', x)



    def Covariance(self, params):
        for i, C in enumerate(self.Clist):
            x = C.GetCovariance(params)
            print ('--'*10)
            print ('mean', x[0], '\n cov', x[1])
            return x



    def color_legend(self, leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
            txt.set_color(line.get_color())


    def draw_frame(self):
            leg = plt.legend(loc='upper right')
            leg.draw_frame(False)
            self.color_legend(leg)