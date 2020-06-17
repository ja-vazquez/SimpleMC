from simplemc.tools.cosmich import cosmochain
import matplotlib.pyplot as plt
from matplotlib import rcParams




class Simple_plots(cosmochain):
    name = 'jav'

    def __init__(self, dir_name, roots, label=None, **kwargs):
        self.dir_name   = dir_name
        self.roots      = roots

        self.label    = label
        self.colors   = kwargs.pop("colors", ['red', 'blue', 'green', 'orange'])
        nchains = kwargs.pop("nchains", None)
        skip = kwargs.pop("skip", 0)

        if (type(roots) == type("ch")):
            self.Clist = cosmochain(dir_name + roots)
        elif len(roots)>0:
            try:
                self.Clist = [cosmochain(dir_name + r) for r in roots]
            except:
                self.Clist = [cosmochain(dir_name + roots[0]+".txt", nchains=nchains, skip_=skip)]



    def Plots1D(self, params, **kwargs):
        smooth = kwargs.pop("smooth", 2)
        normpeak = kwargs.pop("normpeak", True)

        plt.figure(figsize=(15, 3*(len(params)//4+1)))
        for i, C in enumerate(self.Clist):
            for j, param in enumerate(params):
                plt.subplot((len(params)-1)//4+1, 4, j+1)
                xx, yy = C.GetHisto(param, smooth=smooth, NormPeak=normpeak)
                plt.plot(xx, yy, label=self.label[i], color=self.colors[i])

                plt.xlabel(C.latexname(str(param)))
                plt.ylabel('prob.')

        self.draw_frame()
        plt.tight_layout()
        #plt.savefig('Plot_1D.pdf')
        #plt.show()



    def Plots2D(self, params_pair, **kwargs):
        pbest = kwargs.pop("pbest", True)
        solid = kwargs.pop("solid", True)
        plt.figure(figsize=(5*len(params_pair), 4))
        for i, C in enumerate(self.Clist):
            for j, params in enumerate(params_pair):
                plt.subplot(1, len(params_pair), j+1)
                C.Plot2D(params[0], params[1], label=self.label[i], pbest=pbest,
                         filled=self.colors[i], solid=solid, conts=[0.68, 0.95])

                plt.xlabel(C.latexname(params[0]))
                plt.ylabel(C.latexname(params[1]))

        self.draw_frame()
        plt.tight_layout()
        #plt.savefig('Plot_2D.pdf')
        plt.show()



    def triangle(self, parlist, new_style=True, **kwargs):
        color = kwargs.pop("color", "blue")
        rcParams.update({'xtick.labelsize': 12,
                         'ytick.labelsize': 12,})
        for i, C in enumerate(self.Clist):
            C.plotAll(color=color, parlist=parlist, new_style=new_style)
            #plt.savefig('Plot_triangle.pdf')
            #plt.show()




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




    def cornerPlotter(self, params, color='g', show_titles=True, fill_contours=True):
        import corner
        for C in self.Clist:
            parcols = [C.parcol[p] for p in params]
            pchains = C.chain[:, parcols]
            lnames  = [C.latexname(l) for l in params]

            figure = corner.corner(pchains, labels=lnames, \
                               bins             = 50, \
                               weights          =C.chain[:,0], \
                               color            ='b', \
                               quantiles        =[0.5], \
                               show_titles      =True, \
                               title_fmt        ='.4f', \
                               smooth1d         =True, \
                               smooth           =True, \
                               fill_contours    =fill_contours, \
                               plot_contours    =True, \
                               plot_density     =True, \
                               levels           =(0.68,0.95),\
                               title_kwargs={"fontsize": 12})
            figure.savefig('Plot_corner.pdf')



    def fgivenx(self, params, z, func, labels=None):
        from fgivenx import plot_contours

        plist = [self.Clist.slist(p) for p in params]

        cbar1 = plot_contours(func, z, list(zip(*plist)), weights=self.Clist.chain[:, 0],
                              contour_line_levels=[1,2], linewidths = 0.8, colors=plt.get_cmap('Greens'))
        cbar1 = plt.colorbar(cbar1,ticks=[0,1,2])
        cbar1.set_ticklabels(['','$1\sigma$','$2\sigma$'])

        plt.grid()
        if labels:
            plt.ylabel(r'$%s$'%labels[1])
            plt.xlabel(r'$%s$'%labels[0])
        plt.savefig('fgivenx.pdf')
        plt.tight_layout()
        plt.show()


