from scipy import *
from cosmich import *
import pylab
import sys


class plot_aux():

    def colour(self, x):
        if x == 0:
            return 'red'
        if x == 1:
            return 'blue'
        if x == 2:
            return 'green'
        if x == 3:
            return 'black'
        if x == 4:
            return 'orange'
        if x == 5:
            return 'magenta'
        if x == 6:
            return 'cyan'
        if x == 7:
            return 'yellow'
        if x > 9:
            print("Increase colouring")

    def label_m(self, file, name):
        with open(file+'.paramnames') as inf:
            for line in inf:
                parts = line.split()
                if len(parts) > 1:
                    if name in parts:
                        return '$'+str(parts[1])+'$'

    def ldata(self, dataset):
        cosmodata = ''
        if 'BBAO+SN' == dataset:
            cosmodata = 'BAO+SN'
        elif 'BBAO+Planck' == dataset:
            cosmodata = 'BAO+Planck'
        elif 'SN+Planck' == dataset:
            cosmodata = 'SN+Planck'
        elif 'BBAO+SN+Planck' == dataset or 'BBAO+Planck+SN' == dataset:
            cosmodata = 'BAO+SN+Planck'
        elif 'GBAO' == dataset:
            cosmodata = 'Galaxy BAO'
        elif 'LBAO' == dataset:
            cosmodata = 'Lyman-$\\alpha$ BAO'
        elif 'GBAO+PlDa' == dataset:
            cosmodata = 'Galaxy BAO+Planck $D_M$'
        elif 'LBAO+PlDa' == dataset:
            cosmodata = 'Lyman-$\\alpha$ BAO+Planck $D_M$'
        elif 'LBAO+Planck' == dataset:
            cosmodata = 'Lyman-$\\alpha$ BAO+Planck'
        elif 'BBAO' == dataset:
            cosmodata = 'Combined BAO'
        elif 'Planck' == dataset:
            cosmodata = 'Planck'
        elif 'BBAO+PlDa' == dataset:
            cosmodata = 'Combined BAO+Planck $D_M$'
        elif 'BBAO+PlRd' == dataset:
            cosmodata = 'Combined BAO+Planck $r_d$'
        elif 'BBAO+SN+PlRd' == dataset:
            cosmodata = 'BAO+SN+Planck $r_d$'
        elif 'PLA:base_mnu_Alens:planck_lowl_lowLike:post_BAO' == dataset:
            cosmodata = 'PLA Planck+BAO ($A_L$ free)'
        elif 'PLA:base_mnu:planck_lowl_lowLike:post_BAO' == dataset:
            cosmodata = 'PLA Planck+BAO ($A_L=1$)'
        elif 'PLA:base:planck_lowl_lowLike:' == dataset:
            cosmodata = 'Planck (full)'
        elif 'mnu_AL_PL_BBAO' == dataset:
            cosmodata = 'BAO+Planck (full; $A_L$ free)'
        elif 'mnu_PL_BBAO' == dataset:
            cosmodata = 'BAO+Planck (full; $A_L=1$)'
        elif 'NFC:r_PL_LBAO' == dataset:
            cosmodata = 'Lyman-$\\alpha$ BAO + Planck (full, $N_{\\rm eff},r$)'
        elif 'NFC:r_PL_BBAO_JLA' == dataset:
            cosmodata = 'Lyman-$\\alpha$ BAO + Planck (full, $N_{\\rm eff}$)'
        elif 'NFC:r_PL_BBAO' == dataset:
            cosmodata = 'BAO + Planck (full, $N_{\\rm eff},r$)'
        else:
            print('Add dataset name')
        return cosmodata

    def get_filename(self, dire, model, extra, dataset):
        if 'PLA' in dataset:
            tmp, model, data, post = dataset.split(':')
            if len(post) > 0:
                post = "_"+post
            return '/astro/u/anze/work/Planck/PLA/%s/%s/%s' % (model, data, model+"_"+data+post)
        elif 'mnu' in dataset:
            return 'chains_CosmoMC/%s' % (dataset)
        elif 'NFC' in dataset:
            return 'chains_CosmoMC/'+model+'_'+dataset.split(':')[1]
        else:
            return dire + model+'_'+extra+'_'+dataset

    def Plotting_1d(self, dire, modell, extra, datasetl, param, Nb=70, Ol=False, loc=None, legcolor=False, lwidth=None, linecolor=None, legend=None, minmax=None):
        pylab.figure()
        a, b, cc = 0, 0, 0

        for model in modell:
            b += 1
            for dataset in datasetl:
                # a+=1
                fname = self.get_filename(dire, model, extra, dataset)
                C = cosmochain(fname)

                if Ol == 'True':
                    if model == 'oLCDM':
                        C['Ol'] = 1-C['Ok']-C['Om']
                    else:
                        C['Ol'] = 1-C['Om']
                if "Decay" in model:
                    if "01" in model:
                        fx = 0.1
                    elif "05" in model:
                        fx = 0.5
                    else:
                        fx = 1.0
                    C['lambdaf'] = C['lambda']*fx
                    C.lname['lambdaf'] = '\\lambda f_x'
                if type(Nb) == type([]):
                    Nbx = Nb[cc]
                else:
                    Nbx = Nb
                if minmax:
                    xx, yy = C.GetHisto(
                        param[0], NormPeak=False, nbins=Nbx, mnval=minmax[0], mxval=minmax[1])
                else:
                    xx, yy = C.GetHisto(param[0], NormPeak=False, nbins=Nbx)

                if legend:
                    label = legend[cc]
                else:
                    label = self.ldata(dataset) if b == 1 else ""
                # if minmax:
                xx, yy = smline(xx, yy, minmax)
                #  pass
                if linecolor:
                    line = linecolor[cc]
                else:
                    line = self.colour(a)

                if lwidth:
                    lwh = lwidth[cc]
                else:
                    lwh = 4

                pylab.plot(xx, yy, line, lw=lwh, label=label)
                a += 1
                cc += 1
        if Ol == 'True':
            pylab.xlabel('$\\Omega_{\Lambda}$')
        else:
            pylab.xlabel(C.latexname(param[0]))
        pylab.ylabel('$\\rm prob.$')
        if loc:
            leg = pylab.legend(loc=loc)
        else:
            leg = pylab.legend()

            # No box in the legend, and colour legend
        if legcolor:
            leg.draw_frame(False)
            color_legend(leg)

    def pname_translate(self, dataset, p):
        if "PLA" in dataset or "NFC" in dataset:
            if p == 'Nnu':
                p = 'nnu'

        return p

    def Plotting_2d(self, dire, modell, extra, datasetl, param, Nbin=50, Ol=False,
                    lwidth=None, loc=None, fillcolor=None, legcolor=False, labels=None):
        pylab.figure(figsize=(9, 6))
        a = 0
        for model in modell:
            for dataset in datasetl:
                if type(Nbin) == type([1, 2, 3]):
                    Nb = Nbin[a]
                else:
                    Nb = Nbin
                C = cosmochain(self.get_filename(dire, model, extra, dataset))
                if 'NFC' in dataset:
                    Nb = 19
                try:
                    C['h'] = C['H0*']/100.0
                    C.lname['h'] = 'h'
                except:
                    pass
                try:
                    C['Nnu'] = C['nnu']*1.0
                    C.lname['Nnu'] = 'N_{\\rm eff}'
                except:
                    pass
                try:
                    C['Om'] = C['omegam*']*1.0
                    C.lname['Om'] = '\Omega_m'
                except:
                    pass
                if "Decay" in model:
                    if "01" in model:
                        fx = 0.1
                    elif "05" in model:
                        fx = 0.5
                    else:
                        fx = 1.0
                    C['lambdaf'] = C['lambda']*fx
                    C.lname['lambdaf'] = '\\lambda f_x'

                if lwidth:
                    lwh = lwidth[a]
                else:
                    lwh = 3

                if fillcolor:
                    fcol = fillcolor[a]
                    print(a, fillcolor[a])
                else:
                    fcol = self.colour(a)

                if Ol == 'True':
                    C['Ol'] = 1-C['Ok']-C['Om']
                    C.lname['Ol'] = '\Omega_\Lambda'
                if labels:
                    lab = labels[a]
                else:
                    lab = self.ldata(dataset)
                C.Plot2D(param[0], param[1], filled=fcol,
                         lw=lwh, label=lab, N=Nb)
                a += 1

        file = str(dire + model+'_'+extra+'_'+dataset)
        pylab.xlabel(C.latexname(param[0]))
        pylab.ylabel(C.latexname(param[1]))
        if loc:
            leg = pylab.legend(loc=loc)
        else:
            leg = pylab.legend()
            # No box in the legend, and colour legend
        if legcolor:
            leg.draw_frame(False)
            color_legend(leg)

    def Plotting_waDE_2d(self, dire, modell, extra, datasetl, param, lwidth=None, Nb=40, loc=None, legcolor=False):
        pylab.figure(figsize=(9, 6))
        a = 0
        for model in modell:
            for dataset in datasetl:
                C = cosmochain(dire + model+'_'+extra+'_'+dataset)
                C['w_z'] = C['w']+0.21*C['wa']
                if lwidth:
                    lwh = lwidth[a]
                else:
                    lwh = 3
                # print C.GetLimits('w_z')
                C.Plot2D(param[0], param[1], filled=self.colour(
                    a), lw=lwh, label=self.ldata(dataset), N=Nb)
                a += 1
        file = str(dire + model+'_'+extra+'_'+dataset)
        pylab.xlabel('$w_{(z=zp)}$')
        pylab.ylabel(self.label_m(file, param[1]))
        if loc:
            leg = pylab.legend(loc=loc)
        else:
            leg = pylab.legend()

            # No box in the legend, and colour legend
        if legcolor:
            leg.draw_frame(False)
            color_legend(leg)

    def plotH(sefl, y, dire, chainname, color, txt):
        if "PLA" in chainname:
            C = cosmochain(chainname, 'auto')
            C['h'] = C['H0*']  # /100
            m, p1, p2, p3, m1, m2, m3 = C.GetLimits('h')
        elif "N:" in chainname:
            m, e = list(map(float, chainname.split()[1:]))
            p1, p2, p3 = e, 2*e, 3*e
            m1, m2, m3 = p1, p2, p3
        else:
            C = cosmochain(dire+chainname)
            C['h'] = C['h']*100
            m, p1, p2, p3, m1, m2, m3 = C.GetLimits('h')
        pylab.errorbar(m, y, xerr=[[m1], [p1]], fmt='-o' +
                       color, lw=4, capthick=4, capsize=10, ms=10)
        print(txt, (p1+m1)/(2*m), m, p1, p2, p3, m1, m2, m3)
        pylab.text(m+p1+0.5, y-0.05, txt, verticalalignment='baseline')

    def Plotting_EDE_2d(self, dire, modell, extra, datasetl, param, Nb=50, Ol=False, loc=None, legcolor=False):
        pylab.figure(figsize=(9, 6))
        a = 0
        for model in modell:
            for dataset in datasetl:
                # a+=1
                C = cosmochain(self.get_filename(dire, model, extra, dataset))
                a += 1
                C.Plot2D(param[0], param[1], filled='black' if a ==
                         1 else 'blue', label=self.ldata(dataset), N=Nb)

        pylab.xlabel(C.latexname(param[0]))
        pylab.ylabel(C.latexname(param[1]))
        if loc:
            leg = pylab.legend(loc=loc)
        else:
            leg = pylab.legend()
            # No box in the legend, and colour legend
        if legcolor:
            leg.draw_frame(False)
            color_legend(leg)

        Ode = arange(0, 0.6, 0.01)
        Om = [0.302/(1+Odes) for Odes in Ode]
        Om2 = [0.302*(1-Odes) for Odes in Ode]

        # pylab.plot(Ode,Om,'r--')
        pylab.plot(Ode, Om2, 'r-')


def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
    for line, txt in zip(leg.get_lines(), leg.get_texts()):
        txt.set_color(line.get_color())
