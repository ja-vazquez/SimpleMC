#!/usr/bin/env python
#
# Finally a getdist improvement
#
#

import sys
import scipy as sp
from glob import glob
import numpy.fft as fft
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage.filters import gaussian_filter


# chains are col0: weight, col1:Likel, rest:params
def myloadtxt(fname, cols=None):
    # fixed the annoying "setting an array element with a sequence."
    # when chain still runnig.
    da  = []
    bad = []
    for line in open(fname).readlines()[:-2]:
        try:
            cline = list(map(float, line.split()))
        except ValueError:
            continue
        if (cols == None):
            cols = len(cline)
        if len(cline) == cols:
            da.append(cline)
        else:
            bad.append(cline)

    da = sp.array(da)
    if (len(bad)):
        sys.exit("exit: BAD= format")
    return da



class cosmochain:
    def __init__(self, root, nchains='auto', skip_=None, temp=1.0, \
                       balance=True, weightfunc=None, kirkby=False):
        """
        Class to analize chains.
        Parameters
        ----------
        root
        nchains
        skip_
        temp
        balance
        weightfunc
        kirkby

        Returns
        -------

        """

        if (balance and ("PLA" in root or "neffch" in root or "Decay" in root)):
            balance = False

        # get list of chain_files
        if type(nchains) == type('auto'):
            if nchains == 'auto':
                flist = glob(root + '_?.txt') + glob(root + '_??.txt')
        else:
            if (nchains == None):
                flist = [root]
            elif (type(nchains[0]) == type(1)):
                flist = [root + '_' + str(x) + '.txt' for x in nchains]

        if len(flist) == 0:
            sys.exit("exit: Bad chain spec." + root +'_?.txt')

        # get parameter names
        try:
            lines = open(flist[0].replace('.txt', '.paramnames')).readlines()
        except:
            try:
                lines = open(root + '.paramnames').readlines()
            except:
                try:
                    lines = open(root + '.params').readlines()
                except:
                    sys.exit("exit: neither params nor paramnames included")

        self.paramnames = [line.split()[0] for line in lines]
        self.latexnames = [' '.join(line.split()[1:]) for line in lines]
        self.lname      = dict(zip(self.paramnames, self.latexnames))

        self.parcol = {}
        for i, name in enumerate(self.paramnames):
            print(i, name)
            self.parcol[name] = i+2
        print("Got ", len(self.paramnames), "parameters.")


        data = []
        for fname in flist:
            print("Reading", fname, "...", end=' ')
            if (kirkby):
                skip = 3
                print("kirkby style ", end=' ')
                cdata = open(fname).readlines()[skip:-1]
                cdata = sp.array([[1, 0] + list(map(float, x.split()))
                               for x in cdata])
            else:
                da = myloadtxt(fname, cols=len(self.paramnames)+2)
                if len(da) == 0:
                    continue
                if (skip_ == None):
                    finlike = da[-1, 1]
                    ii = 0
                    while (da[ii, 1] > finlike):
                        ii += 1
                    # plus 10
                    skip  = ii+20
                    cdata = da[skip:-1]
                else:
                    skip = skip_
                    cdata = da[skip:-1]
                print("skipping:", skip)
                if (balance):
                    print("balancing...")
                    # make weights the same in average to avoid issues with
                    # temperature
                    cdata[:, 0] /= cdata[:, 0].mean()
                cdata = da[skip:-1]

            data += list(cdata)
        self.chain = sp.array(data)

        if weightfunc != None:
            print("Reweighting")
            for i, line in enumerate(self.chain):
                self.chain[i, 0] = weightfunc(line)

        #print(len(self.chain), len(self.chain[0]))
        del data

        try:
            self.bestarg = self.chain[:, 1].argmin()
            self.best    = self.chain[self.bestarg]
        except:
            print("WTF?")

        if (temp != 1):
            like = self.chain[:, 1]
            like = like-like.min()
            self.chain[:, 0] *= exp(-(temp-1.0)*like)




    def latexname(self, name):
        #print(self.lname[name])
        return '$'+ self.lname[name] +'$'



    #get column of a given parameter
    def slist(self, name):
        return self.chain[:, self.parcol[name]]



    def __getitem__(self, key):
        if (type(key) == type("p")):
            key = self.parcol[key]
        return self.chain[:, key]



    def __setitem__(self, key, res):
        if (type(key) == type("p")):
            if key in self.parcol:
                key = self.parcol[key]
            else:
                # need to add a column
                N = len(self.chain)
                self.chain = sp.concatenate((self.chain, sp.zeros((N, 1))), 1)
                Nc = len(self.chain[0])-1
                self.parcol[key] = Nc
                self.paramnames.append(key)
                key = Nc
        self.chain[:, key] = res




    def BestSample(self):
        ii = self.chain[:, 1].argmin()
        return self.chain[ii, :]



    def GetMarginal(self, param, val):
        param = self.parcol[param]
        return ((self.chain[:, param] > val)*self.chain[:, 0]).sum()/self.chain[:, 0].sum()



    def GetHisto(self, param, nbins=50, ncolumn=None, minval=None, maxval=None, \
                       smooth=None, NormPeak=False, plot=None, lw=4):
        "Returns histogram for plotting."
        if ncolumn != None:
            column = ncolumn
        else:
            if (type(param) == type("st")):
                param = self.parcol[param]
            column = self.chain[:, param]

        if minval  == None:    minval = column.min()
        if maxval  == None:    maxval = column.max()  # (to add the last one)
        if (minval == maxval): return None, None
        maxval *= 1.001

        step = (1.0*maxval-1.0*minval)/nbins
        tmp  = list(map(int, (column - minval)/step))

        histo = sp.zeros((nbins,))

        for ii in range(len(tmp)):
            if (tmp[ii] >= 0) and (tmp[ii] < nbins):
                histo[tmp[ii]] += self.chain[ii, 0]

        xval = sp.array([minval+(x+0.5)*step for x in range(nbins)])
        yval = histo/step
        #print(xval, minval, maxval, step)

        if smooth:
            yvalpad = sp.array([0, 0, 0, 0]+list(yval)+[0, 0, 0, 0])
            if smooth == 1:
                yval = (yvalpad[3:nbins+3] + yvalpad[4:nbins+4] +
                        yvalpad[5:nbins+5])/3.0
            if smooth == 2:
                yval = (yvalpad[2:nbins+2] + yvalpad[3:nbins+3] + yvalpad[4:nbins+4]+
                        yvalpad[5:nbins+5]++yvalpad[6:nbins+6])/5.0

        if (NormPeak):
            yval /= yval.max()
        else:
            area = yval.sum()*step
            yval /= area

        yval = gaussian_filter1d(yval, sigma=2)
        if (plot != None): plt.plot(xval, yval, plot, lw=lw)
        return xval, yval*1.01



    def Plot1D(self, p1, sty='r-', label="", N=50):
        xx, yy = self.GetHisto(p1, nbins=N)
        plt.plot(xx, yy, sty, label="", lw=2)



    def Plot2D(self, param1, param2, N=60, lims=None, conts=[0.68, 0.95, 0.997], filled=True, lw=2,\
                     pbest=False, blur=None, ncolumn1=None, ncolumn2=None, label="", solid=True):

        if (ncolumn1 == None):
            if (type(param1) == type("p")):
                param1 = self.parcol[param1]
            xx = self.chain[:, param1]
        else:
            xx = ncolumn1

        if (ncolumn2 == None):
            if (type(param2) == type("p")):
                param2 = self.parcol[param2]
            yy = self.chain[:, param2]
        else:
            yy = ncolumn2

        we = self.chain[:, 0]

        if (lims == None):
            xmin, xmax = xx.min(), xx.max()
            de = (xmax-xmin)/100.
            if (de == 0):
                de = xmax/100.
            xmin -= de
            xmax += de

            ymin, ymax = yy.min(), yy.max()
            de = (ymax-ymin)/100.
            if (de == 0):
                de = ymax/100.
            ymin -= de
            ymax += de
        else:
            xmin, xmax, ymin, ymax = lims

        out   = 0
        grid  = sp.zeros((N, N))
        for x, y, w in zip(xx, yy, we):
            i1 = int((x-xmin)/(xmax-xmin)*N)
            i2 = int((y-ymin)/(ymax-ymin)*N)
            try:
                grid[i2, i1] += w
            except:
                out += w

        if (out > 0):
            print("warn: out =", out/we.sum())

        b = grid.flatten()
        b = b.tolist()
        b.sort(reverse=True)
        b = sp.array(b)

        c = b*1.0
        c = c.cumsum()
        c /= c[-1]

        l1, l2, l3 = 0, 0, 0

        print (conts, len(conts) >= 3)
        for val, cum in zip(b, c):
            if (cum > conts[0]) and (l1 == 0):
                l1 = val
            if (cum > conts[1]) and (l2 == 0):
                l2 = val
            if (len(conts) >= 3) and (cum > conts[2]) and (l3 == 0):
                l3 = val

        print('contours', l1, l2, l3)
        #grid =smline(grid)

        limits = ( xmin, xmax, ymin, ymax)
        lcont  = [l3, l2, l1] if (len(conts) >= 3) else [l2, l1]
        grid   = gaussian_filter(grid, sigma= 0.9)

        mcolor ={'blue':'Blues', 'red':'Reds'}

        if type(filled) == type('string'):
            print('lw=', lw)
            plt.contour(grid, lcont, extent=limits, origin='lower', \
                            aspect='auto', colors=filled, linewidths=lw)
            if (label != ""):
                plt.plot([], [], color=filled, linewidth=lw, label=label)

            if solid:
                plt.contourf(grid, lcont + [max(b)], extent=limits, origin='lower', \
                                aspect='auto', cmap=plt.get_cmap(mcolor[filled]))

        if (pbest):
            plt.plot(self.best[param1], self.best[param2], 'bo')

        return grid, lcont+[max(b)], limits




    def GetLimits(self, param, ML=False, nch=None,
                  limlist=[0.5-0.997/2, 0.5-0.95/2, 0.5-0.68/2, 0.5, 0.5+0.68/2, 0.5+0.95/2, 0.5+0.997/2],
                  returnlims=False):
        if nch != None:
            lis = list(zip(nch, self.chain[:, 0]))
        else:
            if (type(param) == type("st")):
                param = self.parcol[param]
            lis = list(zip(self.chain[:, param], self.chain[:, 0]))

        lis.sort()
        lis  = sp.array(lis)
        pars = sp.array(lis[:, 0])
        wei  = sp.array(lis[:, 1])
        wei  = wei.cumsum()
        wei  = wei/wei[-1]
        lims = []
        for lim in limlist:
            a = sp.where(wei > lim)[0][0]
            lims.append(pars[a])
        if (ML):
            print("Central moved:", lims[3], end=' ')
            if nch != None:
                lims[3] = nch[self.bestarg]
            else:
                lims[3] = self.best[param]
            print("->", lims[3])

        if returnlims:
            return lims
        m, ep1, ep2, ep3, em1, em2, em3 = lims[3], lims[4]-lims[3], lims[5] - \
            lims[3], lims[6]-lims[3], lims[3] - \
            lims[2], lims[3]-lims[1], lims[3]-lims[0]
        return (m, ep1, ep2, ep3, em1, em2, em3)



    def GetCovariance(self, parlist):
        N   = len(parlist)
        cov = sp.zeros((N, N))
        mn  = sp.zeros(N)
        sw  = 0.0
        for el in self.chain:  # [:100]:
            sw  += el[0]
            vec  = sp.array([el[self.parcol[i]] for i in parlist])
            mn  += vec*el[0]
            cov += el[0]*sp.array([[v1*v2 for v1 in vec] for v2 in vec])
            # print cov[0,0], vec[0],sw
        mn  /= sw
        cov /= sw

        cov -= sp.array([[v1*v2 for v1 in mn] for v2 in mn])
        return mn, cov



    def plotAll(self, color, justlines=False, parlist=None, new_style=False):
        cc = 0
        if not parlist:
            N = len(self.paramnames)
            parlist = list(range(N))
        else:
            N = len(parlist)
            for i, el in enumerate(parlist):
                if type(el) == type('string'):
                    parlist[i] = self.parcol[el]-2


        if new_style:
            fig, axs = plt.subplots(N, N, sharex='col',
                                    gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(12,10))
            for ic, i in enumerate(parlist):
                for jc, j in enumerate(parlist):
                    cc += 1
                    if (ic < jc):
                        axs[ic, jc].axis('off')
                        continue
                    if (i < 0) or (j < 0):
                        continue

                    if jc>0: axs[ic, jc].tick_params(labelleft=False)
                    if (ic == jc):
                        xv, yv = self.GetHisto(i+2, smooth=2, NormPeak=True)
                        axs[ic,jc].plot(xv, yv, '-', color=color)
                    elif (ic > jc):
                        smt =self.Plot2D(j+2, i+2, filled=1, conts=[0.68, 0.95])
                        axs[ic, jc].contour(smt[0], smt[1], extent=smt[2], origin='lower', \
                                    aspect='auto', colors='blue', linewidths=1)
                        axs[ic, jc].contourf(smt[0], smt[1], extent=smt[2], origin='lower', \
                                    aspect='auto', cmap=plt.get_cmap('Blues'))

                    if (jc == 0 and ic >0):
                        if not justlines:
                            axs[ic, jc].set_ylabel(self.latexname(self.paramnames[i]))

                    if (ic == N-1) and (not justlines):
                        axs[ic, jc].set_xlabel(self.latexname(self.paramnames[j]), fontsize=10)
        else:
            for ic, i in enumerate(parlist):
                for jc, j in enumerate(parlist):
                    cc += 1
                    if (ic < jc):
                        continue
                    if (i < 0) or (j < 0):
                        continue

                    plt.subplot(N, N, cc)
                    if (ic == jc):
                        xv, yv = self.GetHisto(i+2, smooth=2, NormPeak=True)
                        plt.plot(xv, yv, '-', color=color)
                    elif (ic > jc):
                        print(i, j, 'aaa', N)
                        self.Plot2D(j+2, i+2, filled=color, conts=[0.68, 0.95])
                        if (jc == 0):
                            if not justlines:
                                plt.ylabel(self.latexname(self.paramnames[i]))

                    if (ic == N-1) and (not justlines):
                        plt.xlabel(self.latexname(self.paramnames[j]), fontsize=10)




    def GetLimitsOld(self, param, lims=[0.6826894920, 0.9544997360, 0.9973002039]):
        "returns median and pairs corresponding to lims"
        line = [(x[param+1], x[0]) for x in self.chain]
        # print line
        line.sort()
        sw = self.chain[:, 0].sum()
        lowl = [None]*len(lims)
        highl = [None]*len(lims)
        lowtrig = 0.5- sp.array(lims)/2
        hightrig = 0.5+ sp.array(lims)/2

        # print sw, hightrig

        sum = 0
        mean = None
        for qq in line:
            sum += qq[1]
            if (not mean) and (sum >= 0.5*sw):
                mean = qq[0]
            for ii in range(len(lims)):
                if (not lowl[ii]) and (sum >= lowtrig[ii]*sw):
                    lowl[ii] = qq[0]
                if (not highl[ii]) and (sum >= hightrig[ii]*sw):
                    highl[ii] = qq[0]
        lowl -= mean
        highl -= mean
        str = "%2.0f^{+%2.0f+%2.0f+%2.0f}_{%2.0f%2.0f%2.0f}" % (
            mean, highl[0], highl[1], highl[2], lowl[0], lowl[1], lowl[2])
        print(str)

        return mean, lowl, highl




def smline(x, y, mnmx):
    # first lets.pad
    y = sp.log(y+1e-30)
    N = len(y)
    y = sp.array([y[0]]*N+list(y)+[y[-1]]*N)
    rft = fft.rfft(y)
    Nx = len(rft)
    k = sp.linspace(0, 1, Nx)
    rft *= sp.exp(-k*k/(2*0.2**2))
    y = fft.irfft(rft)
    y = y[N:2*N]
    y = sp.exp(y)
    return x, y


def smline2(pic):
    Nx = len(pic[0])
    Ny = len(pic)
    picp = sp.zeros((3*Nx, 3*Ny))
    picp[Nx:2*Nx, Ny:2*Ny] = sp.log(pic+1e-10)
    rft = fft.rfft2(picp)
    Nxf = len(rft[0])
    Nyf = len(rft)
    print(Nx, Nxf)
    print(Ny, Nyf)
    print(rft.sum())
    for i in range(Nxf):
        for j in range(Nyf):
            kx = i*1.0/(Nxf)
            ky = j*1.0/(Nyf)
            k = sp.sqrt(kx*kx+ky*ky)
            rft[j, i] *= sp.exp(-k*k/(0.1*3.0**2))
    print(rft.sum())
    picp = fft.irfft2(rft)
    pic = sp.exp(picp[Nx:2*Nx, Ny:2*Ny])
    return pic
