#
# This is the MCMC module.
# it spits out chains that are compatible with CosmoMC
# it calculates cov matrix during burn-in.
# chain_num tells it to spit out multi-node chains.
# optional temperature makes it sample at a higher temperature but note that
# this guy, as opposed to cosmomc, reweights the weights on the fly.
#

from random import *
from math import *
from scipy import *
import scipy.linalg as la
import copy
import random
import sys
import os.path as path


class MCMCAnalyzer:
    def __init__(self, like, outfile, skip=5000, nsamp=100000, temp=1.0, cov=None, chain_num=None):

        self.like = like
        self.outfile = outfile
        self.nsamp = nsamp
        self.skip = skip
        self.temp = float(temp)  # temperature
        self.chain_num = chain_num
        self.cpars = like.freeParameters()
        minvals, maxvals = [], []
        for lb, hb in [p.bounds for p in self.cpars]:
            minvals.append(lb)
            maxvals.append(hb)
        self.minvals = array(minvals)
        self.maxvals = array(maxvals)
        print("Bounds:", self.minvals, self.maxvals)

        self.N = len(self.cpars)

        if (like.name() == "Composite"):
            self.sublikenames = like.compositeNames()
            self.composite = True
        else:
            self.composite = False

        if (cov == None):
            # make initial cov matrix from diagonal "errors"
            errs = [0.01*p.error**2 for p in self.cpars]
            self.init_pcov(diag(errs))
        else:
            self.init_pcov(cov)

        self.RunChain()

    def RunChain(self):
        self.openFiles()
        self.cloglike, self.cloglikes = self. getLikes()
        # set up logofs based on the first log like which should be
        # the same for all chains. Better than nothing.
        # self.logofs=self.cloglike
        # Actually, above doesn't seem to work very well. Instead, use zero, as our likelihoods never became very large
        self.logofs = 0
        # current weight
        self.cw = 0
        # current counter
        self.co = 0
        # mean for burin
        self.swx = 0
        self.meanx = zeros(self.N)
        self.meanxx = zeros((self.N, self.N))
        # max loglike
        self.maxloglike = -1e30
        # are we done
        self.done = False
        print("Starting chain...")

        while not (self.done):
            ppars, numout = self.GetProposal()
            self.cw += numout  ## things hitting outside the prior are formally rejected samples
            self.like.updateParams(ppars)
            ploglike, ploglikes = self.getLikes()
            if (isnan(ploglike)):
                print("Something bad has happened, nan in loglike, assuming zero log")
                ploglike = -1e50
            # print cloglike, ploglike, [p.value for p in like.freeParameters()], [p.value for p in self.cpars]
            if (ploglike > self.cloglike):
                accept = True
            else:
                accept = (exp((ploglike-self.cloglike)/self.temp)
                          > uniform(0., 1.))

            # print [p.value for p in ppars], accept, ploglike
            # stop
            if (accept):
                self.ProcessAccepted(ppars, ploglike, ploglikes)
            else:
                self.cw += 1 
        self.closeFiles()

    def GetProposal(self):
        vec = zeros(self.N)
        numreject=0
        while True:
            ppars = copy.deepcopy(self.cpars)
            step = self.draw_pcov()
        # print step# [p.value for p in  step]
            for i, p in enumerate(ppars):
                p.value += step[i]
                vec[i] = p.value

            if all(vec > self.minvals) and all(vec < self.maxvals):
                return ppars, numreject
            numreject+=1
            
    def init_pcov(self, mat):
        self.chol = la.cholesky(mat)

    def draw_pcov(self):
        a = array([random.gauss(0., 1,) for i in range(self.N)])
        return dot(a, self.chol)

    def openFiles(self):
        outfile = self.outfile
        if self.chain_num in [None, 1]:
            fpar = open(outfile+".paramnames", 'w')
            for p in self.cpars:
                fpar.write(p.name+"\t\t\t"+p.Ltxname+"\n")
            if self.composite:
                for name in self.sublikenames:
                    fpar.write(name+"_like \t\t\t"+name+"\n")
                fpar.write("theory_prior \t\t\t None \n")
            fpar.close()

        formstr = '%g '+'%g '*(self.N+1)
        if (self.composite):
            formstr += '%g '*(len(self.sublikenames)+1)
        formstr += '\n'

        if (self.chain_num == None):
            cfname = outfile+".txt"
            mlfname = outfile+".maxlike"
        else:
            cfname = outfile+"_%i.txt" % (self.chain_num)
            mlfname = outfile+"_%i.maxlike" % (self.chain_num)

        if (path.isfile(cfname)):
            print("Due to bad habits in the past, won't open existing file.", cfname)
            sys.exit(1)
        self.fout = open(cfname, 'w')
        self.mlfout = open(mlfname, 'w')

        self.formstr = formstr

    def closeFiles(self):
        self.fout.close()
        self.mlfout.close()

    def getLikes(self):
        if (self.composite):
            cloglikes = self.like.compositeLogLikes_wprior()
            cloglike = cloglikes.sum()

        else:
            cloglikes = []
            cloglike = self.like.loglike_wprior()

        return cloglike, cloglikes

    def ProcessAccepted(self, ppars, ploglike, ploglikes):
        self.co += 1
        if (self.co % 1000 == 0):
            print("Accepted samples", self.co, self.cw)
        vec = [p.value for p in self.cpars]

        if (self.co > self.skip):
            # weight rescaled
            wers = self.cw*exp((self.cloglike-self.logofs)
                               * (self.temp-1.0)/self.temp)

            if (self.composite):
                outstr = self.formstr % tuple(
                    [wers, -self.cloglike]+vec + self.cloglikes.tolist())
            else:
                outstr = self.formstr % tuple([wers, -self.cloglike]+vec)

            self.fout.write(outstr)
            # Flush file on regular basis
            if (self.co % 1000 == 0):
                self.fout.flush()

            if (self.cloglike > self.maxloglike):
                self.maxloglike = self.cloglike
                print("New maxloglike", self.maxloglike)
                self.mlfout.seek(0)
                self.mlfout.write(outstr)
                self.mlfout.flush()

            if self.co > self.nsamp:
                self.done = True

        elif (self.co < self.skip):
            self.swx += self.cw
            v = array(vec)
            self.meanx += v*self.cw
            self.meanxx += outer(v, v)*self.cw
            if (self.cw > 30):
                print("Still burning in, weight too large")
                self.chol *= 0.9
                print(self.cw)
        else:  # co==skip
            self.meanx /= self.swx
            self.meanxx /= self.swx
            self.meanxx -= outer(self.meanx, self.meanx)
            print("Re-initializing covariance matrix after burn-in")
            print(self.meanxx)
            for i, p in enumerate(self.cpars):
                print(p.name, p.value, sqrt(self.meanxx[i, i]))

            self.init_pcov(self.meanxx)

        self.cw = 1
        self.cpars = ppars
        self.cloglike = ploglike
        if self.composite:
            self.cloglikes = ploglikes
