#
#
# This module importance samples a cosmomc chain.
#
#
from random import *
from math import *
import numpy as np
import os.path as path
from Parameter import *


class CosmoMCImportanceSampler:
    def __init__(self, inroot, outroot, like):
        self.dtype = []
        self.dtype.append(('weight', 'f8'))
        self.dtype.append(('nloglike', 'f8'))
        op = open(outroot+".paramnames", 'w')
        for line in open(inroot+".paramnames").readlines():
            op.write(line)
            name = line.split()[0]
            self.dtype.append((name, 'f8'))
        op.close()
        self.inroot = inroot
        self.outroot = outroot
        self.ofset = None
        self.names = [d[0] for d in self.dtype]
        self.like = like

    def run(self):
        c = 1
        while True:
            fname = self.inroot+"_%i.txt" % c
            ofname = self.outroot+"_%i.txt" % c
            if (path.isfile(fname)):
                print("Processing ", fname)
                self.importanceSample(fname, ofname)
                c += 1
            else:
                break
        print("Done.")

    def importanceSample(self, fname, ofname):
        chain = np.loadtxt(fname, dtype=self.dtype)
        llike = []
        N = len(chain)
        for i, line in enumerate(chain):
            if i % 1000 == 0:
                print("%i/%i..." % (i, N))
            self.like.updateParams(self.cosmomc2april(line))
            llike.append(self.like.loglike())
        llike = np.array(llike)
        if self.ofset is None:
            self.ofset = llike.max()
        rewe = np.exp(llike-self.ofset)
        chain['weight'] *= rewe
        print("min/mean/max weight=", rewe.min(), rewe.mean(), rewe.max())
        np.savetxt(ofname, chain)

    def cosmomc2april(self, line):
        plist = [Parameter("Obh2", line['omegabh2']),
                 Parameter("Om", line['omegam*']),
                 Parameter("h", line['H0*']/100.)]
        if "nnu" in self.names:
            plist.append(Parameter("Nnu", line['nnu']))
        if "w" in self.names:
            plist.append(Parameter("w", line['w']))
        if "wa" in self.names:
            plist.append(Parameter("wa", line['wa']))
        if "omegak" in self.names:
            plist.append(Parameter("Ok", line['omegak']))
        if "mnu" in self.names:
            plist.append(Parameter("mnu", line['mnu']))
        return plist
