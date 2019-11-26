#!/usr/bin/env python

#
# This little code 
# tests Anderson vs EH approximation
#


## infrastructure
from RunBase import *


### Which datasets do we include in the calculation
for approx in ["EH"]:#"Anderson","EH":
    print "Doing ",approx
    LCDMCosmology.rd_approx=approx
    L=CompositeLikelihood([
        DR11LOWZ(),
        DR11CMASS(),
        DR11LyaAuto(),
        DR11LyaCross()
        ])
    T=oLCDMCosmology()
    L.setTheory(T)
    M=MCMCAnalyzer(L,"chains/test_"+approx)


