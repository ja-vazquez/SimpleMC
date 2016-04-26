#!/usr/bin/env python

#
# This little code demonstrates how to
# use this stuff

## infrastructure
from RunBase import *


### Which datasets do we include in the calculation
L=CompositeLikelihood([
    DR11LOWZ(),
    DR11CMASS(),
    DR11LyaAuto(),
    DR11LyaCross(),
    PlanckLikelihood()
    ])
# note that I could also say something like
# L=ParseDataset('BBAO+Planck')

### My basic theory is LCDM,
T=LCDMCosmology()
# note taht I could also say something like
# T=ParseModel('LCDM')

## Let's see what are my free parameters
T.printFreeParameters()

## Let likelihood know about the theory 
L.setTheory(T)
# so now we can ask it for a value
print "loglike at default parameters = ",L.loglike()

## Then run the MCMC chain
M=MCMCAnalyzer(L,"test_chain", skip=100,nsamp=10000,temp=2)
