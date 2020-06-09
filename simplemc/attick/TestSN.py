#!/usr/bin/env python
from RunBase import *
import warnings
warnings.simplefilter("ignore")

## we need to add hubble parameter likelihood
## otherwise h will be unconstrained
L=CompositeLikelihood([
    BetouleSN(),
    HubbleParameterLikelihood()])

T=LCDMCosmology()

L.setTheory(T)
T.printFreeParameters()
A=MaxLikeAnalyzer(L)
T.printParameters(A.params)



for pl in ['WW','PLA1','PLA2']:
    T=wLCDMCosmology()
    L=CompositeLikelihood([
        BetouleSN(),
        PlanckLikelihood(matrices=pl)])
    L.setTheory(T)
    A=MaxLikeAnalyzer(L)
    T.printParameters(A.params)

stop()
T=wLCDMCosmology()
L=CompositeLikelihood([
    BetouleSN(),
    PlanckLikelihood()])
L.setTheory(T)
A=MaxLikeAnalyzer(L)
T.printParameters(A.params)
#MCMCAnalyzer(L,"chains/test_bet",skip=-1, cov=A.cov)



