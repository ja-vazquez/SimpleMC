#!/usr/bin/env python
from RunBase import *

T=LCDMCosmology(Om=1.0,disable_radiation=True)
T2=LCDMCosmology()
zlist=[1/a-1 for a in arange(0.01, 1, 0.01)]
gf1=T.growth_list(zlist)
gf2=T2.growth_list(zlist)
gf2*=gf1[0]/gf2[0]
for z,g1,g2 in zip(zlist, gf1,gf2):
   a=1/(1+z)
   print a,z, g1,g2


