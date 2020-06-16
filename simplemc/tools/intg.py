#!/usr/bin/env python
from RunBase import *
C=LCDMCosmology()

Df=C.Da_z(1090)

lims=[0.25,0.50, 0.95, 0.99, 0.995]

l=lims.pop(0)
z=0.5;
while True:
    f=C.Da_z(z)/Df
    if (f>l):
        print z,f, "XXX"
        if len(lims)>0:
            l=lims.pop(0)
        else:
            break
    z*=1.03


