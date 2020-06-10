#!/usr/bin/env python

from RunBase import *
import time

T=ParseModel("LCDM")
#T.setMnu(0.0)
L=ParseDataset("BBAO+CMBP+SN")#+CMBP")
T.printFreeParameters()

L.setTheory(T)
print T.WangWangVec()

t0 = time.time()
for i in range(30):
    print i
    loglike=L.loglike()
t= time.time()-t0
print loglike,t
