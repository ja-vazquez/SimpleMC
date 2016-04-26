#!/usr/bin/env python
from RunBase import *
L=ParseModel('LCDM')
D=ParseModel('Decay')

print D.Hinv_z(0), L.Hinv_z(0)
print D.Hinv_z(1000), L.Hinv_z(1000)

daL=ParseDataset('BBAO+Planck')
daD=ParseDataset('BBAO+Planck')

daL.setTheory(L)
daD.setTheory(D)


daD.updateParams([Parameter("lambda",0.005)])

print daL.loglike(), daD.loglike()

L.printFreeParameters()
D.printFreeParameters()

