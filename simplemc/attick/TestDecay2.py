#!/usr/bin/env python
from RunBase import *
L=ParseModel('LCDM')
D=ParseModel('Decay')

print D.Hinv_z(0), L.Hinv_z(0)
print D.Hinv_z(1000), L.Hinv_z(1000)

daL=ParseDataset('Planck')
daD=ParseDataset('Planck')

daL.setTheory(L)
daD.setTheory(D)

daL.updateParams([Parameter("lambda",0.0), Parameter("Om",0.30)])
daD.updateParams([Parameter("lambda",0.2), Parameter("Om",0.30)])
print daL.loglike(), daD.loglike()

for lam in [0.0,0.1, 0.2,0.3]:
    daD.updateParams([Parameter("lambda",lam), Parameter("Om",0.30-0.18*lam)])
    print lam, daD.loglike(),'XXX'


stop

L.printFreeParameters()
D.printFreeParameters()
print D.rx(0.0), D.rr(0.0)
print D.rx(-0.17), D.rr(-0.17)
print D.rx(-0.47), D.rr(-0.47)
def w(z):
    a=1/(1.+z)
    loga=log(a)
    f=(-3*D.rx(loga)/a**3-4*D.rr(loga)/a**4)/(D.rx(loga)/a**3+D.rr(loga)/a**4)
    print f, f/(-3.0)-1.0

w(0)
w(0.5)
w(2.0)
w(100)
print L.Ocb-L.Obh2/(L.h**2)
