#!/usr/bin/env python

from RunBase import *
from ChainIterator import *

C=ChainIterator('chains_140411_155701','LCDM','phy','BBAO+CMBP')

grlist=[]
wlist=[]
for i in range(0,C.N,1000):
    T=C.theory(i)
    grlist.append(T.growth(10.0))
    wlist.append(C.weight(i))

grlist=array(grlist)
wlist=array(wlist)
mn=(grlist*wlist).sum()/wlist.sum()
er=sqrt((grlist**2*wlist).sum()/wlist.sum()-mn**2)
print "growth z=10/z=0 = ",mn,"+/-",er

