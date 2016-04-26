#!/usr/bin/env python

from RunBase import *
from ChainIterator import *

modell=['Decay', 'Decay01', 'Decay05']

dire = 'chains_SimpleMC/decay/'
data = 'BBAO+Planck'
extra = 'phy'
z = 'z2'
zval=2.0

for model in modell:

  outfile = 'chains/'+model+'_'+z+'_'+extra+'_'+data
  fout=open(outfile+"_1.txt",'w')
  formstr ='%g '*(6) + '\n'

  C=ChainIterator(dire,model,extra,data, balance=False)

  grlist=[]
  wlist=[]
  for i in range(0,C.N,400):
      T=C.theory(i)
      outstr=formstr%tuple([C.weight(i),0.5*C.chi2(i),C.pvalue(i,'h'), C.pvalue(i,'Om'),C.pvalue(i,'lambda') , T.Om_z(1.0/(1+zval))])
      fout.write(outstr)
      fout.flush()	
      #print C.chi2(i), C.weight(i), C.pvalue(i,'h'), T.Om_z(0.5) 
      grlist.append(T.Om_z(0.5))
      wlist.append(C.weight(i))

  fout.close()
  grlist=array(grlist)
  wlist=array(wlist)
  mn=(grlist*wlist).sum()/wlist.sum()
  er=sqrt((grlist**2*wlist).sum()/wlist.sum()-mn**2)
  print "mean =/- 1sigma = ",mn,"+/-",er
