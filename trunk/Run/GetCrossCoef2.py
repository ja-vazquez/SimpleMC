#!/usr/bin/env python
# Input arguments are dataset combination or all

import os, sys
from RunBase import *
from scipy.integrate import *
from TabulatedBAOLikelihood import *

datall=['CMASS', 'Auto', 'Cross', 'CM+Au', 'CM+Cro', 'Lya', 'CM+Ly']

S=[]
if len(sys.argv)>1:
   if 'all' in sys.argv[1]:
      datal=datall
   else:
      for i in arange(1,len(sys.argv),1):
          S.append(sys.argv[i])
      datal=S
else:
   sys.exit('Add a dataset(s)')
print sys.argv[1] 

L=DR11CMASS()
J=DR11LyaAuto()
K=DR11LyaCross()


def fun(x,y,p1,p2,data):
    loglike = 1e-50
    if 'CMASS' in data:
       loglike= L.loglike_aperp_apar(x,y) 
    if 'Auto' in data:
       loglike= J.loglike_aperp_apar(x,y) 
    if 'Cross' in data:
       loglike= K.loglike_aperp_apar(x,y)
    if 'CM+Au' in data:
       loglike= L.loglike_aperp_apar(x,y) + J.loglike_aperp_apar(x,y) 
    if 'CM+Cro' in data:
       loglike= L.loglike_aperp_apar(x,y) + K.loglike_aperp_apar(x,y)  
    if 'Lya' in data:
       loglike= J.loglike_aperp_apar(x,y) + K.loglike_aperp_apar(x,y)
    if 'CM+Ly' in data:
       loglike= L.loglike_aperp_apar(x,y) + J.loglike_aperp_apar(x,y)+K.loglike_aperp_apar(x,y) 
 
    if (loglike< -2):
        P=0
    else:
        P=exp(loglike)
    return P*x**p1*y**p2



def intg(p1,p2,data):
    return dblquad(fun, 0.8, 1.2, lambda x:0.8, lambda x:1.2, args=(p1,p2,data))



for data in datal:
   print 'Data considered=', data, 'from', datal
   norm = intg(0,0,data)[0]
   meanaper = intg(1,0,data)[0]/norm
   meanapar = intg(0,1,data)[0]/norm
   varaper  = intg(2,0,data)[0]/norm-meanaper**2
   varapar  = intg(0,2,data)[0]/norm-meanapar**2
   varcross = intg(1,1,data)[0]/norm-meanaper*meanapar

   print "------------"
   print "dataset(s)", data
   print "norm=", norm
   print "means=",meanaper, meanapar
   print "var=",sqrt(varaper), sqrt(varapar)
   print "r=", varcross/sqrt(varaper*varapar)
   print "------------" 


