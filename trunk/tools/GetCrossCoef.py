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


def fun(x,y,p1,p2,loglikef):
    loglike = loglikef(x,y)
 
    #if (loglike< -4.5):
    if (loglike<lcut):
        P=0
    elif (loglike>0):
       print "loglike>0?"
       stop()
    else:
        P=exp(loglike)
    return P*x**p1*y**p2



def intg(p1,p2,loglikef):
   #return dblquad(fun, 0.8, 1.2, lambda x:0.8, lambda x:1.2, args=(p1,p2,loglikef))[0]
   ar=arange(0.8,1.2,0.0005)
   toret= array([[fun(x,y,p1,p2,loglikef) for x in ar] for y in ar]).sum()
   print toret
   return toret

#import pylab
#ar=arange(0.8,1.2,0.01)
#pl2d=[[fun(x,y,0,0,datal) for x in ar] for y in ar]
#pylab.imshow(pl2d, extent=(0.8,1.2,0.8,1.2), interpolation='nearest')
#pylab.plot(1.0,1.0,'ro')
#pylab.show()

sigcut=1
lcut=-sigcut**2/2.0
print fun(1.,1.,0.,0.,L.loglike_aperp_apar)

for data in datal:
   print 'Data considered=', data, 'from', datal
   if 'CMASS'==data:
      loglike= L.loglike_aperp_apar
   if 'Auto'==data:
      loglike= J.loglike_aperp_apar
   if 'Cross'==data:
      loglike= K.loglike_aperp_apar
   if 'CM+Au'==data:
      loglike= lambda x,y:L.loglike_aperp_apar(x,y) + J.loglike_aperp_apar(x,y) 
   if 'CM+Cro'==data:
      loglike= lambda x,y:L.loglike_aperp_apar(x,y) + K.loglike_aperp_apar(x,y)  
   if 'Lya'==data:
      loglike= lambda x,y:J.loglike_aperp_apar(x,y) + K.loglike_aperp_apar(x,y)
   if 'CM+Ly'==data:
      loglike= lambda x,y:L.loglike_aperp_apar(x,y) + J.loglike_aperp_apar(x,y)+K.loglike_aperp_apar(x,y) 

   norm = intg(0,0,loglike)
   meanaper = intg(1,0,loglike)/norm
   meanapar = intg(0,1,loglike)/norm
   varaper  = intg(2,0,loglike)/norm-meanaper**2
   varapar  = intg(0,2,loglike)/norm-meanapar**2
   varcross = intg(1,1,loglike)/norm-meanaper*meanapar

   print "------------"
   print "dataset(s)", data
   print "norm=", norm
   print "means=",meanaper, meanapar
   print "var=",sqrt(varaper), sqrt(varapar)
   print "r=", varcross/sqrt(varaper*varapar)
   print "------------" 


