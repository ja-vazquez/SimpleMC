#!/usr/bin/env python

#
# This little plots 
# results of TestAndVSEH.py
#

import pylab, cosmich



C1=cosmich.cosmochain('chains/test_EH.txt',[''])
C2=cosmich.cosmochain('chains/test_Anderson.txt',[''])
nams=[x[:-1] for x in open('chains/test_EH.params').readlines()]
N=len(nams)
cc=0
pylab.figure(figsize=(4*N,4*N))
for i in range(N):
    for j in range(N):
        cc+=1
        if (i<j):
            continue
        pylab.subplot(N,N,cc)
        if (i==j):
            xv,yv=C1.GetHisto(i+2)
            xv2,yv2=C2.GetHisto(i+2)
            pylab.plot(xv,yv,'r-')
            pylab.plot(xv2,yv2,'b-')
        elif (i>j):
            C1.Plot2D(j+2,i+2,filled=2)
            C2.Plot2D(j+2,i+2,filled=0)
            if (j==0):
                pylab.ylabel(nams[i])
        if (i==N-1):
            pylab.xlabel(nams[j])

pylab.savefig("EHvsAnd.pdf")

