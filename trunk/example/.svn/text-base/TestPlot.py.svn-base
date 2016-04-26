#!/usr/bin/env python

from cosmich import *
import pylab

C=cosmochain('test_chain.txt',None)

pylab.subplot(121)
C.Plot2D("h","Om",N=25)
pylab.xlabel("$h$")
pylab.ylabel("$\Omega_m$")

pylab.subplot(122)
C.Plot1D('h',N=40)
pylab.xlabel("$h$")
pylab.ylabel("prob")

print "Plotting to testplot.pdf"
pylab.savefig("testplot.pdf")

