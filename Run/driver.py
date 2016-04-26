#!/usr/bin/env python
import sys
from RunBase import *

#
# This little code allows running something
# based on command line. To be used in conjuction wiht 
# the wqdriver that subtmis to BNL queue
#

if (len(sys.argv)<4):
    print """Usage: 

Run/driver.py pre/phy model dataset [chain #] [# to skip] [# to take] [directory]

where model is one of 

%s 

and datasets is one of 

%s

Example:

Run/driver.py phy LCDM BBAO+Planck

"""%(model_list, data_list)
    sys.exit(1)


prefact,model,datasets=sys.argv[1:4]
if len(sys.argv)>4:
    chainno=int(sys.argv[4])
else:
    chainno=1

if len(sys.argv)>5:
    skip=int(sys.argv[5])
else:
    skip=5000

if len(sys.argv)>6:
    nsamp=int(sys.argv[6])
else:
    nsamp=100000

if len(sys.argv)>7:
    chainsdir=sys.argv[7]
else:
    chainsdir = 'chains/'

temp=2.0 ## temperature at which to sample, weights get readjusted on the fly

print "Running ",model, prefact,datasets, skip, nsamp, temp, chainno
T=ParseModel(model)
L=ParseDataset(datasets)
if prefact=="pre":
    T.setVaryPrefactor()
T.printFreeParameters()
L.setTheory(T)

if chainno>0:
    M=MCMCAnalyzer(L, chainsdir + "/" + model+"_"+prefact+"_"+datasets,skip=skip,nsamp=nsamp, temp=temp, chain_num=chainno)
else:
    A=MaxLikeAnalyzer(L)
    T.printParameters(A.params)
