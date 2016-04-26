#!/usr/bin/env python

import sys

sys.path.append('Run')
#sys.path.append('../Run')

from RunBase import *
from ChainIterator import *
from math import *

def contour(chaindir,cosmo,data,compute_param_function) :

    C=ChainIterator(chaindir,cosmo,'phy',data)
    z=numpy.arange(0,3,0.01)

    # find best model 
    min_chi2 = numpy.min(C.chain.chain[:,1])
    model_best_index = numpy.argmin(C.chain.chain[:,1])
    print "minimum chi2 = ", min_chi2
    T=C.theory(model_best_index)

    # compute param for this model
    param_best=numpy.zeros(z.shape)
    for i in range(z.shape[0]) :
        param_best[i]=compute_param_function(T,z[i])

    # allocate param for +- 1 sigma
    param_1sigma_plus=param_best.copy()
    param_1sigma_minus=param_best.copy()


    # selection of models within 1sig
    model_1sig_indices=numpy.where(C.chain.chain[:,1]<=min_chi2+1)[0]
    print "number of models at 1 sig = ",model_1sig_indices.shape[0]
    # if too many, look only at those close to 1 sigma (this is for LCDM, there are too many)
    chi2_floor=min_chi2
    while (model_1sig_indices.shape[0]>3000) :
        chi2_floor+=0.01
        model_1sig_indices=numpy.where((C.chain.chain[:,1]<=min_chi2+1)&(C.chain.chain[:,1]>=chi2_floor))[0]
    print "number of selected models for the 1sig contour plot = ",model_1sig_indices.shape[0]
    


    # loop on models to find max excursion 
    count=0
    for index in model_1sig_indices :
        count+=1
        T=C.theory(index)
        for i in range(z.shape[0]) :
            parami=compute_param_function(T,z[i])
            param_1sigma_minus[i]=min(param_1sigma_minus[i],parami)
            param_1sigma_plus[i]=max(param_1sigma_plus[i],parami)
        if count%100==0 :
            print "done %d/%d"%(count,model_1sig_indices.shape[0])
            
    # return -1sig,best,1sig
    return z,param_1sigma_minus,param_best,param_1sigma_plus


