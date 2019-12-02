"""
This module processes the samples from a nested sampler and prints, saves chains in a text file 
and creates a param file.
"""
import sys
sys.path=["pybambi","../pybambi"]+sys.path
import os.path as path
from os import remove
import numpy as np
import nestle


def printNestle(result,dims,listParameters):
    """
    This method prints the results of parameter estimation with Nestle samplers.

    Parameters:

    result:         object of nestle
    dims:           Integer. Free parameters of the model
    listParameters: list [] of the free parameters
    """
    p, cov = nestle.mean_and_cov(result.samples, result.weights)
    outtext = []
    for i in range(dims):
        outtext.append( "%s = %5f +/- %5f"%(listParameters[i],p[i], np.sqrt(cov[i, i])))
        print(outtext[i])
    return (outtext)

def saveChainNestle(result,outputname):
    """
    This generates an output Simple(cosmo)MC style for Nestle Samplers.
    
    Parameters:

    result:     Nestle object
    outputname: str name of the output file

    """
    if(path.isfile(outputname+'.txt')):
    	print("An existing file with the same name has been deleted.", outputname+'.txt')
    	remove(outputname+'.txt')

    f = open(outputname+'.txt','a+')

    for i in range(len(result.samples)):
        strweights = str(result.weights[i]).lstrip('[').rstrip(']')
        strlogl=str(-1*result.logl[i]).lstrip('[').rstrip(']')
        strsamples=str(result.samples[i]).lstrip('[').rstrip(']')
        row = strweights+' '+strlogl+' '+strsamples
        nrow = " ".join( row.split() )
        f.write(nrow+'\n')
        #f.write(strweights+' '+strlogl+' '+strsamples+'\n')
    f.close()

def paramFiles(T,L,outputname):
    """
    This method writes the .paramnames file with theirs LaTeX names.

    Parameters: 

    T:          T is an instance of ParseModel(model)
    L:          L is an instance of ParseDataset(datasets)
    outputname  str, name of the output file

    """
    cpars = T.freeParameters()
    parfile = outputname.strip(".txt")+".paramnames"
    
    if (path.isfile(parfile)):
        print("Existing parameters file!")
    
    fpar=open(parfile,'w')   
    for p in cpars:
        fpar.write(p.name+"\t\t\t"+p.Ltxname+"\n")
    fpar.close()