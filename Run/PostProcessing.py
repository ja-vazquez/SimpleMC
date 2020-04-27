"""
This module processes the samples from a nested sampler and prints, saves chains in a text file 
and creates a param file.
"""
import sys
sys.path=["pybambi" , "../pybambi"]+sys.path
import os.path as path
from os import remove
import numpy as np
import dynesty

class PostProcessing():
    def __init__(self, list_result, paramList, outputname,\
                 olambda = True, chainsdir = 'chains', skip = 0.1, engine = 'dynesty'):
        self.analyzername = list_result[0]
        self.result = list_result[1]
        self.paramList = paramList
        self.filename = chainsdir + '/' + outputname
        self.args = []
        self.list_result = list_result
        self.skip = skip
        self.olambda = olambda
        self.engine = engine

        for i in range(2, len(self.list_result)):
            self.args.append(self.list_result[i])
              

    def saveNestedChain(self):
        """
        This generates an output Simple(cosmo)MC style for dynesty Samplers.

        """

        f = open(self.filename + '_1.txt', 'w+')
        if self.engine == 'dynesty':
            weights = np.exp(self.result['logwt'] - self.result['logz'][-1])
<<<<<<< HEAD
            postsamples = self.result.samples
       
=======
            postsamples = dynesty.utils.resample_equal(self.result.samples, weights)
            #postsamples =self.result.samples
>>>>>>> 33667434ec3900d2c4380db33a366dd782e1a47b
            print('\n Number of posterior samples is {}'.format(postsamples.shape[0]))

        
            for i, sample in enumerate(postsamples):
                strweights = str(weights[i])
                strlogl = str(self.result['logl'][i])
                strsamples = str(sample).lstrip('[').rstrip(']')
                #if self.olambda:
                #    OLambda = 1 - sample[0]
                #    strOLambda = ' '+str(OLambda)
                #else:
                #    strOLambda = ""
                row = strweights + ' ' + strlogl + ' ' + strsamples #+ strOLambda
                nrow = " ".join( row.split() )
                f.write(nrow+'\n')

        elif self.engine == 'nestle':
            for i in range(len(self.result.samples)):
                strweights = str(self.result.weights[i]).lstrip('[').rstrip(']')
                strlogl=str(-1*self.result.logl[i]).lstrip('[').rstrip(']')
                strsamples=str(self.result.samples[i]).lstrip('[').rstrip(']')
                row = strweights+' '+strlogl+' '+strsamples
                nrow = " ".join( row.split() )
                f.write(nrow+'\n')

        f.close()

    #AJUSTAR!
    def saveEmceeSamples(self):
        dims = len(self.paramList)
        postsamples = self.result.chain[:, self.skip:, :].reshape((-1, dims))
        f = open(self.filename + '.txt', 'a+')
        for i, row in enumerate(postsamples):
            strsamples = str(row).lstrip('[').rstrip(']')
            f.write(strsamples+'\n')
        # tau = self.result.get_autocorr_time()
        # burnin = int(2 * np.max(tau))
        # thin = int(0.5 * np.min(tau))
        # samples = self.result.get_chain(discard=burnin, flat=True, thin=thin)
        # log_prob_samples = self.result.get_log_prob(discard=burnin, flat=True, thin=thin)
        # log_prior_samples = self.result.get_blobs(discard=burnin, flat=True, thin=thin)
        # all_samples = np.concatenate((samples, log_prob_samples[:, None],
        #               log_prior_samples[:, None]), axis=1)
        # for i in all_samples:
        #     print(i)
        f.close()
    def paramFiles(self, T, L):
        """
        This method writes the .paramnames file with theirs LaTeX names.

        Parameters: 

        T:          T is an instance of ParseModel(model)
        L:          L is an instance of ParseDataset(datasets)

        """
        cpars = T.freeParameters()
        parfile = self.filename + ".paramnames"
        
        if (path.isfile(parfile)):
            print("Existing parameters file!")
        
        fpar=open(parfile, 'w')   
        for p in cpars:
            fpar.write(p.name+"\t\t\t"+p.Ltxname+"\n")
        #if self.olambda:
        #    fpar.write("OLambda\t\t\t \Omega_{\Lambda}\n")
        #fpar.close()

    def writeSummary(self, time, *args):
        file = open(self.filename + "_Summary" + ".txt",'w')
        file.write('SUMMARY\n')
        
        for item in self.args:
            if type(item) is list:
                for element in item:
                    file.write(str(element)+'\n')
            else:
                file.write(str(item)+'\n')

        for item in args:
            if type(item) is list:
                for element in item:
                    file.write(str(element)+'\n')
            else:
                file.write(str(item)+'\n')


        print("\nElapsed time: %.3f minutes = %.3f seconds"%(time/60,time))  
        file.write('\nElapsed time: %.3f minutes = %.3f seconds \n'%(time/60,time))  
        file.close()

    def getdistAnalyzer(self, cov=False):
        #import getdist
        from getdist import mcsamples
        
        mcsamplefile = mcsamples.loadMCSamples(self.filename, settings={'ignore_rows':self.skip})

        if cov:
            cov = mcsamplefile.cov(pars = self.paramList)
            np.savetxt(self.filename + '_' + 'cov.txt', cov)
            print("Covariance matrix:\n")
            print(cov)
            print("\n")

        means = mcsamplefile.getMeans()

        stddev = mcsamplefile.std(self.paramList)
        
        summaryResults = []

        for i, param in enumerate(self.paramList):
            print(self.paramList[i] + " : " + str(round(means[i], 4)) + \
                "+/-" + str(round(stddev[i], 4)))
            summaryResults.append(self.paramList[i] + " : " + str(round(means[i], 4)) + \
                                     "+/-" + str(round(stddev[i], 4)))

        #if self.olambda:
        #    summaryResults.append("OLambda : " + str(round(1 - means[0], 4)) \
        #        + "+/-" + str(round(stddev[0], 4)))
        
        #print("summaryResults", summaryResults[0])
        #mcstats.likeSummary()
<<<<<<< HEAD
        return summaryResults
=======
        return summaryResults


    def saveChainNestle(self):
        """
        This generates an output Simple(cosmo)MC style for Nestle Samplers.

        Parameters:

        result:     Nestle object
        outputname: str name of the output file

        """
        if(path.isfile(outputname+'.txt')):
            print("There's a file with the same name.", outputname+'.txt')
            #outputname+= '_tmp'
            #remove(outputname+'.txt')
            #sys.exit(1)

        f = open(outputname+'_1.txt','a+')

        for i in range(len(result.samples)):
            strweights = str(result.weights[i]).lstrip('[').rstrip(']')
            strlogl=str(-1*result.logl[i]).lstrip('[').rstrip(']')
            strsamples=str(result.samples[i]).lstrip('[').rstrip(']')
            row = strweights+' '+strlogl+' '+strsamples
            nrow = " ".join( row.split() )
            f.write(nrow+'\n')
            #f.write(strweights+' '+strlogl+' '+strsamples+'\n')
        f.close()
>>>>>>> 33667434ec3900d2c4380db33a366dd782e1a47b
