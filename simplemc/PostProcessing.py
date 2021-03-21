"""
This module processes the samples from a nested sampler and prints, saves chains in a text file
and creates a param file.
"""
# from simplemc.tools.Simple_Plots import Simple_plots
from simplemc.cosmo.Derivedparam import AllDerived
from simplemc.analyzers.dynesty import utils as dyfunc
from simplemc import logger
import numpy as np
import re

#TODO EMCEE

class PostProcessing:
    """
       In this class we...
    """
    def __init__(self, list_result, paramList, filename, \
                 skip=0.1, engine=None, addDerived=True, loglike=None):
        self.analyzername = list_result[0]
        self.result    = list_result[1]
        self.paramList = paramList
        self.filename  = filename
        self.skip      = skip
        self.engine    = engine
        self.derived   = addDerived
        self.loglike   = loglike
        self.args = []
        if addDerived:
            self.AD = AllDerived()

        for i in range(2, len(list_result)):
            self.args.append(list_result[i])



    def writeSummary(self, time, *args):
        file = open(self.filename + "_Summary" + ".txt", 'w')
        file.write('SUMMARY\n-------\n')

        for item in self.args:
            if type(item) is list:
                for element in item:
                    if element is not None:
                        file.write(str(element) + '\n')
            else:
                if item is not None:
                    file.write(str(item) + '\n')

        for item in args:
            if type(item) is list:
                for element in item:
                    if element is not None:
                        file.write(str(element) + '\n')
            else:
                if item is not None:
                    file.write(str(item) + '\n')

        if self.engine =='dynesty':
            pars = self.loglike.freeParameters()
            samples, weights = self.result.samples, np.exp(self.result.logwt - self.result.logz[-1])
            means, cov = dyfunc.mean_and_cov(samples, weights)
            stdevs = np.sqrt(np.diag(cov))
            
            for i, p in enumerate(pars):
                mean = means[i]
                std = stdevs[i]
                print("{}: {:.4f} +/- {:.4f}".format(p.name, mean, std))
                file.write("{}: {:.4f} +/- {:.4f}\n".format(p.name, mean, std))

            file.write("nlive: {:d}\nniter: {:d}\nncall: {:d}\n"
                       "eff(%): {:6.3f}\nlogz: "
                       "{:6.3f} +/- {:6.3f}".format(self.result.nlive, self.result.niter,
                                               sum(self.result.ncall), self.result.eff,
                                               self.result.logz[-1], self.result.logzerr[-1]))

        logger.info("\nElapsed time: {:.3f} minutes = {:.3f} seconds".format(time / 60, time))
        file.write('\nElapsed time: {:.3f} minutes = {:.3f} seconds \n'.format(time / 60, time))
        file.close()


    def getdistAnalyzer(self, cov=False):
        from getdist import mcsamples

        mcsamplefile = mcsamples.loadMCSamples(self.filename, settings={'ignore_rows': self.skip})

        if cov:
            cov = mcsamplefile.cov(pars=self.paramList)
            np.savetxt(self.filename + '_' + 'cov.txt', cov)
            logger.info("Covariance matrix:\n")
            logger.info(cov)
            logger.info("\n")

        means  = mcsamplefile.getMeans()
        stddev = mcsamplefile.std(self.paramList)
        summaryResults = []

        for i, param in enumerate(self.paramList):
            logger.info(self.paramList[i] + " : " + str(round(means[i], 4)) + \
                        "+/-" + str(round(stddev[i], 4)))
            summaryResults.append(self.paramList[i] + " : " + str(round(means[i], 4)) + \
                                  "+/-" + str(round(stddev[i], 4)))
        return summaryResults



    # AJUSTAR!
    def saveEmceeSamples(self, thin=1):
        f = open(self.filename + '.txt', 'w+')
        logprobs = self.result.get_log_prob(discard=self.skip, flat=True, thin=thin)
        postsamples = self.result.get_chain(discard=self.skip, flat=True, thin=thin)

        for i, row in enumerate(postsamples):
            strsamples = str(row).lstrip('[').rstrip(']')
            strsamples = "{} {} {}\n".format(1, -2 * (logprobs[i]), strsamples)
            strsamples = re.sub(' +', ' ', strsamples)
            strsamples = re.sub('\n ', ' ', strsamples)
            if self.derived:
                for pd in self.AD.listDerived(self.loglike):
                    strsamples = "{} {}".format(strsamples, pd.value)
            f.write(strsamples)
        f.close()


