from simplemc.cosmo.Derivedparam import AllDerived
from simplemc.analyzers.dynesty import utils as dyfunc
from simplemc import logger
import numpy as np
import re
import sys
from .analyzers import MCEvidence


class PostProcessing:
    """
    This class makes postprocessing such as generate a summary with some statistics.

    Parameters
    ---------
    list_result : list
        List with results from sampling.
    paramList : list
        List with Parameter objects.
    filename : str.
        File name.
    skip : float
        Burn-in.
    addDerived : bool
        Derived parameters?
    loglike : object
        Likelihood object.

    """
    def __init__(self, list_result, paramList, filename,
                 skip=0.1, addDerived=True, loglike=None):
        self.analyzername = list_result[0]
        self.result    = list_result[1]
        self.paramList = paramList
        self.N = len(paramList)
        self.filename  = filename
        self.skip      = skip
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

        pars = self.loglike.freeParameters()
        if self.analyzername == 'nested':
            samples, weights = self.result.samples, np.exp(self.result.logwt - self.result.logz[-1])
        elif self.analyzername == 'mcmc':
            samples, weights = self.result[0], self.result[1]
        elif self.analyzername == 'emcee':
            samples = self.result.get_chain(flat=True)
            weights = np.ones(len(samples))
        else:
            samples = None
            weights = None

        if self.analyzername in ['mcmc', 'nested', 'emcee']:
            means, cov = dyfunc.mean_and_cov(samples, weights)
            stdevs = np.sqrt(np.diag(cov))

            for i, p in enumerate(pars):
                mean = means[i]
                std = stdevs[i]
                print("{}: {:.4f} +/- {:.4f}".format(p.name, mean, std))
                file.write("{}: {:.4f} +/- {:.4f}\n".format(p.name, mean, std))

            if self.analyzername == 'nested':
                file.write("nlive: {:d}\nniter: {:d}\nncall: {:d}\n"
                           "eff(%): {:6.3f}\nlogz: "
                           "{:6.3f} +/- {:6.3f}\n".format(self.result.nlive, self.result.niter,
                                                   sum(self.result.ncall), self.result.eff,
                                                   self.result.logz[-1], self.result.logzerr[-1]))


        logger.info("\nElapsed time: {:.3f} minutes = {:.3f} seconds".format(time / 60, time))
        file.write('\nElapsed time: {:.3f} minutes = {:.3f} seconds \n'.format(time / 60, time))
        file.close()

    def mcevidence(self):
        if self.analyzername not in ['mcmc', 'nested', 'emcee']:
            sys.exit('MCEvidence only work on Bayesian samplers (mcmc, nested, '
                     'emcee) not in optimizers')

        mcev = MCEvidence('{}'.format(self.filename))
        mcevres = mcev.evidence(covtype='all')

        burn_frac = 0.0

        if (mcevres == np.inf).all():
            logger.info("MCEvidence failed to calculate Bayesian evidence,\n"
                        " it is trying again.")
            valid = False
            while not valid:
                burn_frac += 0.1
                logger.info("Burn-in: {}%".format(burn_frac*100))
                mcev = MCEvidence('{}'.format(self.filename),
                                  burnlen=burn_frac)

                mcevres = mcev.evidence(covtype='all')
                if not (mcevres == np.inf).all():
                    valid = True
                if burn_frac > 0.8:
                    print("MCEvidence can't estimate the evidence to your samples")

        return '\nlog-Evidence with mcevidence: {}\n' \
                   'Burn-in fraction: {:.1}\n'.format(mcevres, burn_frac)

    def plot(self, chainsdir, show=False):
        """
        Simple connection with the plotters.

        Parameters
        -----------
        show : bool
            Default False
        """
        from .plots.SimplePlotter import SimplePlotter
        figure = SimplePlotter(chainsdir, self.paramList, path=self.filename, show=show)

        return figure
