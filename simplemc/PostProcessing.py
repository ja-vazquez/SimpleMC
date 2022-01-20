from simplemc.cosmo.Derivedparam import AllDerived
from simplemc.analyzers.dynesty import utils as dyfunc
from simplemc import logger
import numpy as np
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

    """
    def __init__(self, dict_result, paramList, filename,
                 skip=0.1, addDerived=True):
        self.dict_result = dict_result
        self.analyzername = dict_result['analyzer']
        self.result = dict_result['result']
        self.time = dict_result['time']
        self.paramList = paramList
        self.N = len(paramList)
        self.filename = filename
        self.skip = skip
        self.derived = addDerived
        self.args = []
        if addDerived:
            self.AD = AllDerived()

        if self.analyzername in ['mcmc', 'nested', 'emcee']:
            self.maxlogl = np.max(self.result['loglikes'])
        else:
            self.maxlogl = self.result['maxlike']

    def writeSummary(self):
        file = open(self.filename + "_Summary" + ".txt", 'w')
        file.write('SUMMARY\n-------\n')

        for key in self.dict_result:
            if key not in ['result', 'time']:
                file.write('{}: {}\n'.format(key, self.dict_result[key]))

        for key in self.result:
            if key not in ['param_fit', 'samples', 'cov', 'logwt', 'logzerr', 'weights', 'loglikes']:
                if isinstance(self.result[key], float):
                    if key == 'logz':
                        file.write('{}: {:.4f} +/- {:.4f}\n'.format(key, self.result[key], self.result['logzerr']))
                    else:
                        file.write('{}: {:.4f}\n'.format(key, self.result[key]))
                else:
                    file.write('{}: {}\n'.format(key, self.result[key]))

        if self.analyzername in ['mcmc', 'nested', 'emcee']:
            loglikes, samples, weights = self.result['loglikes'], self.result['samples'], self.result['weights']
            means, cov_dy = dyfunc.mean_and_cov(samples, weights)
            stdevs = np.sqrt(np.diag(cov_dy))
            param_fits = means
            print("\nCovariance matrix saved in .covmat file\n {} \n".format(cov_dy))
            np.savetxt(self.filename + ".covmat", cov_dy, delimiter=',')
            try:
                from getdist import mcsamples
                getdistsamples = mcsamples.loadMCSamples(self.filename)
                cov_getdist = getdistsamples.cov(self.paramList)
                print("\ngetdist cov\n {} \n".format(cov_getdist))
            except:
                pass
        else:
            try:
                stdevs = np.sqrt(np.diag(self.result['cov']))
            except:
                stdevs = np.zeros(self.N)
            param_fits = self.result['param_fit']

        for i, parname in enumerate(self.paramList):
            param_fit = param_fits[i]
            std = stdevs[i]
            print("{}: {:.4f} +/- {:.4f}".format(parname, param_fit, std))
            file.write("{}: {:.4f} +/- {:.4f}\n".format(parname, param_fit, std))

        print("\nInformation criterions:\n")
        print("\tAIC: {:.4f}".format(self.aic_criterion()))

        logger.info("\nElapsed time: {:.3f} minutes = {:.3f} seconds".format(self.time / 60, self.time))
        file.write('\nElapsed time: {:.3f} minutes = {:.3f} seconds \n'.format(self.time / 60, self.time))
        file.close()

    def writeMaxlike(self):
        file = open(self.filename + ".maxlike", 'w')
        maxlogl_idx = np.argmax(self.result['loglikes'])

        maxsamp = str(self.result['samples'][maxlogl_idx]).lstrip('[').rstrip(']')
        maxw = self.result['weights'][maxlogl_idx]
        file.write('{} {} {}'.format(maxw, self.maxlogl, maxsamp))
        file.close()

    def mcevidence(self, k):
        if self.analyzername not in ['mcmc', 'nested', 'emcee']:
            sys.exit('MCEvidence only work on Bayesian samplers (mcmc, nested, '
                     'emcee) not in optimizers')

        mcev = MCEvidence('{}'.format(self.filename), kmax=k, dims=self.N)
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
                                  burnlen=burn_frac, kmax=k, dims=self.N)

                mcevres = mcev.evidence(covtype='all')
                if not (mcevres == np.inf).all():
                    valid = True
                if burn_frac > 0.8:
                    print("MCEvidence can't estimate the evidence to your samples")

        return '\nlog-Evidence with mcevidence: {}\n' \
                   'Burn-in fraction: {:.1}\n'.format(mcevres, burn_frac)

    def aic_criterion(self):
        aic = -2*self.maxlogl + 2*(len(self.paramList))
        return aic

    def bic_criterion(self):
        pass

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
