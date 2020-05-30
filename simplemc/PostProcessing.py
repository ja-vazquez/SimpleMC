"""
This module processes the samples from a nested sampler and prints, saves chains in a text file
and creates a param file.
"""
# from simplemc.tools.Simple_Plots import Simple_plots
import sys
import os.path as path
import numpy as np
import re
import scipy as sp


class PostProcessing():
    """
    In this class we...
    """

    def __init__(self, list_result, paramList, filename, \
                 skip=0.1, engine='dynesty', addDerived=True, loglike=None):
        self.analyzername = list_result[0]
        self.result = list_result[1]
        self.paramList = paramList
        self.filename = filename
        self.args = []
        self.list_result = list_result
        self.skip = skip
        self.engine = engine
        self.derived = addDerived
        self.loglike = loglike
        if self.derived:
            self.AD = AllDerived()

        for i in range(2, len(self.list_result)):
            self.args.append(self.list_result[i])

    def saveNestedChain(self):
        """
        This generates an output Simple(cosmo)MC style for dynesty Samplers.

        """
        if path.isfile(self.filename + '_1.txt'):
            logger.critical("Output file exists! Please choose another"
                            " name or move the existing file.")
            sys.exit(1)
        else:
            f = open(self.filename + '_1.txt', 'w+')
        if self.engine == 'dynesty':
            weights = np.exp(self.result['logwt'] - self.result['logz'][-1])

            postsamples = self.result.samples

            logger.info('\n Number of posterior samples is {}'.format(postsamples.shape[0]))

            for i, sample in enumerate(postsamples):
                strweights = str(weights[i])
                strlogl = str(self.result['logl'][i])
                strsamples = str(sample).lstrip('[').rstrip(']')
                row = strweights + ' ' + strlogl + ' ' + strsamples  # + strOLambda
                nrow = " ".join(row.split())
                if self.derived:
                    for pd in self.AD.listDerived(self.loglike):
                        nrow += " " + str(pd.value)
                f.write(nrow + '\n')

        elif self.engine == 'nestle':
            for i in range(len(self.result.samples)):
                strweights = str(self.result.weights[i]).lstrip('[').rstrip(']')
                strlogl = str(-1 * self.result.logl[i]).lstrip('[').rstrip(']')
                strsamples = str(self.result.samples[i]).lstrip('[').rstrip(']')
                row = strweights + ' ' + strlogl + ' ' + strsamples
                nrow = " ".join(row.split())
                f.write(nrow + '\n')
        f.close()

    # AJUSTAR!
    def saveEmceeSamples(self):
        dims = len(self.paramList)
        # postsamples = self.result.chain[:, self.skip:, :].reshape((-1, dims))
        tau = self.result.get_autocorr_time()
        logger.info("Autocorrelation time: {}".format(tau))
        tau = np.mean(tau)
        burnin = int(0.5 * np.max(tau))
        thin = int(0.5 * np.min(tau))
        f = open(self.filename + '.txt', 'w+')
        logprobs = self.result.get_log_prob(discard=burnin, flat=True, thin=thin)
        flat_log_prior_samps = self.result.get_blobs(flat=True)
        postsamples = self.result.get_chain(discard=burnin, flat=True, thin=thin)

        for i, row in enumerate(postsamples):
            strsamples = str(row).lstrip('[').rstrip(']')
            # strsamples = "{} {} {}\n".format(1, -2 * (logprobs[i] - flat_log_prior_samps[i]), strsamples)
            strsamples = "{} {} {}\n".format(1, -2 * (logprobs[i]), strsamples)
            strsamples = re.sub(' +', ' ', strsamples)
            strsamples = re.sub('\n ', ' ', strsamples)
            f.write(strsamples)
        f.close()

    def paramFiles(self):
        """
        This method writes the .paramnames file with theirs LaTeX names.

        Parameters:

        T:          T is an instance of ParseModel(model)
        L:          L is an instance of ParseDataset(datasets)

        """
        cpars = self.loglike.freeParameters()
        parfile = self.filename + ".paramnames"

        if (path.isfile(parfile)):
            logger.info("Existing parameters file!")

        fpar = open(parfile, 'w')
        for p in cpars:
            fpar.write(p.name + "\t\t\t" + p.Ltxname + "\n")
        if self.derived:
            for pd in self.AD.list:
                fpar.write(pd.name + "\t\t\t" + pd.Ltxname + "\n")

    def writeSummary(self, time, *args):
        file = open(self.filename + "_Summary" + ".txt", 'w')
        file.write('SUMMARY\n')

        for item in self.args:
            if type(item) is list:
                for element in item:
                    file.write(str(element) + '\n')
            else:
                file.write(str(item) + '\n')

        for item in args:
            if type(item) is list:
                for element in item:
                    file.write(str(element) + '\n')
            else:
                file.write(str(item) + '\n')

        logger.info("\nElapsed time: %.3f minutes = %.3f seconds" % (time / 60, time))
        file.write('\nElapsed time: %.3f minutes = %.3f seconds \n' % (time / 60, time))
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

        means = mcsamplefile.getMeans()

        stddev = mcsamplefile.std(self.paramList)

        summaryResults = []

        for i, param in enumerate(self.paramList):
            logger.info(self.paramList[i] + " : " + str(round(means[i], 4)) + \
                        "+/-" + str(round(stddev[i], 4)))
            summaryResults.append(self.paramList[i] + " : " + str(round(means[i], 4)) + \
                                  "+/-" + str(round(stddev[i], 4)))
        return summaryResults


class AllDerived:
    def __init__(self):
        # self.cpars = cpars
        self.Ol = Derivedparam('Ol', 0, '\Omega_\Lambda*')
        self.H0 = Derivedparam('H0', 0, 'H_0*')
        self.Age = Derivedparam('Age', 0, 'Age[Gyr]*')
        self.list = [self.Ol, self.H0, self.Age]

    def listDerived(self, like):
        self.like = like
        self.cpars = like.freeParameters()
        self.Ol.setValue(self.computeDerived('Ol'))
        self.H0.setValue(self.computeDerived('H0'))
        self.Age.setValue(self.computeDerived('Age'))
        return self.list

    def computeDerived(self, parname):
        import scipy.integrate as integrate
        if parname == 'Ol':
            for par in self.cpars:
                if par.name == 'Om':
                    return 1 - par.value
        elif parname == 'H0':
            for par in self.cpars:
                if par.name == 'h':
                    return par.value * 100
        elif parname == 'Age':
            return integrate.quad(self.compuAge, 0, 10 ** 5)[0] / 3.24076E-20 / (3.154E7 * 1.0E9)
        else:
            sys.exit('Define derived parameter', parname)

    def compuAge(self, z):
        return 1.0 / ((1 + z) * 100.0 * self.like.theory_.h * sp.sqrt(self.like.theory_.RHSquared_a(1.0 / (1 + z))))


class Derivedparam:
    def __init__(self, name, value, Ltxname=None):
        self.name = name
        if Ltxname:
            self.Ltxname = Ltxname
        else:
            self.Ltxname = name
        self.value = value

    def setLatexName(self, Ltx):
        self.Ltxname = Ltx

    def setValue(self, val):
        self.value = val