

from simplemc.analyzers.MaxLikeAnalyzer import MaxLikeAnalyzer
from simplemc.analyzers.MCMCAnalyzer import MCMCAnalyzer
from simplemc.runbase import ParseModel, ParseDataset
from simplemc.PostProcessing import PostProcessing
from simplemc.cosmo.Parameter import Parameter
from scipy.special import ndtri
import multiprocessing as mp
from simplemc import logger
import numpy as np
import sys, os
import time


#TODO #from simplemc.analyzers.SimpleGenetic import SimpleGenetic
#TODO #import corner
#TODO ##import emcee

class DriverMC:
    """
        This class is the manager and wrapper between all
        the analyzers and the pertinent functions.
    """
    def __init__(self, iniFile=None, **kwargs):
        """
        Read the input parameters or ini file
        Parameters
        ----------
        iniFile
        kwargs

        Returns
        -------

        """
        self.iniFile = iniFile
        if self.iniFile:
            self.iniReader(iniFile)
        else:
            self.chainsdir    = kwargs.pop('chainsdir', 'simplemc/chains')
            self.model        = kwargs.pop('model', None)
            self.prefact      = kwargs.pop('prefact', 'phy')
            self.varys8       = kwargs.pop('varys8', False)
            self.datasets     = kwargs.pop('datasets', 'HD')
            self.analyzername = kwargs.pop('analyzername', 'mcmc')
            self.addDerived   = kwargs.pop('addDerived', False)

            ## Next two are for custom model
            self.custom_parameters = kwargs.pop('custom_parameters', None)
            self.custom_function   = kwargs.pop('custom_function', None)
            ## Following two are for custom data
            self.path_to_data = kwargs.pop('path_to_data', None)
            self.path_to_cov  = kwargs.pop('path_to_cov', None)

            if os.path.exists(os.path.join(self.chainsdir)):
                self.chainsdir = os.path.join(self.chainsdir)
            else:
                logger.info("Your chains directory does not exist. Create a new one and try again.")
                sys.exit(1)
            if kwargs:
                logger.critical('Unexpected **kwargs for DriverMC: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'DriverMC **kwargs are:\n\tmodel\n\t'
                            'datasets\n\t'
                            'analyzername {"nested", "mcmc", "maxlike", "emcee" , "genetic"} Default: mcmc'
                            '\n\tchainsdir Default: SimpleMC_chains\n\t')
                sys.exit(1)


        #Initialize the Theory & Datasets
        T = ParseModel(self.model, custom_parameters=self.custom_parameters,
                                   custom_function=self.custom_function)
        L = ParseDataset(self.datasets, path_to_data=self.path_to_data,
                                        path_to_cov=self.path_to_cov)

        if self.prefact == "pre":  T.setVaryPrefactor()
        if self.varys8  == "True": T.setVarys8()
        T.printFreeParameters()

        #set the likelihood for a model
        L.setTheory(T)

        self.T, self.L = T, L

        pars_info = self.L.freeParameters()
        self.bounds     = [p.bounds for p in pars_info]
        self.means      = [p.value  for p in pars_info]
        self.paramsList = [p.name   for p in pars_info]
        self.dims       = len(self.paramsList)

        self.outputname = "{}_{}_{}_{}".format(self.model, self.prefact,
                                    self.datasets, self.analyzername)
        self.outputpath = "{}/{}".format(self.chainsdir, self.outputname)




    def executer(self, **kwargs):
        ti = time.time()
        if self.analyzername == 'mcmc':
            self.result = self.mcmcRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'nested':
            self.result = self.nestedRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'maxlike':
            self.result = self.maxLikeRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'emcee':
            self.result = self.emceeRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'genetic':
            self.result = self.geneticRunner(iniFile=self.iniFile, **kwargs)
        else:
            sys.exit("Sampler/Analyzer name invalid")
        self.ttime = time.time() - ti
        return self.outputname, self.T.freeParameters()


######################################Initialization#############################

    def iniReader(self, iniFile):
        """
         It reads the ini file.
        Parameters
            ini. file
        ----------
        iniFile

        Returns
        -------
            settings
        """
        import configparser
        self.config = configparser.ConfigParser()

        self.config.read(iniFile)
        self.chainsdir    = self.config.get('custom', 'chainsdir',
                                         fallback=os.path.join('simplemc/chains'))
        self.model        = self.config.get('custom', 'model')
        self.prefact      = self.config.get('custom', 'prefact',      fallback='phy')
        self.varys8       = self.config.get('custom', 'varys8',       fallback=False)
        self.datasets     = self.config.get('custom', 'datasets',     fallback='HD')
        self.analyzername = self.config.get('custom', 'analyzername', fallback='mcmc')
        self.addDerived   = self.config.get('custom', 'addDerived',   fallback=False)

        self.custom_parameters = self.config.get('custom', 'custom_parameters', fallback=None)
        self.custom_function   = self.config.get('custom', 'custom_function',   fallback=None)
        ## Following two are for custom data
        self.path_to_data = self.config.get('custom', 'path_to_data', fallback=None)
        self.path_to_cov  = self.config.get('custom', 'path_to_cov',  fallback=None)
        return True







    def mcmcRunner(self, iniFile=None, **kwargs):
        """
            This method calls MCMCAnalyzer.
            Returns [MCMCAnalyzer object, time, evidence via MCEvidence (if it is possible)]

        """
        if iniFile:
            nsamp = self.config.getint('mcmc', 'nsamp', fallback=15000)
            skip = self.config.getint('mcmc', 'skip', fallback=0)
            ## temperature at which to sample, weights get readjusted on the fly
            temp = self.config.getint('mcmc', 'temp', fallback=2)
            chainno = self.config.getint('mcmc', 'chainno', fallback=1)
            # self.GRcriteria = float(config['mcmc']['GRcriteria'])
            self.addDerived = self.config.getboolean('mcmc', 'addDerived', fallback=False)
            evidence = self.config.getboolean('mcmc', 'addDerived', fallback=False)
        else:
            nsamp = kwargs.pop('nsamp', 15000)
            skip = kwargs.pop('skip', 0)
            temp = kwargs.pop('temp', 2)
            chainno = kwargs.pop('chainno', 1)
            self.addDerived = kwargs.pop('addDerived', False)
            evidence = kwargs.pop('evidence', False)
            if kwargs:
                logger.critical('Unexpected **kwargs for MCMC: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'MCMC executer kwargs are:\n\tnsamp (int) Default: 15000\n\t'
                            'skip (int) Default 1000\n\ttemp (float) Default: 2.0'
                            '\n\tchainno (int) Default: 1\n\t'
                            'addDerived (boolean) Default: False\n\tevidence (boolean) Default: False')
                sys.exit(1)
                #raise TypeError('Unexpected **kwargs: {}'.format(kwargs))
        logger.info("\n\tnsamp: {}\n\tskip: {}\n\t"
                    "temp: {}\n\tchain num: {}\n\tevidence: {}".format(
                    nsamp, skip, temp, chainno, evidence))
        self.outputChecker()
        M = MCMCAnalyzer(self.L, self.outputpath,
                        skip=skip, nsamp=nsamp, temp = temp,
                        chain_num=chainno, addDerived=self.addDerived)
        if evidence:
            try:
                from MCEvidence import MCEvidence
                logger.info("Aproximating bayesian evidence with MCEvidence (arXiv:1704.03472)\n")
                MLE = MCEvidence(self.outputpath + ".txt" ).evidence()
                return ['mcmc', M, "Evidence with MCEvidence : {}\n".format(MLE)]
            except:
                #writeSummary(self.chainsdir, outputname, ttime)
                # print("Warning!")
                # print("MCEvidence could not calculate the Bayesian evidence [very small weights]\n")
                logger.error("MCEvidence could not calculate the Bayesian evidence [very small weights]")
        else:
            return ['mcmc', M]




    def setPriors(self):
        """
        After setTheory, you can set the priors.

        Returns a list of lists of the object of the Parameter class.
        """
        self.parameters = self.T.freeParameters()
        names      = []
        values     = []
        errorlist  = []
        boundlist  = []
        latexnames = []
        for parameter in self.parameters:
            names.append(parameter.name)
            values.append(parameter.value)
            errorlist.append(parameter.error)
            boundlist.append(parameter.bounds)
            latexnames.append(parameter.Ltxname)
        return [names, values, errorlist, boundlist, latexnames]


#################################### logLike and prior Transform function ###########################

    def instantiatePars(self, value):
        """
        This method returns an instance of
        Parameter objects with the sampler values

        -- value : it is the generated vector by the sampler
        """
        aux=[]
        for item in value:
            aux.append(item)
        names, values, errorlist, boundlist, latexnames = self.setPriors()
        instances = []
        for i in range(len(names)):
            instances.append(Parameter(names[i],
                aux[i], errorlist[i], boundlist[i], latexnames[i]))
        return instances

    def logLike(self, values):
        """
            If the sampler used isn't the MCMC of MCMCAnalyzer
            then, we need to set other types of likelihoods and
            priors objects. This method allows that. Therefore, it is a
            likelihood defined for an external samplers and it is used
            as parameter of the sampler run function.

            Parameters:

            values:     implicit values, they are generated by the sampler.

        """
        listPars = self.instantiatePars(values)
        self.T.updateParams(listPars)
        self.L.setTheory(self.T)
        if (self.L.name()=="Composite"):
            cloglikes=self.L.compositeLogLikes_wprior()
            loglike=cloglikes.sum()
        else:
            loglike = self.L.loglike_wprior()
        return loglike

    #priorsTransform
    def priorTransform(self, theta):
        """prior Transform for gaussian and flat priors"""
        priors = []
        n = self.nsigma

        if self.priortype == 'g':
            for c, bound in enumerate(self.bounds):
                mu = self.means[c]
                sigma = (bound[1]-bound[0])/n
                priors.append(mu+sigma*(ndtri(theta[c])))
        else:
            for c, bound in enumerate(self.bounds):
            # When theta 0-> append bound[0], if theta 1-> append bound[1]
                priors.append(theta[c]*(bound[1]-bound[0])+bound[0])
                # At this moment, np.array(priors) has shape (dims,)
        return np.array(priors)

###########loglike for genetic
    def logLikeGenetic(self, *v):
        values = []
        # print("Parameters:")
        for i, element in enumerate(v):
            values.append(element)
        listPars = self.instantiatePars(values)
        self.T.updateParams(listPars)
        self.L.setTheory(self.T)
        if (self.L.name() == "Composite"):
            cloglikes = self.L.compositeLogLikes_wprior()
            loglike = cloglikes.sum()
        else:
            loglike = self.L.loglike_wprior()

        return loglike

############# for emcee: logPosterior and logPrior
    def logPosterior(self, theta):
        """
        The natural logarithm of the joint posterior.

        Args:
            theta (tuple): a sample containing individual parameter values
            data (list): the set of data/observations
            sigma (float): the standard deviation of the data points
            x (list): the abscissa values at which the data/model is defined
        """

        lp = self.logPrior(theta)  # get the prior

        # if the prior is not finite return a probability of zero (log probability of -inf)
        if not np.isfinite(lp):
            return -np.inf

        # return the likeihood times the prior (log likelihood plus the log prior)
        return lp + self.logLike(theta)

    def logPrior(self, theta):
        """
        The natural logarithm of the prior probability.

        Args:
            theta (tuple): a sample containing individual parameter values

        Note:
            We can ignore the normalisations of the prior here.
            Uniform prior on all the parameters.
        """

        # set prior to 1 (log prior to 0) if in the range and zero (-inf) outside the range

        for i, bound in enumerate(self.bounds):
            if bound[0] < theta[i] < bound[1]:
                flag = True
            else:
                flag = False
                break

        if flag == True:
            return 0.0
        else:
            return -np.inf
#################################################Samplers ################################



    def nestedRunner(self, iniFile=None, **kwargs):
        """
            This method calls Dynesty samplers:
            --
            -- MULTINEST
            bound : {'none', 'single', 'multi', 'balls', 'cubes'},
        """
        if iniFile:
            nlivepoints = self.config.getint('nested', 'nlivepoints', fallback=1024)
            accuracy = self.config.getfloat('nested', 'accuracy', fallback=0.01)
            self.priortype = self.config.get('nested', 'priortype', fallback='u')
             #nsigma is the default value for sigma in gaussian priors
            self.nsigma  = self.config.get('nested', 'sigma', fallback=2)
            nestedType = self.config.get('nested', 'nestedType', fallback='multi')
            neuralNetwork = self.config.getboolean('nested', 'neuralNetwork', fallback=False)
            dynamic = self.config.getboolean('nested', 'dynamic', fallback=False)
            nproc = self.config.getint('nested', 'nproc', fallback=1)
            self.engine = self.config.get('nested', 'engine', fallback='dynesty')
            self.addDerived = self.config.getboolean('nested', 'addDerived', fallback=False)
            showfiles = self.config.getboolean('nested', 'showfiles', fallback=True)
        else:
            nlivepoints = kwargs.pop('nlivepoints', 1024)
            accuracy = kwargs.pop('accuracy', 0.01)
            self.priortype = kwargs.pop('priortype', 'u')
            nestedType = kwargs.pop('nestedType', 'multi')
            neuralNetwork = kwargs.pop('neuralNetwork', False)
            dynamic = kwargs.pop('dynamic', False)
            nproc = kwargs.pop('nproc', 1)
            self.engine = kwargs.pop('engine', 'dynesty')
            self.addDerived = kwargs.pop('addDerived', False)
            showfiles = kwargs.pop('showfiles', True)

            if kwargs:
                logger.critical('Unexpected **kwargs for nested sampler: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'Nested executer options are:\n\tnlivepoints (int) Default: 1024\n\t'
                            'accuracy (float) Default: 0.01\n\tpriortype ({"u", "g"}) Default: "u"\n\t'
                            'nestedType {"multi", "single", "balls", "cubes"} Default: "multi"\n\t'
                            'neuralNetwork (boolean) Default: True\n\t'
                            'dynamic (boolean) Default: False\n\t'
                            'addDerived (boolean) Default: True\n\t'
                            'engine {"dynesty", "nestle"} Default: "dynesty"\n\t'
                            'addDerived (boolean). Default: False\n\t'
                            'showfiles (boolean). Default: True.')
                sys.exit(1)

        self.outputpath += '_{}_{}'.format(self.engine, nestedType)
        if neuralNetwork:
            self.outputpath = "{}+ANN_".format(self.outputpath)
        self.outputChecker()

        pool, nprocess = self.mppool(nproc)
        logger.info("\n\tnlivepoints: {}\n"
                    "\taccuracy: {}\n"
                    "\tnested type: {}\n"
                    "\tengine: {}".format(nlivepoints, accuracy,
                                          nestedType, self.engine))
        if self.engine == 'dynesty':
            if showfiles:
                from simplemc.analyzers.dynesty import dynesty
            else:
                try:
                    import dynesty
                except:
                    from simplemc.analyzers.dynesty import dynesty

        if neuralNetwork:
            logger.info("\tUsing neural network.")
            #from bambi import run_pyBAMBI
            split = self.config.getfloat('neural', 'split', fallback=0.8)
            numNeurons = self.config.getint('neural', 'numNeurons', fallback=100)
            from simplemc.analyzers.pybambi.pybambimanager import BambiManager
            learner = 'keras'
            kerasmodel = None
            thumper = BambiManager(self.logLike, learner, proxy_tolerance=0.1,
                                   failure_tolerance=0.5, ntrain=nlivepoints, model=kerasmodel)
            self.logLike = thumper.loglikelihood
            # M = run_pyBAMBI(self.logLike, self.priorTransform, self.dims,
            #                 eff=accuracy, nestedType=nestedType, dynamic=dynamic,
            #                 nlive=nlivepoints, learner='keras')
        if dynamic:
            logger.info("\nUsing dynamic nested sampling...")
            sampler = dynesty.DynamicNestedSampler(self.logLike, self.priorTransform,
                      self.dims, bound=nestedType, pool=pool, queue_size=nprocess)

            sampler.run_nested(nlive_init=nlivepoints, dlogz_init=0.05, nlive_batch=100,
                            maxiter_init=10000, maxiter_batch=1000, maxbatch=10,
                               outputname=self.outputpath,
                               addDerived=self.addDerived, smcloglike=self.L)

            M = sampler.results

        elif self.engine == 'dynesty':
            sampler = dynesty.NestedSampler(self.logLike, self.priorTransform, self.dims,
                        bound=nestedType, sample = 'unif', nlive = nlivepoints,
                        pool = pool, queue_size=nprocess)
            sampler.run_nested(dlogz=accuracy, outputname=self.outputpath,
                               addDerived=self.addDerived, smcloglike=self.L)
            M = sampler.results


        elif self.engine == 'nestle':
            import nestle
            M = nestle.sample(self.logLike, self.priorTransform, ndim=self.dims, method=nestedType,
            npoints=nlivepoints, dlogz=accuracy, callback=nestle.print_progress,
            pool = pool, queue_size = nprocess)
        else:
            logger.critical('wrong selection')
            sys.exit(1)
        try:
            pool.close()
        except:
            pass

        return ['nested', M, M.summary(), 'nested :{}'.format(nestedType),
                'dynamic : {}'.format(dynamic), 'ANN :{}'.format(neuralNetwork)]


    def emceeRunner(self, iniFile=None, **kwargs):
        if iniFile:
            walkers = self.config.getint('emcee', 'walkers', fallback=500)
            nsamp = self.config.getint('emcee', 'nsamp', fallback=20000)
            burnin = self.config.getint('emcee', 'burnin', fallback=0)
            nproc = self.config.getint('emcee', 'nproc', fallback=1)
        else:
            walkers = kwargs.pop('walkers', 30)
            nsamp = kwargs.pop('nsamp', 20000)
            burnin = kwargs.pop('burnin',0)
            nproc = kwargs.pop('nproc', 1)
            if kwargs:
                logger.critical('Unexpected **kwargs for emcee sampler: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'Emcee executer options are:'
                            '\n\twalkers (int) Default: 30\n\t'
                            'nsamp (int) Default: 20000\n\t'
                            'burnin (int) Default: 0\n\t'
                            'nproc (int) Default: 1')
                sys.exit(1)
        logger.info("\n\twalkers: {}\n\tnsamp: {}\n"
                    "\tburnin: {}\n\t"
                    "nproc: {}".format(walkers, nsamp, burnin, nproc))
        self.outputpath = "{}_{}_walkers".format(self.outputpath, walkers)
        self.outputChecker()
        pool, _ = self.mppool(nproc)
        # initial_state = None
        ini = []
        for bound in self.bounds:
            ini.append(np.random.uniform(bound[0], bound[1], walkers))
        inisamples = np.array(ini).T # initial samples

        # set up the sampler
        sampler = emcee.EnsembleSampler(walkers, self.dims, self.logPosterior, pool=pool)

        # pass the initial samples and total number of samples required
        sampler.run_mcmc(inisamples, nsamp+burnin, progress=True)

        # extract the samples (removing the burn-in)
        postsamples = sampler.chain[:, burnin:, :].reshape((-1, self.dims))
        self.burnin = burnin

        try:
            pool.close()
        except:
            pass
        # fig = corner.corner(postsamples)
        # fig.savefig('emcee.png')

        return ['emcee', sampler, 'walkers : {}'.format(walkers), 'samples: {}'.format(nsamp)]
###################################Optimizers##############################################

    def maxLikeRunner(self, iniFile=None, **kwargs):
        self.outputpath += 'optimization'
        self.outputChecker()
        if iniFile:
            withErrors = self.config.getboolean('MaxLike', 'withErrors', fallback=False)
        else:
            withErrors = kwargs.pop('withErrors', False)
            if kwargs:
                logger.critical('Unexpected **kwargs for MaxLike: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'MaxLikeAnalyzer executer options are:'
                            '\n\twithErrors (boolean) Default: False')
                sys.exit(1)

        A = MaxLikeAnalyzer(self.L, withErrors=withErrors)
        self.T.printParameters(A.params)
        return A

    def geneticRunner(self, iniFile=None, **kwargs):
        self.outputpath = '{}_optimization'.format(self.outputpath)
        self.outputChecker()

        if iniFile:
            n_individuals = self.config.getint('genetic', 'n_individuals', fallback=100)
            n_generations = self.config.getint('genetic', 'n_generations' , fallback=1000)
            selection_method = self.config.get('genetic', 'selection_method', fallback='tournament')
            mut_prob = self.config.getfloat('genetic', 'mut_prob', fallback=0.4)
        else:
            n_individuals = kwargs.pop('n_individuals', 100)
            n_generations = kwargs.pop('n_generations', 1000)
            selection_method = kwargs.pop('selection_method', 'tournament')
            mut_prob = kwargs.pop('mut_prob', 0.4)
            if kwargs:
                logger.critical('Unexpected **kwargs for genetic optimizer: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'genetic executer options are:'
                            '\n\tn_individuals (int) Default: 100\n\t'
                            'n_generations (int) Default: 1000\n\t'
                            'selection_method {"tournament","rank","roulette"} Default: "tournament"'
                            '\n\tmut_prob (float) Default: 0.4')
                sys.exit(1)
        logger.info("\n\tn_individuals: {}\n\tn_generations_ {}"
                    "\n\tselection method: {}\n\t"
                    "mut prob: {}".format(n_individuals, n_generations,
                                        selection_method,mut_prob))
        M = SimpleGenetic(self.logLikeGenetic, self.dims, self.bounds,
                          n_individuals=n_individuals,
                          n_generations=n_generations,
                          prob_mut=mut_prob, method_selection=selection_method,
                          outputname=self.outputpath)
        return ['genetic', M]


###############################Post-processing############################################
    def outputChecker(self):
        while True:
            if os.path.isfile(self.outputpath+".txt"):
                logger.info("{0} file already exists, {0}_new was created".format(self.outputpath))
                self.outputpath = "{}_new".format(self.outputpath)
            else:
                break
        while True:
            flag = False
            for i in range(1,10):
                if os.path.isfile("{}_{}.txt".format(self.outputpath, i)):
                    flag = True
            if flag:
                logger.info("{0}_{1} file already exists, {0}_new was created".format(self.outputpath, i))
                self.outputpath = "{}_new".format(self.outputpath)
            else:
                break
        if self.analyzername == 'nested':
            while True:
                if os.path.isfile(self.outputpath + "_live.txt") or os.path.isfile(self.outputpath + "_dead.txt"):
                    logger.info("{0} file already exists, {0}_new was created".format(self.outputpath))
                    self.outputpath = "{}_new".format(self.outputpath)
                else:
                    break

        return True

    def postprocess(self, summary=True, stats=False, addtxt=None):
        if addtxt:
            self.result.extend(addtxt)
        if self.analyzername == 'nested':
            pp = PostProcessing(self.result, self.paramsList, self.outputpath,
                engine=self.engine, addDerived=self.addDerived, loglike=self.L)
            pp.paramFiles()
            pp.saveNestedChain()
        elif self.analyzername == 'emcee':
            pp = PostProcessing(self.result, self.paramsList, self.outputpath,
                                skip=self.burnin, addDerived=self.addDerived, loglike=self.L)
            pp.paramFiles()
            pp.saveEmceeSamples()
        else:
            pp = PostProcessing(self.result, self.paramsList, self.outputpath)
        if summary:
            pp.writeSummary(self.ttime)

    #def plotter(self):
    #    from PlotterMC import PlotterMC
    #    plot = PlotterMC(self.dims, chainsdir = self.chainsdir, chainsfile = self.outputname)

# ### pool from multiprocessing

    def mppool(self, nproc):
        if nproc <= 0:
            ncores = mp.cpu_count()
            nprocess = ncores//2
            logger.info("Using  {} processors of {}.".format(nprocess, ncores))
            pool = mp.Pool(processes=nprocess)
        elif nproc == 1:
            logger.info("Using 1 processor")
            nprocess = None
            pool = None
        else:
            nprocess = nproc
            logger.info("Using {} processors.".format(nprocess))
            pool = mp.Pool(processes=nprocess)

        return pool, nprocess
