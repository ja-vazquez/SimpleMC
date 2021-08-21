
#TODO check: if self.analyzername is None

from .analyzers import MaxLikeAnalyzer
from .analyzers import GA_deap
from .analyzers import MCMCAnalyzer
from .analyzers import DynamicNestedSampler, NestedSampler
from .analyzers import EnsembleSampler
from .cosmo.Derivedparam import AllDerived
from . import ParseDataset, ParseModel
from . import PostProcessing
from scipy.special import ndtri
from simplemc import logger
import numpy as np
import sys, os
import time


class DriverMC:
    """
        This class is the manager and wrapper between all
        the analyzers and the pertinent functions.
        It reads the input parameters or ini file.

        Parameters
        -----------
        iniFile: ini file
            Text file with ini extension that contains all the settings
            to SimpleMC. If use this option the following kwargs not are necessary.

        chainsdir : str
            Directory for the outputs.
        model : str
            Choose the model {LCDOM, LCDMasslessnu, nuLCDM, NeffLCDM, noradLCDM, nuoLCDM, nuwLCDM, oLCDM, wCDM, waCDM, owCDM,"\
            owaCDM, JordiCDM, WeirdCDM, TLight, StepCDM, Spline, PolyCDM, fPolyCDM, Decay, Decay01, Decay05,"\
            EarlyDE, EarlyDE_rd_DE, SlowRDE}

        prefact :str
            {phy, pre}

        vary8 : bool
            Default False

        datasets str:
            Default HD (Hubble distance, i. e. Cosmic Chronometers).
            You can combine HD+SN+BBAO+Planck+UnionSN+...

        analyzername : str
            The name of the analyzer. It can be a sampler: {mcmc, nested, emcee}
            or a optimizer: {maxlike, ga_deap}

        compute_derived : bool
            True generates at the flight some derived parameters (such as
            Omega_Lambda or Universe Age, and save them in the output text file.

        custom_parameters : list
            List of Parameter instances.

        custom_function : method
            Custom method that reads a parameter list and a vector x, unzip the list,
            and return a f(x) in terms of the parameters.

        path_to_data : str
            path of a dataset text file.

        path_to_cov : str
            path of a covariance matrix text file.

        fn : str
            Type of function to use in the likelihood due a custom data {"generic", "hz", ...}.
    """
    def __init__(self, iniFile=None, **kwargs):

        self.iniFile = iniFile
        if self.iniFile:    self.iniReader(iniFile)
        else:
            self.chainsdir    = kwargs.pop('chainsdir', 'simplemc/chains')
            self.model        = kwargs.pop('model',   None)
            self.prefact      = kwargs.pop('prefact', 'phy')
            self.varys8       = kwargs.pop('varys8',  False)
            self.datasets     = kwargs.pop('datasets','HD')
            self.analyzername = kwargs.pop('analyzername', None)
            self.addDerived   = kwargs.pop('addDerived', False)
            self.useNeuralLike = kwargs.pop('useNeuralLike', False)
            self.mcevidence = kwargs.pop('mcevidence', False)
            self.overwrite = kwargs.pop('overwrite', True)


            ## Next two are for custom model
            self.custom_parameters = kwargs.pop('custom_parameters', None)
            self.custom_function   = kwargs.pop('custom_function', None)
            ## Following two are for custom data
            self.path_to_data = kwargs.pop('path_to_data', None)
            self.path_to_cov  = kwargs.pop('path_to_cov', None)
            self.fn = kwargs.pop("fn", "generic")

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
                            'analyzername {"nested", "mcmc", "maxlike", "emcee" , "ga_deap"} Default: mcmc'
                            '\n\tchainsdir Default: SimpleMC_chains\n\t')
                sys.exit(1)


        #Initialize the Theory & Datasets
        T = ParseModel(self.model, custom_parameters=self.custom_parameters,
                                   custom_function=self.custom_function)

        L = ParseDataset(self.datasets, path_to_data=self.path_to_data,
                                        path_to_cov=self.path_to_cov, fn=self.fn)

        if self.prefact == "pre":  T.setVaryPrefactor()
        if self.varys8  == True: T.setVarys8()
        T.printFreeParameters()

        #set the likelihood for a model
        L.setTheory(T)

        self.T, self.L = T, L

        self.pars_info  = self.L.freeParameters()
        self.bounds     = [p.bounds for p in self.pars_info]
        self.means      = [p.value  for p in self.pars_info]
        self.paramsList = [p.name   for p in self.pars_info]
        self.dims       = len(self.paramsList)
        self.result     = None

        self.root = "{}_{}_{}".format(self.model, self.prefact,
                                    self.datasets)
        self.outputpath = "{}/{}".format(self.chainsdir, self.root)

        if self.useNeuralLike:
            neural_model = self.neuralLike(iniFile=self.iniFile)
            self.logLike = neural_model.loglikelihood



    def executer(self, **kwargs):
        """
        This is a wrapper of the runners of the analyzer in order to make
        easier the execution, mainly if is through an ini file.
        **kwargs from mcmcRunner, nestedRunner, emceeRunner, genetic_deap and maxlikeRunner.

        """
        if self.analyzername == 'mcmc':
            self.mcmcRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'nested':
            self.nestedRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'emcee':
            self.emceeRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'maxlike':
            self.maxLikeRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'ga_deap':
            self.geneticdeap(iniFile=self.iniFile, **kwargs)
        else:
            sys.exit("{}: Sampler/Analyzer name invalid".format(self.analyzername))
        self.postprocess()
        return True


##----------------------Initialization ------------------------

    def iniReader(self, iniFile):
        """
        It reads the ini file through configparser.

        Parameters
        -----------
        iniFile : file .ini
            Text file with settings

        """
        import configparser
        self.config = configparser.ConfigParser()

        self.config.read(iniFile)
        self.chainsdir    = self.config.get(        'custom', 'chainsdir',\
                                         fallback=os.path.join('simplemc/chains'))
        self.model        = self.config.get(        'custom', 'model')
        self.prefact      = self.config.get(        'custom', 'prefact',      fallback='phy')
        self.datasets     = self.config.get(        'custom', 'datasets',     fallback='HD')
        self.analyzername = self.config.get(        'custom', 'analyzername', fallback=None)
        self.varys8       = self.config.getboolean( 'custom', 'varys8',       fallback=False)
        self.addDerived   = self.config.getboolean( 'custom', 'addDerived',   fallback=False)
        self.useNeuralLike = self.config.getboolean('custom', 'useNeuralLike', fallback=False)
        self.mcevidence = self.config.getboolean('custom', 'mcevidence', fallback=False)
        self.overwrite = self.config.getboolean('custom', 'overwrite', fallback=True)

        self.custom_parameters = self.config.get(   'custom', 'custom_parameters', fallback=None)
        self.custom_function   = self.config.get(   'custom', 'custom_function',   fallback=None)
        ## Following two are for custom data
        self.path_to_data = self.config.get(        'custom', 'path_to_data', fallback=None)
        self.path_to_cov  = self.config.get(        'custom', 'path_to_cov',  fallback=None)
        self.fn = self.config.get('custom', 'fn',  fallback="generic")
        return True




##======== ======== ======== Samplers ======== ======== ========


##---------------------- MCMC ----------------------

    def mcmcRunner(self, iniFile=None, **kwargs):
        """
        This method calls MCMCAnalyzer.

        Parameters
        ------------
        nsamp : int
            Number of mcmc steps.

        skip : int
            Burn-in.

        temp : float
            Temperature for the weights.

        GRstop : float
            Gelman Rubin criteria for stopping (0, 0.1].

        evidence : bool
            True if after the mcmc chain was generated,
            estimates bayesian evidence throug MCEvidence (arXiv:1704.03472).

        """
        if iniFile:
            nsamp    = self.config.getint(      'mcmc', 'nsamp',   fallback=50000)
            skip     = self.config.getint(      'mcmc', 'skip',    fallback=300)
            ## temperature at which to sample, weights get readjusted on the fly
            temp     = self.config.getfloat(    'mcmc', 'temp',    fallback=2)
            GRstop   = self.config.getfloat(    'mcmc', 'GRstop',  fallback=0.01)
            checkGR  = self.config.getfloat(    'mcmc', 'checkGR', fallback=500)
            evidence = self.config.getboolean(  'mcmc', 'evidence',fallback=False)
        else:
            nsamp    = kwargs.pop('nsamp', 50000)
            skip     = kwargs.pop('skip',  300)
            temp     = kwargs.pop('temp',  2)
            GRstop   = kwargs.pop('GRstop', 0.01)
            checkGR  = kwargs.pop('checkGR', 500)
            evidence = kwargs.pop('evidence', False)
            if kwargs:
                logger.critical('Unexpected **kwargs for MCMC: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use default values.\n'
                            'MCMC executer kwargs are:\n\tnsamp (int) Default: 50000\n\t'
                            'skip (int) Default 300\n\ttemp (float) Default: 2.0'
                            '\n\tevidence (bool) Default: False')
                sys.exit(1)
                #raise TypeError('Unexpected **kwargs: {}'.format(kwargs))
        logger.info("\n\tnsamp: {}\n\tskip: {}\n\t"
                    "temp: {}\n\tevidence: {}".format(
                    nsamp, skip, temp, evidence))
        if self.analyzername is None: self.analyzername = 'mcmc'
        self.outputpath = "{}_{}".format(self.outputpath, self.analyzername)
        #Check whether the file already exists
        self.outputChecker()
        ti = time.time()

        #Main process
        M = MCMCAnalyzer(self.L, self.outputpath, skip=skip, nsamp=nsamp, temp = temp,
                         addDerived=self.addDerived, GRstop=GRstop, checkGR=checkGR)

        self.ttime = time.time() - ti

        self.result = ['mcmc', M.get_results()[:2], "Maxlike: {}".format(M.maxloglike),
                       "Gelman-Rubin diagnostic: {}".format(M.get_results()[2])]

        return True



##---------------------- Nested ----------------------

    def nestedRunner(self, iniFile=None, **kwargs):
        """
        This method calls Dynesty nested samplers.

        Parameters
        ___________
        dynamic : bool
            Default `False`

        neuralNetwork : bool
            If True use a pybambi neural network.
            Default: False.

        nestedType : str
            {single, multi, balls, cubes}

        nlivepoints : int
            Number of live points.

        accuracy : float
            Stopping criteria in terms of logz.

        nproc : int
            Number of processors to parallelize.
            Use 1 or 0 if you don't want parallelise.

        priortype : str
            Gaussian or uniform prior {'g', 'u'}.

        nsigma : float
            Sigma for gaussian priors.


        """
        if iniFile:
            dynamic     = self.config.getboolean(   'nested', 'dynamic',      fallback=False)
            neuralNetwork = self.config.getboolean( 'nested', 'neuralNetwork',fallback=False)
            nestedType  = self.config.get(          'nested', 'nestedType',   fallback='multi')
            nlivepoints = self.config.getint(       'nested', 'nlivepoints',  fallback=1024)
            accuracy    = self.config.getfloat(     'nested', 'accuracy',     fallback=0.01)
            nproc       = self.config.getint(       'nested', 'nproc',        fallback=1)

            self.priortype = self.config.get('nested', 'priortype', fallback='u')
            #nsigma is the default value for sigma in gaussian priors
            self.nsigma = self.config.get('nested', 'sigma', fallback=2)

            # Neural network settings
            split      = self.config.getfloat('neural', 'split', fallback=0.8)
            numNeurons = self.config.getint('neural', 'numNeurons', fallback=100)
            epochs = self.config.getint('neural', 'epochs', fallback=100)
            model  = self.config.get( 'model', 'model',   fallback=None)
            savedmodelpath = self.config.get('neural', 'savedmodelpath', fallback=None)
            it_to_start_net = self.config.getint('neural', 'it_to_start_net', fallback=None)
            dlogz_start = self.config.getfloat('neural', 'dlogz_start', fallback=5)
            updInt = self.config.getint('neural', 'updInt', fallback=nlivepoints)
            proxy_tolerance = self.config.getfloat('neural', 'proxy_tolerance', fallback=0.3)
            failure_tolerance = self.config.getfloat('neural', 'failure_tolerance', fallback=0.5)

        else:
            dynamic     = kwargs.pop('dynamic',    False)
            neuralNetwork = kwargs.pop('neuralNetwork', False)
            nestedType  = kwargs.pop('nestedType', 'multi')
            nlivepoints = kwargs.pop('nlivepoints', 1024)
            accuracy    = kwargs.pop('accuracy',    0.01)
            nproc       = kwargs.pop('nproc', 1)

            self.priortype = kwargs.pop('priortype', 'u')
            self.nsigma = kwargs.pop('sigma', 2)

            # For neural networks
            split = kwargs.pop('split', 0.8)
            numNeurons = kwargs.pop('numNeurons', 100)
            epochs = kwargs.pop('epochs', 100)
            model = kwargs.pop('model', None)
            savedmodelpath = kwargs.pop('savedmodelpath', None)
            it_to_start_net = kwargs.pop('it_to_start_net', 10000)
            dlogz_start = kwargs.pop('dlogz_start', 5)
            updInt = kwargs.pop('updInt', nlivepoints)
            proxy_tolerance = kwargs.pop('proxy_tolerance', 0.3)
            failure_tolerance = kwargs.pop('failure_tolerance', 0.5)

            if kwargs:
                logger.critical('Unexpected **kwargs for nested sampler: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'Nested executer options are:\n\tnlivepoints (int) Default: 1024\n\t'
                            'accuracy (float) Default: 0.01\n\tpriortype ({"u", "g"}) Default: "u"\n\t'
                            'nestedType {"multi", "single", "balls", "cubes"} Default: "multi"\n\t'
                            'neuralNetwork (bool) Default: True\n\t'
                            'dynamic (bool) Default: False\n\t'
                            'addDerived (bool) Default: True')
                sys.exit(1)

        #stored output files
        if self.analyzername is None: self.analyzername = 'nested'
        self.outputpath = '{}_{}_{}'.format(self.outputpath, self.analyzername, nestedType)
        if neuralNetwork:
            self.outputpath = "{}_ANN".format(self.outputpath)
        self.outputChecker()

        self.neuralNetwork = neuralNetwork
        #paralel run
        pool, nprocess = self.mppool(nproc)
        logger.info("\n\tnlivepoints: {}\n"
                    "\taccuracy: {}\n"
                    "\tnested type: {}".format(nlivepoints, accuracy, nestedType))

        ti = time.time()
        if neuralNetwork:
            logger.info("\tUsing neural network.")
            from simplemc.analyzers.pybambi.bambi import bambi
            # self.logLike =
            thumper = bambi(self.logLike, self.dims,
                            split=split, numNeurons=numNeurons,
                            epochs=epochs, model=model,
                            savedmodelpath=savedmodelpath,
                            it_to_start_net=it_to_start_net,
                            updInt=updInt, dlogz_start=dlogz_start,
                            proxy_tolerance=proxy_tolerance,
                            failure_tolerance=failure_tolerance)

            self.logLike = thumper.loglikelihood
            dumper = thumper.dumper
        else:
            dumper = None


        if dynamic:
            logger.info("\nUsing dynamic nested sampling...")
            sampler = DynamicNestedSampler(self.logLike, self.priorTransform,
                                           self.dims, bound=nestedType, pool=pool,
                                           queue_size=nprocess)

            sampler.run_nested(nlive_init=nlivepoints, dlogz_init=0.05, nlive_batch=100,
                            maxiter_init=10000, maxiter_batch=1000, maxbatch=10,
                            outputname=self.outputpath, addDerived=self.addDerived, simpleLike=self.L)
            M = sampler.results


        else:
            sampler = NestedSampler(self.logLike, self.priorTransform, self.dims,
                        bound=nestedType, sample = 'unif', nlive = nlivepoints,
                        pool = pool, queue_size=nprocess, use_pool={'loglikelihood': False})
            sampler.run_nested(dlogz=accuracy, outputname=self.outputpath,
                               addDerived=self.addDerived, simpleLike=self.L, dumper=dumper)
            M = sampler.results

        try:
            pool.close()
        except:
            pass
        self.ttime = time.time() - ti
        self.result = ['nested', M, M.summary(), 'nested :{}'.format(nestedType),
                       'dynamic : {}'.format(dynamic), 'ANN :{}'.format(neuralNetwork)]
        return True


##---------------------- EMCEE ----------------------

    def emceeRunner(self, iniFile=None, **kwargs):
        """
        This method calls the emcee library to use ensamble sampler.

        Parameters
        -----------
        walkers : int
            Number of walkers or ensambles.

        nsamp : int
            Number of mcmc steps for each walker.

        burnin : int
            skip steps.

        nproc : int
            Number of processors in order to parallelise.

        """
        if iniFile:
            walkers = self.config.getint('emcee', 'walkers', fallback=self.dims*2+2)
            nsamp   = self.config.getint('emcee', 'nsamp', fallback=10000)
            burnin  = self.config.getint('emcee', 'burnin', fallback=0)
            nproc   = self.config.getint('emcee', 'nproc', fallback=1)
        else:
            walkers = kwargs.pop('walkers', self.dims*2+2)
            nsamp   = kwargs.pop('nsamp', 10000)
            burnin  = kwargs.pop('burnin', 0)
            nproc   = kwargs.pop('nproc', 1)
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

        if self.analyzername is None: self.analyzername = 'emcee'
        self.outputpath = "{}_{}_{}_walkers".format(self.outputpath, self.analyzername, walkers)
        self.outputChecker()
        pool, _ = self.mppool(nproc)
        # initial_state = None
        ini = []
        for bound in self.bounds:
            ini.append(np.random.uniform(bound[0], bound[1], walkers))
        inisamples = np.array(ini).T # initial samples

        ti = time.time()
        sampler = EnsembleSampler(walkers, self.dims,
                                  self.logPosterior, pool=pool)
        #testing
        sampler.sample(initial_state=self.means, tune=True, thin_by=3)
        # pass the initial samples and total number of samples required
        sampler.run_mcmc(inisamples, nsamp + burnin,
                         progress=True, outputname=self.outputpath,
                         addDerived=self.addDerived, simpleLike=self.L)
        self.ttime = time.time() - ti
        self.burnin = burnin
        try:
            pool.close()
        except:
            pass
        self.result = ['emcee', sampler, 'walkers : {}'.format(walkers), 'samples: {}'.format(nsamp)]
        return True



##======== ======== ======== Optimizers ======== ======== ========


##---------------------- MaxLikeAnalizer ----------------------

    def maxLikeRunner(self, iniFile=None, **kwargs):
        """
        It calls MaxLikeAnalyzer class.

        Parameters
        ----------
        withErrors : bool
            Plot errors.

        plot_par1 : bool
            First parameter to plot.

        plot_par2 : bool
            Second parameter to plot.

        """
        if self.analyzername is None:
            self.analyzername = 'maxlike'
        self.outputpath = '{}_{}_optimization'.format(self.outputpath, self.analyzername)
        self.outputChecker()
        if iniFile:
            compute_errors = self.config.getboolean('maxlike', 'compute_errors', fallback=False)
            show_contours = self.config.getboolean('maxlike', 'show_contours', fallback=False)
            plot_param1 = self.config.get('maxlike', 'plot_param1', fallback=None)
            plot_param2 = self.config.get('maxlike', 'plot_param2', fallback=None)
            compute_derived = self.config.getboolean('maxlike', 'compute_derived', fallback=False)
        else:
            compute_errors = kwargs.pop('compute_errors', False)
            show_contours = kwargs.pop('show_contours', False)
            plot_param1 = kwargs.pop('plot_param1', None)
            plot_param2 = kwargs.pop('plot_param2', None)
            compute_derived = kwargs.pop('compute_derived ', False)
            if kwargs:
                logger.critical('Unexpected **kwargs for MaxLike: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'MaxLikeAnalyzer executer options are:'
                            '\n\twithErrors (bool) Default: False')
                sys.exit(1)
        ti = time.time()
        A = MaxLikeAnalyzer(self.L, self.model, compute_errors=compute_errors,
                            compute_derived=compute_derived, show_contours=show_contours,\
                            plot_param1=plot_param1, plot_param2=plot_param2)
        params = self.T.printParameters(A.params)
        self.ttime = time.time() - ti
        self.result = ['maxlike', A, params, 'Optimal loglike: {}'.format(A.opt_loglike)]
        return True


##---------------------- Genetic Algorithms ----------------------

    def geneticdeap(self, iniFile=None, **kwargs):
        """
        Genetic algorithms from Deap library.

        """
        if self.analyzername is None: self.analyzername = 'ga_deap'
        self.outputpath = '{}_{}'.format(self.outputpath, self.analyzername)
        self.outputChecker()
        if iniFile:
            plot_fitness = self.config.getboolean('ga_deap', 'plot_fitness', fallback=False)
            compute_errors = self.config.getboolean('ga_deap', 'compute_errors', fallback=False)
            show_contours = self.config.getboolean('ga_deap', 'show_contours', fallback=False)
            plot_param1 = self.config.get('ga_deap', 'plot_param1', fallback=None)
            plot_param2 = self.config.get('ga_deap', 'plot_param2', fallback=None)

            population = self.config.getint('ga_deap','population', fallback=20)
            crossover = self.config.getfloat('ga_deap', 'crossover', fallback=0.7)
            mutation = self.config.getfloat('ga_deap', 'mutation', fallback=0.3)
            max_generation = self.config.getint('ga_deap', 'max_generation', fallback=100)
            hof_size = self.config.getint('ga_deap','hof_size', fallback=1)
            crowding_factor = self.config.getfloat('ga_deap', 'crowding_factor',fallback=1)
        else:
            plot_fitness = kwargs.pop('plot_fitness', False)
            compute_errors = kwargs.pop('compute_errors', False)
            show_contours = kwargs.pop('show_contours', False)
            plot_param1 = kwargs.pop('plot_param1', None)
            plot_param2 = kwargs.pop('plot_param2', None)

            population = kwargs.pop('population', 20)
            crossover = kwargs.pop('crossover', 0.7)
            mutation = kwargs.pop('mutation', 0.3)
            max_generation = kwargs.pop('max_generation', 100)
            hof_size = kwargs.pop('hof_size', 1)
            crowding_factor = skwargs.pop('crowding_factor', 1)

        ti = time.time()
        M = GA_deap(self.L, self.model, outputname=self.outputpath,
                    population=population, crossover=crossover,
                    mutation=mutation, max_generation=max_generation,
                    hof_size=hof_size, crowding_factor=crowding_factor,
                    plot_fitness=plot_fitness, compute_errors=compute_errors,
                    show_contours=show_contours, plot_param1=plot_param1,
                    plot_param2=plot_param2)
        result = M.main()
        self.ttime = time.time() - ti
        #M.plotting()
        self.result = ['genetic', result, 'Population: {}'.format(population),
                       'Max number of generations: {}'.format(max_generation),
                       'Mutation: {}'.format(mutation), 'Crosover: {}'.format(crossover),
                       result[3]]
        return True

##---------------------- logLike and prior Transform function ----------------------
##---------------------- for nested samplers ----------------------

    def logLike(self, values):
        """
        If the sampler used isn't the MCMC of MCMCAnalyzer then, we need to set
        other types of likelihoods and priors objects. This method allows that. It is a
        likelihood defined for an external samplers and is used
        as parameter of the sampler run function.

        Parameters
        -----------

        values : n-dim vector
            implicit values, they are generated by the sampler.

        """

        assert len(self.pars_info) == len(values)
        for pars, val in zip(self.pars_info, values):
            pars.setValue(val)

        self.T.updateParams(self.pars_info)
        self.L.setTheory(self.T)
        if (self.L.name()=="Composite"):
            cloglikes=self.L.compositeLogLikes_wprior()
            loglike=cloglikes.sum()
        else:
            loglike = self.L.loglike_wprior()
        return loglike


    #priorsTransform
    def priorTransform(self, theta):
        """
        Prior Transform for gaussian and flat priors

        Parameters
        -----------

        theta : array
            Vector of the parameter space
        """
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

############# for emcee: logPosterior and logPrior
    def logPosterior(self, theta):
        """
        The natural logarithm of the joint posterior.

        Parameters
        ------------
        theta : tuple
            A sample containing individual parameter values

        data : list
            The set of data/observations

        sigma : float
            The standard deviation of the data points

        x : list
            The abscissa values at which the data/model is defined
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

        Parameters
        -----------
            theta : tuple
                A sample containing individual parameter values
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




###############################Post-processing############################################
    def outputChecker(self):
        """
        This method check if the name of the outputfile exists, if it already exists creates a
        new one with extension _new in its name.

        """
        self.paramFiles()
        i = 1
        if self.overwrite:
            self.outputpath = "{}".format(self.outputpath)
        else:
            if os.path.isfile("{}_1.txt".format(self.outputpath)):
                sys.exit('File with outputname {}_1.txt already exists.\n'
                         'Please move your files or set overwrite=True to g'
                         'overwrite outputs'.format(self.outputpath))
        return True

    def paramFiles(self):
        """
        This method writes the .paramnames file with theirs LaTeX names.

        Parameters
        -----------

        T : model
            T is result of ParseModel(model)
        L : likelihood
            L is result of ParseDataset(datasets)

        """
        cpars   = self.L.freeParameters()
        parfile = "{}.paramnames".format(self.outputpath)
        fpar = open(parfile, 'w')
        for p in cpars:
            fpar.write(p.name + "\t\t\t" + p.Ltxname + "\n")
        if self.addDerived:
            AD = AllDerived()
            for pd in AD.list:
                fpar.write(pd.name + "\t\t\t" + pd.Ltxname + "\n")
        if self.analyzername in ['mcmc', 'nested', 'emcee']:
            if (self.L.name() == "Composite"):
                self.sublikenames = self.L.compositeNames()
                for name in self.sublikenames:
                    fpar.write(name + "_like \t\t\t" + name + "\n")
                fpar.write("theory_prior \t\t\t None \n")


    def postprocess(self, addtxt=None):
        """
        It calls the PostProcessing class.

        Parameters
        ----------
         summary : bool
            True for save summary.

         addtxt : list
            A list with strings to save with the summary.
        """
        if addtxt:
            self.result.extend(addtxt)
        pp = PostProcessing(self.result, self.paramsList, self.outputpath,
                            addDerived=self.addDerived, loglike=self.L)
        if self.mcevidence:
            try:
                ev = pp.mcevidence()
                pp.writeSummary(self.ttime, ev)
            except:
                pp.writeSummary(self.ttime)
        else:
            pp.writeSummary(self.ttime)


    def plot(self, show=False):
        """
        Simple connection with the plotters.

        Parameters
        -----------
        show : bool
            Default False
        """
        from .tools.SimplePlotter import SimplePlotter
        figure = SimplePlotter(self.chainsdir, self.paramsList, path=self.outputpath, show=show)

        return figure

# ### pool from multiprocessing

    def mppool(self, nproc):
        """
        It creates a multiprocessing objet to parallelise nested and emcee samplers.

        Parameters
        ------------
         nproc : int
            number of processors to use.

        Returns
        ---------
        pool : multiprocessing.Pool
            object

        nproc : int
            Number of processes

        """
        import multiprocessing as mp
        from multiprocessing.pool import ThreadPool
        if nproc <= 0:
            ncores = mp.cpu_count()
            nprocess = ncores//2
            logger.info("Using  {} processors of {}.".format(nprocess, ncores))
        elif nproc == 1:
            logger.info("Using 1 processor")
            nprocess = None
            pool = None
        else:
            nprocess = nproc
            ncores = mp.cpu_count()
            logger.info("Using {} processors of {} .".format(nprocess, ncores))

        if nprocess != None:
                pool = mp.Pool(processes=nprocess)

        return pool, nprocess

    def neuralLike(self, iniFile=None, **kwargs):
        """
        Under construction.
        This method trains a neural network in order to learn the likelihood function.
        """
        from simplemc.analyzers.neuralike.NeuralManager import NeuralManager
        self.outputpath = '{}_neuralike'.format(self.outputpath)
        if iniFile:
            ndivsgrid = self.config.getint('neuralike', 'ndivsgrid', fallback=50)
            epochs = self.config.getint('neuralike', 'epochs', fallback=500)
            learning_rate = self.config.getfloat('neuralike', 'learning_rate', fallback=5e-4)
            batch_size = self.config.getint('neuralike', 'batch_size', fallback=32)
            psplit = self.config.getfloat('neuralike', 'psplit', fallback=0.8)
            hidden_layers_neurons = [int(x) for x in self.config.get('neuralike', 'hidden_layers_neurons',
                                                                     fallback=[100, 100, 200]).split(',')]
            nproc = self.config.getint('neuralike', 'nproc', fallback=3)
        else:
            ndivsgrid = kwargs.pop('ndivsgrid', 50)
            epochs = kwargs.pop('epochs', 500)
            learning_rate = kwargs.pop('learning_rate', 5e-4)
            batch_size = kwargs.pop('batch_size', 32)
            psplit = kwargs.pop('psplit', 0.8)
            hidden_layers_neurons = kwargs.pop('hidden_layers_neurons', [100, 100, 200])
            nproc = kwargs.pop('nproc', 3)
        if nproc > 1:
            import multiprocessing as mp
            pool = mp.Pool(processes=nproc)
        else:
            pool = None

        return NeuralManager(self.logLike, self.bounds, self.root, ndivsgrid=ndivsgrid,
                             epochs=epochs, hidden_layers_neurons=hidden_layers_neurons, psplit=psplit,
                             learning_rate=learning_rate, batch_size=batch_size, pool=pool)