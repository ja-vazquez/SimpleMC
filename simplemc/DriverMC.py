

from simplemc.analyzers.MaxLikeAnalyzer import MaxLikeAnalyzer
from simplemc.analyzers.SimpleGenetic import SimpleGenetic
from simplemc.analyzers.MCMCAnalyzer import MCMCAnalyzer
from simplemc.analyzers.dynesty import dynesty
from simplemc.cosmo.Derivedparam import AllDerived
from simplemc.runbase import ParseModel, ParseDataset
from simplemc.PostProcessing import PostProcessing
from scipy.special import ndtri
import multiprocessing as mp
from simplemc import logger
import numpy as np
import sys, os
import time


#TODO import corner
#TODO keras neural network
#TODO GR criteria in .ini file
#TODO Remove Nestle?
#TODO join genetic likelihood with the other one
class DriverMC:
    """
        This class is the manager and wrapper between all
        the analyzers and the pertinent functions.
    """
    def __init__(self, iniFile=None, **kwargs):
        """
        Read the input parameters or ini file.

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
            or a optimizer: {maxlike, genetic}

        addDerived : bool
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
                            'analyzername {"nested", "mcmc", "maxlike", "emcee" , "genetic"} Default: mcmc'
                            '\n\tchainsdir Default: SimpleMC_chains\n\t')
                sys.exit(1)


        #Initialize the Theory & Datasets
        T = ParseModel(self.model, custom_parameters=self.custom_parameters,
                                   custom_function=self.custom_function)
        L = ParseDataset(self.datasets, path_to_data=self.path_to_data,
                                        path_to_cov=self.path_to_cov, fn=self.fn)

        if self.prefact == "pre":  T.setVaryPrefactor()
        if self.varys8  == "True": T.setVarys8()
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



    def executer(self, **kwargs):
        """
        This is a wrapper of the runners of the analyzer in order to make
        easier the execution, mainly if is through an ini file.

        Parameters
        ----------
        All the parameters of the methods:

            - mcmcRunner
            - nestedRunner
            - emceeRunner
            - geneticRunner
            - maxlikeRunner

        """
        if self.analyzername == 'mcmc':
            self.mcmcRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'nested':
            self.nestedRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'emcee':
            self.emceeRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'maxlike':
            self.maxLikeRunner(iniFile=self.iniFile, **kwargs)
        elif self.analyzername == 'genetic':
            self.geneticRunner(iniFile=self.iniFile, **kwargs)
        else:
            sys.exit("{}: Sampler/Analyzer name invalid".format(self.analyzername))
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

        self.custom_parameters = self.config.get(   'custom', 'custom_parameters', fallback=None)
        self.custom_function   = self.config.get(   'custom', 'custom_function',   fallback=None)
        ## Following two are for custom data
        self.path_to_data = self.config.get(        'custom', 'path_to_data', fallback=None)
        self.path_to_cov  = self.config.get(        'custom', 'path_to_cov',  fallback=None)
        self.fn = self.config.get('custom', 'path_to_cov',  fallback="generic")
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

        chainno : int
            Number of chains in parallel.

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
            chainno  = self.config.getint(      'mcmc', 'chainno', fallback=1)
            GRstop   = self.config.getfloat(    'mcmc', 'GRstop',  fallback=0.01)
            evidence = self.config.getboolean(  'mcmc', 'evidence',fallback=False)
        else:
            nsamp    = kwargs.pop('nsamp', 50000)
            skip     = kwargs.pop('skip',  300)
            temp     = kwargs.pop('temp',  2)
            chainno  = kwargs.pop('chainno', 1)
            GRstop   = kwargs.pop('GRstop', 0.01)
            evidence = kwargs.pop('evidence', False)
            if kwargs:
                logger.critical('Unexpected **kwargs for MCMC: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use default values.\n'
                            'MCMC executer kwargs are:\n\tnsamp (int) Default: 50000\n\t'
                            'skip (int) Default 300\n\ttemp (float) Default: 2.0'
                            '\n\tchainno (int) Default: 1\n\t'
                            'evidence (bool) Default: False')
                sys.exit(1)
                #raise TypeError('Unexpected **kwargs: {}'.format(kwargs))
        logger.info("\n\tnsamp: {}\n\tskip: {}\n\t"
                    "temp: {}\n\tchain num: {}\n\tevidence: {}".format(
                    nsamp, skip, temp, chainno, evidence))
        if self.analyzername is None: self.analyzername = 'mcmc'
        self.outputpath = "{}_{}".format(self.outputpath, self.analyzername)
        #Check whether the file already exists
        self.outputChecker()
        ti = time.time()

        #Main process
        M = MCMCAnalyzer(self.L, self.outputpath, skip=skip, nsamp=nsamp, temp = temp,
                        chain_num=chainno, addDerived=self.addDerived, GRstop=GRstop)

        self.ttime = time.time() - ti

        #Compute Bayesian Evidence
        if evidence:
            try:
                from MCEvidence import MCEvidence
                logger.info("Aproximating bayesian evidence with MCEvidence (arXiv:1704.03472)\n")
                MLE = MCEvidence(self.outputpath + ".txt" ).evidence()
                self.result = ['mcmc', M, "Evidence with MCEvidence : {}\n".format(MLE), strresult]
            except:
                #writeSummary(self.chainsdir, outputname, ttime)
                # print("Warning!")
                # print("MCEvidence could not calculate the Bayesian evidence [very small weights]\n")
                logger.error("MCEvidence could not calculate the Bayesian evidence [very small weights]")
        else:
            self.result = ['mcmc', M, "Maxlike: {}".format(M.maxloglike)]

        return True



##---------------------- Nested ----------------------

    def nestedRunner(self, iniFile=None, **kwargs):
        """
        This method calls Dynesty samplers.

        Parameters
        ___________
        engine : str
            Use dynesty or nestle library

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
            self.engine = self.config.get(          'nested', 'engine',       fallback='dynesty')
            dynamic     = self.config.getboolean(   'nested', 'dynamic',      fallback=False)
            neuralNetwork = self.config.getboolean( 'nested', 'neuralNetwork',fallback=False)
            nestedType  = self.config.get(          'nested', 'nestedType',   fallback='multi')
            nlivepoints = self.config.getint(       'nested', 'nlivepoints',  fallback=1024)
            accuracy    = self.config.getfloat(     'nested', 'accuracy',     fallback=0.01)
            nproc       = self.config.getint(       'nested', 'nproc',        fallback=1)

            self.priortype = self.config.get('nested', 'priortype', fallback='u')
            #nsigma is the default value for sigma in gaussian priors
            self.nsigma = self.config.get('nested', 'sigma', fallback=2)

            if neuralNetwork:
                split      = self.config.getfloat('neural', 'split', fallback=0.8)
                numNeurons = self.config.getint(  'neural', 'numNeurons', fallback=100)
        else:
            self.engine = kwargs.pop('engine',    'dynesty')
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


            if kwargs:
                logger.critical('Unexpected **kwargs for nested sampler: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'Nested executer options are:\n\tnlivepoints (int) Default: 1024\n\t'
                            'accuracy (float) Default: 0.01\n\tpriortype ({"u", "g"}) Default: "u"\n\t'
                            'nestedType {"multi", "single", "balls", "cubes"} Default: "multi"\n\t'
                            'neuralNetwork (bool) Default: True\n\t'
                            'dynamic (bool) Default: False\n\t'
                            'addDerived (bool) Default: True\n\t'
                            'engine {"dynesty", "nestle"} Default: "dynesty"')
                sys.exit(1)

        #stored output files
        if self.analyzername is None: self.analyzername = 'nested'
        self.outputpath = '{}_{}_{}_{}'.format(self.outputpath, self.analyzername, self.engine, nestedType)
        if neuralNetwork: self.outputpath = "{}_ANN".format(self.outputpath)
        self.outputChecker()

        #paralel run
        pool, nprocess = self.mppool(nproc)
        logger.info("\n\tnlivepoints: {}\n"
                    "\taccuracy: {}\n"
                    "\tnested type: {}\n"
                    "\tengine: {}".format(nlivepoints, accuracy, nestedType, self.engine))

        ti = time.time()
        if neuralNetwork:
            logger.info("\tUsing neural network.")
            #from bambi import run_pyBAMBI
            from simplemc.analyzers.pybambi.bambi import loglike_thumper
            learner = 'keras'
            # kerasmodel = None
            self.logLike = loglike_thumper(self.logLike, learner=learner, \
                                             ntrain=nlivepoints, nDims=self.dims)

        if dynamic:
            logger.info("\nUsing dynamic nested sampling...")
            sampler = dynesty.DynamicNestedSampler(self.logLike, self.priorTransform,
                      self.dims, bound=nestedType, pool=pool, queue_size=nprocess)

            sampler.run_nested(nlive_init=nlivepoints, dlogz_init=0.05, nlive_batch=100,
                            maxiter_init=10000, maxiter_batch=1000, maxbatch=10,
                            outputname=self.outputpath, addDerived=self.addDerived, simpleLike=self.L)
            M = sampler.results


        elif self.engine == 'dynesty':
            sampler = dynesty.NestedSampler(self.logLike, self.priorTransform, self.dims,
                        bound=nestedType, sample = 'unif', nlive = nlivepoints,
                        pool = pool, queue_size=nprocess)
            sampler.run_nested(dlogz=accuracy, outputname=self.outputpath,
                               addDerived=self.addDerived, simpleLike=self.L)
            M = sampler.results


        elif self.engine == 'nestle':
            try:
                import nestle
                M = nestle.sample(self.logLike, self.priorTransform, ndim=self.dims, method=nestedType,
                npoints=nlivepoints, dlogz=accuracy, callback=nestle.print_progress,
                pool = pool, queue_size = nprocess)
            except ImportError as error:
                sys.exit("{}: Please install nestle module"
                         "or use dynesty engine".format(error.__class__.__name__))
        else:
            sys.exit('wrong selection')
        try:
            pool.close()
        except:
            pass
        self.ttime = time.time() - ti
        self.result = ['nested', M, M.summary(), 'nested :{}'.format(nestedType),
                       'dynamic : {}'.format(dynamic), 'ANN :{}'.format(neuralNetwork),
                       'engine: {}'.format(self.engine)]
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
            walkers = self.config.getint('emcee', 'walkers', fallback=500)
            nsamp   = self.config.getint('emcee', 'nsamp', fallback=20000)
            burnin  = self.config.getint('emcee', 'burnin', fallback=0)
            nproc   = self.config.getint('emcee', 'nproc', fallback=1)
        else:
            walkers = kwargs.pop('walkers', 30)
            nsamp   = kwargs.pop('nsamp', 20000)
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
        try:
            import emcee
            ti = time.time()
            sampler = emcee.EnsembleSampler(walkers, self.dims, self.logPosterior, pool=pool)
            # pass the initial samples and total number of samples required
            sampler.run_mcmc(inisamples, nsamp + burnin, progress=True)
            # extract the samples (removing the burn-in)
            # postsamples = sampler.chain[:, burnin:, :].reshape((-1, self.dims))
            self.ttime = time.time() - ti
        except ImportError as error:
            sys.exit("{}: Please install this module"
                     "or try using other sampler".format(error.__class__.__name__))
        # set up the sampler
        self.burnin = burnin
        try:
            pool.close()
        except:
            pass
        # fig = corner.corner(postsamples)
        # fig.savefig('emcee.png')
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

        plot_par1 : bool

        plot_par2 : bool

        """
        if self.analyzername is None: self.analyzername = 'maxlike'
        self.outputpath = '{}_{}_optimization'.format(self.outputpath, self.analyzername)
        self.outputChecker()
        if iniFile:
            withErrors = self.config.getboolean('maxlike', 'withErrors', fallback=False)
            showplot   = self.config.getboolean('maxlike', 'showplot', fallback=False)
            plot_par1  = self.config.get('maxlike', 'plot_par1', fallback=None)
            plot_par2  = self.config.get('maxlike', 'plot_par2', fallback=None)
            showderived= self.config.getboolean('maxlike', 'showderived', fallback=False)

        else:
            withErrors  = kwargs.pop('withErrors', False)
            showplot    = kwargs.pop('showplot', False)
            plot_par1   = kwargs.pop('plot_par1', None)
            plot_par2   = kwargs.pop('plot_par2', None)
            showderived = kwargs.pop('showderived', False)
            if kwargs:
                logger.critical('Unexpected **kwargs for MaxLike: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'MaxLikeAnalyzer executer options are:'
                            '\n\twithErrors (bool) Default: False')
                sys.exit(1)
        ti = time.time()
        A = MaxLikeAnalyzer(self.L, self.model, withErrors=withErrors, showderived=showderived,\
                            showplot=showplot, param1=plot_par1, param2=plot_par2)
        params = self.T.printParameters(A.params)
        self.ttime = time.time() - ti
        self.result = ['maxlike', A, params]
        return True


##---------------------- Genetic Algorithms ----------------------


    def geneticRunner(self, iniFile=None, **kwargs):
        """
        It calls SimpleGenetic class.

        Parameters
        -----------
        n_individuals : int
            Number of individuals.

        n_generations : int
            Number of generations.

        selection_method : str
            Selection method {tournament, roulette, rank}

        mut_prob : float
            Probability of mutation.
        """
        if self.analyzername is None: self.analyzername = 'genetic'
        self.outputpath = '{}_{}_optimization'.format(self.outputpath, self.analyzername)
        self.outputChecker()

        if iniFile:
            n_individuals = self.config.getint('genetic', 'n_individuals', fallback=400)
            n_generations = self.config.getint('genetic', 'n_generations' , fallback=1000)
            selection_method = self.config.get('genetic', 'selection_method', fallback='tournament')
            mut_prob = self.config.getfloat('genetic', 'mut_prob', fallback=0.6)
        else:
            n_individuals = kwargs.pop('n_individuals', 400)
            n_generations = kwargs.pop('n_generations', 1000)
            selection_method = kwargs.pop('selection_method', 'tournament')
            mut_prob = kwargs.pop('mut_prob', 0.6)
            if kwargs:
                logger.critical('Unexpected **kwargs for genetic optimizer: {}'.format(kwargs))
                logger.info('You can skip writing any option and SimpleMC will use the default value.\n'
                            'genetic executer options are:'
                            '\n\tn_individuals (int) Default: 400\n\t'
                            'n_generations (int) Default: 1000\n\t'
                            'selection_method {"tournament","rank","roulette"} Default: "tournament"'
                            '\n\tmut_prob (float) Default: 0.4')
                sys.exit(1)
        logger.info("\n\tn_individuals: {}\n\tn_generations: {}"
                    "\n\tselection method: {}\n\t"
                    "mut prob: {}".format(n_individuals, n_generations,
                                        selection_method,mut_prob))
        ti = time.time()

        M = SimpleGenetic(self.logLikeGenetic, self.dims, self.bounds,
                          n_individuals=n_individuals,
                          n_generations=n_generations,
                          prob_mut=mut_prob, method_selection=selection_method,
                          outputname=self.outputpath)

        self.ttime = time.time() - ti
        self.result = ['genetic', M]
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
        """Prior Transform for gaussian and flat priors"""
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
        for i, element in enumerate(v):
            values.append(element)
        assert len(self.pars_info) == len(values)
        for pars, val in zip(self.pars_info, values):
            pars.setValue(val)

        self.T.updateParams(self.pars_info)
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
        self.paramFiles()

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
        parfile = self.outputpath + ".paramnames"
        fpar = open(parfile, 'w')
        for p in cpars:
            fpar.write(p.name + "\t\t\t" + p.Ltxname + "\n")
        if self.addDerived:
            AD = AllDerived()
            for pd in AD.list:
                fpar.write(pd.name + "\t\t\t" + pd.Ltxname + "\n")
        if self.analyzername == 'mcmc' or (self.analyzername == 'nested' and self.engine=='dynesty'):
            if (self.L.name() == "Composite"):
                self.sublikenames = self.L.compositeNames()
                for name in self.sublikenames:
                    fpar.write(name + "_like \t\t\t" + name + "\n")
                fpar.write("theory_prior \t\t\t None \n")


    def postprocess(self, summary=True, stats=False, addtxt=None):
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
        if self.analyzername == 'nested':
            pp = PostProcessing(self.result, self.paramsList, self.outputpath,
                engine=self.engine, addDerived=self.addDerived, loglike=self.L)
            if self.engine == 'nestle':
                pp.saveNestedChain()
        elif self.analyzername == 'emcee':
            pp = PostProcessing(self.result, self.paramsList, self.outputpath,
                                skip=self.burnin, addDerived=self.addDerived, loglike=self.L)
            pp.saveEmceeSamples()
        else:
            pp = PostProcessing(self.result, self.paramsList, self.outputpath)
        if summary:
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
        figure = SimplePlotter(self.chainsdir, self.outputpath, self.paramsList, show=show)

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


