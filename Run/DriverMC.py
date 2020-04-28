import sys
sys.path = ["Analizers", "Cosmo", "pybambi", "py"] + sys.path

from MaxLikeAnalyzer import MaxLikeAnalyzer
from MCMCAnalyzer import MCMCAnalyzer
from Parameter import Parameter
from RunBase import ParseModel, ParseDataset
from scipy.special import ndtri
from SimpleGenetic import SimpleGenetic
from PostProcessing import PostProcessing 

import dynesty
import emcee

import multiprocessing as mp
import numpy as np
import time
import corner

class DriverMC():
    """
        This class is the manager and wrapper between all 
        the analyzers and the perninent functions.
    """
    def __init__(self, baseConfig):
        #nsigma is the default value for sigma in gaussian priors
        self.nsigma  = 4.
        self.inifile = baseConfig
        self.iniReader()   
        self.T, self.L = self.TLinit()
        #if you like to change the bounds, you need edit ParamDefs.py file.
        self.bounds, self.means    = self.priorValues()
        self.dims, self.paramsList = self.getDims()
        #self.n_sigma = self.dims * self.nsigma
        self.outputname = self.model + "_" + self.prefact + \
            "_" + self.datasets + "_" + self.samplername
            # \ + \ "[" + time.strftime("%H:%M:%S") + "]_"


        if self.prefact == "pre":
            self.T.setVaryPrefactor()
        self.T.printFreeParameters()

        self.result = self.executer()
        self.postprocess()

    def executer(self):
        ti = time.time()
        if self.samplername == 'mcmc':
            s = self.MCMCRunner()
        elif self.samplername == 'nested' and self.neuralNetwork == 'no':
            s = self.NestedRunner()
        elif self.samplername == 'MaxLikeAnalyzer':
            s = self.MaxLikeRunner()
        elif self.samplername == 'nested' and self.neuralNetwork == 'yes':
            s = self.BambiRunner()
        elif self.samplername == 'emcee':
            s = self.emceeRunner()
        elif self.samplername == 'genetic':
            s = self.geneticRunner()
        else:
            sys.exit("Sampler/Analyzer name invalid")
        self.ttime = time.time()- ti
        #print("\n Time elapsed : {} seconds = {} minutes".format(self.ttime, self.ttime/60.))
        return s

######################################Initialization#############################

    def iniReader(self):
        """
        It reads the ini file. 
        """

        if sys.version_info > (3, 0):
            import configparser
            config = configparser.ConfigParser()
            config.read(self.inifile)
            self.chainsdir   = config['custom']['chainsdir']
            self.model       = config['custom']['model']
            self.prefact     = config['custom']['prefact']
            self.datasets    = config['custom']['datasets']
            self.samplername = config['custom']['sampler']
            
            if self.samplername == 'mcmc':
                self.nsamp      = int(config['mcmc']['nsamp'])
                self.skip       = int(config['mcmc']['skip'])
                ## temperature at which to sample, weights get readjusted on the fly
                self.temp       = int(config['mcmc']['temp'])
                self.chainno    = int(config['mcmc']['chainno'])
                #self.GRcriteria = float(config['mcmc']['GRcriteria'])
                self.derived    = config['mcmc']['addderived']

            elif self.samplername in ['nested']:
                print('\nUsing dynesty library')
                self.nlivepoints   = int(config['nested']['nlivepoints'])
                self.accuracy      = float(config['nested']['accuracy'])
                self.priortype     = config['nested']['priortype']
                self.nestedType    = config['nested']['nestedType']
                self.neuralNetwork = config['nested']['neuralNetwork']
                self.dynamic = config['nested']['dynamic'] 
                self.nproc   = int(config['nested']['nproc'])
                self.engine  = config['nested']['engine']

                print("The number of live points is: %d"%(self.nlivepoints))
                if self.nestedType == 'multi':
                    print('\nThe sampler used is MULTINEST. Feroz et al (2009)\n')
                elif self.nestedType == 'single':
                    print('\nMukherjee, Parkinson & Liddle (2006) nested sampling\n')
                if self.neuralNetwork == 'yes': 
                    print("\nANN based on pybambi. Graff et al (2012).\n")
                if self.priortype == 'g': 
                    print("Using %d-sigmas for the gaussian prior.\n"%(self.nsigma))
                else: 
                    print("Using flat priors...")

                if self.neuralNetwork == 'yes':
                    self.split      = float(config['neural']['split'])
                    self.numNeurons = int(config['neural']['numNeurons'])

            elif self.samplername == 'emcee':

                self.ensambles = int(config['emcee']['ensambles'])
                self.samples = int(config['emcee']['samples'])
                self.burnin = int(config['emcee']['burnin'])
                self.nproc   = int(config['emcee']['nproc'])

            elif self.samplername in ['MaxLikeAnalyzer']:
                self.withErrors      = config['MaxLikeAnalyzer']['withErrors']

        else:
            print("Please use Python 3 and try again!")
            exit(1)

        
    def TLinit(self):
        """
        Returns T evaluated at the model, and L in at the datasets.
        """
        T = ParseModel(self.model)
        L = ParseDataset(self.datasets)        
        L.setTheory(T)
        return T, L


    def priorValues(self):
        """Returns some values for the priorTransform """
        bounds = self.SetPriors()[3]
        means  = np.array(self.SetPriors()[1])
        return bounds, means

    def SetPriors(self):
        """
        After setTheory, you can set the priors.

        Returns a list of lists of the object of the Parameter class.
        """
        parameters = self.T.freeParameters()
        names      = []
        values     = []
        errorlist  = []
        boundlist  = []
        latexnames = []
        for parameter in parameters:
            names.append(parameter.name)
            values.append(parameter.value)
            errorlist.append(parameter.error)
            boundlist.append(parameter.bounds)
            latexnames.append(parameter.Ltxname)
        return [names, values, errorlist, boundlist, latexnames]



## que hacer aqui?
    def getDims(self):
        """
        Returns the numbers of dimensions and a parameters list.
        """
        #IGV: We need the names of parameters in a list and, on the other hand,
        #    the dimensions. I don't found a fancy way, probably it exists.  
        freeP = self.T.freeParameters()
        listP = []
        for _, item in enumerate(freeP):
            listP.append(item.name)
        return len(listP), listP



####################################logLike and prior Transform function###########################

    def instantiatePars(self, value):
        """
        This method returns an instance of
        Parameter objects with the sampler values

        -- value : it is the generated vector by the sampler
        """
        aux=[]
        for item in value:
            aux.append(item)
        names, values, errorlist, boundlist, latexnames = self.SetPriors()
        instances = []
        #-1 is for emcee -1
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


    def priorsTransform1(self, theta):
        """prior Transform for gaussian and flat priors"""
        priors = []
        n = 2.

        if self.priortype == 'g':
            for c, bound in enumerate(self.bounds):
                mu = self.means[c]
                sigma = (bound[1]-bound[0])/n
                priors.append(mu+sigma*(ndtri(theta[c])))
        else:
            for c, bound in enumerate(self.bounds):
            #When theta 0-> append bound[0], if theta 1-> append bound[1]
                priors.append(theta[c]*(bound[1]-bound[0])+bound[0])
                #At this moment, np.array(priors) has shape (dims,)
        return np.array(priors)


    def priorTransform(self, theta):
        """
            Prior transform maps the values of the samping method theta to values between 0 and 1.
            
            Parameter:
            ----------

            theta:  it is the hypervolume definied by the priors
        """
        return self.priorsTransform1(theta)


#################################################Samplers ################################

    def MCMCRunner(self):
        """
            This method calls MCMCAnalyzer.
            Returns [MCMCAnalyzer object, time, evidence via MCEvidence (if it is possible)]

        """
        #weights = None
        #self.outputname += str(self.nsamp)
        
        M = MCMCAnalyzer(self.L, self.chainsdir + "/" + self.outputname,\
                        skip=self.skip, nsamp=self.nsamp, temp = self.temp,
                        chain_num=self.chainno, derived=self.derived)
        #,\        GRcriteria = self.GRcriteria)

        try:
            from MCEvidence import MCEvidence 
            print("Aproximating bayesian evidence with MCEvidence (arXiv:1704.03472)\n")
            MLE = MCEvidence(self.chainsdir + "/" + self.outputname + ".txt" ).evidence()
            return ['mcmc', M, "Evidence with MCEvidence : " + str(MLE) + "\n"]

        except:
            #writeSummary(self.chainsdir, outputname, ttime)
            print("Warning!") 
            print("MCEvidence could not calculate the Bayesian evidence [very small weights]\n")
            return ['mcmc', M]


    def NestedRunner(self):
        """
            This method calls Dynesty samplers:
            -- 
            -- MULTINEST
            bound : {'none', 'single', 'multi', 'balls', 'cubes'},
        """
        self.outputname += '_'+self.engine+'_'+self.nestedType
        pool, nprocess = self.mppool()

        showfiles = True
        if self.engine == 'dynesty':
            if showfiles == False:
                from dynesty import dynesty
            else:
                # Check it later
                sys.path =  ['dynesty', '/dynesty'] + sys.path
                from dynesty import dynesty


        if self.dynamic == 'yes':
            print("Using dynamic nested sampling...")

            sampler = dynesty.DynamicNestedSampler(self.logLike, self.priorTransform, \
                      self.dims, bound = self.nestedType, pool = pool, queue_size = nprocess)
            
            sampler.run_nested(nlive_init=self.nlivepoints, dlogz_init=0.05, nlive_batch=100,\
                            maxiter_init=10000, maxiter_batch=1000, maxbatch=10)

            M = sampler.results
            M.summary()

        elif self.engine == 'dynesty':
            sampler = dynesty.NestedSampler(self.logLike, self.priorTransform, self.dims,
                        bound=self.nestedType, sample = 'rwalk', nlive = self.nlivepoints,\
                        pool = pool, queue_size = nprocess)
            sampler.run_nested(dlogz=self.accuracy)
            M = sampler.results
            M.summary()

        elif self.engine == 'nestle':
            import nestle
            M = nestle.sample(self.logLike, self.priorTransform, ndim=self.dims, method=self.nestedType,
            npoints=self.nlivepoints, dlogz=self.accuracy, callback=nestle.print_progress,
            pool = pool, queue_size = nprocess)
            #M.summary()
            print('\n'+'--'*10)
            print('logz=%2.2f +/- %2.2f \n'%(M.logz, M.logzerr))
        else:
            print('wrong selection')
            sys.exit(1)
        try:
            pool.close()
        except:
            pass

        return ['nested', M, M.summary(), 'nested : ' + self.nestedType, 'dynamic : ' +\
                 self.dynamic, 'ANN : ' + self.neuralNetwork]


    def BambiRunner(self):
        from bambi import run_pyBAMBI
        self.outputname += self.nestedType + '_nested+ANN_'+str(self.nlivepoints)

        
        M = run_pyBAMBI(self.logLike, self.priorTransform, self.dims,\
            eff = self.accuracy, nestedType = self.nestedType, dynamic = self.dynamic,\
            nlive = self.nlivepoints, learner = 'keras')
        
        return ['bambi', M, 'nested : ' + self.nestedType, 'dynamic : ' +\
                 self.dynamic, 'ANN : ' + self.neuralNetwork]

    def emceeRunner(self):
        self.outputname += '_ensambles_'+str(self.ensambles)
        pool, _ = self.mppool()

        Nens = self.ensambles   # number of ensemble points

        ini = []

        for bound in self.bounds:
            ini.append(np.random.uniform(bound[0], bound[1], Nens))

        inisamples = np.array(ini).T # initial samples

        Nburnin = self.burnin   # number of burn-in samples
        Nsamples = self.samples  # number of final posterior samples

        # set up the sampler
        sampler = emcee.EnsembleSampler(Nens, self.dims, self.logPosterior, pool = pool)

        # pass the initial samples and total number of samples required
        sampler.run_mcmc(inisamples, Nsamples+Nburnin, progress=True)

        # extract the samples (removing the burn-in)
        postsamples = sampler.chain[:, Nburnin:, :].reshape((-1, self.dims))

        
        print('Number of posterior samples is {}'.format(postsamples.shape[0]))
        try:
            pool.close()
        except:
            pass 
        fig = corner.corner(postsamples)
        fig.savefig('emcee.png')
        return ['emcee', sampler, 'ensambles : ' + str(self.ensambles)]
###################################Optimizers##############################################

    def MaxLikeRunner(self):
        self.outputname += 'optimization' 
        print("Using MaxLikeAnalyzer")
        A = MaxLikeAnalyzer(self.L, withErrors = self.withErrors)
        self.T.printParameters(A.params)
        return True
       

    def geneticRunner(self):
        self.outputname += 'optimization_genetic'
        self.fgenetic = open(self.chainsdir+'/'+self.outputname + '.txt', 'w+')
        print("Using Simple Genetic Algorithm")
        M = SimpleGenetic(self.logLikeGenetic, self.dims, self.bounds)
        self.fgenetic.close()
        return ['genetic', M]


###############################Post-processing############################################
    def postprocess(self):
        if self.samplername == 'nested':
            pp = PostProcessing(self.result, self.paramsList , self.outputname,\
             chainsdir=self.chainsdir, engine=self.engine)
            pp.paramFiles(self.T, self.L)
            pp.saveNestedChain()

        elif self.samplername == 'emcee':
            pp = PostProcessing(self.result, self.paramsList, self.outputname,\
             chainsdir=self.chainsdir, skip=self.burnin)
            pp.paramFiles(self.T, self.L)
            pp.saveEmceeSamples()
        elif self.samplername == 'genetic' or self.samplername == MaxLikeAnalyzer:
            pass
        else:
            pp = PostProcessing(self.result, self.paramsList, self.outputname,\
             chainsdir=self.chainsdir)
            pp.paramFiles(self.T, self.L)
            pp.saveEmceeSamples() 
    
#Do it later -- plots and stats analis

        #if self.samplername in ['mcmc', 'nested']:
        #    stats = pp.getdistAnalyzer(cov = True)
        #    pp.writeSummary(self.ttime, stats)
        #else:
        #    pp.writeSummary(self.ttime)

        

    #def plotter(self):
    #    from PlotterMC import PlotterMC
    #    plot = PlotterMC(self.dims, chainsdir = self.chainsdir, chainsfile = self.outputname)

# ###########for genetic
    def logLikeGenetic(self, *v):
        values = []
        print("-"*10)
        print("Parameters:")
        for i, element in enumerate(v):
            values.append(element)
            print("{}".format(element),  end=' ')
        # values = np.array(values)
        strvalues = str(values).lstrip('[').rstrip(']')
        strvalues = strvalues.replace(',', '')
        listPars = self.instantiatePars(values)
        self.T.updateParams(listPars)
        self.L.setTheory(self.T)
        if (self.L.name()=="Composite"):
            cloglikes=self.L.compositeLogLikes_wprior()
            loglike=cloglikes.sum()
        else:
            loglike = self.L.loglike_wprior()
        print("\nloglike: {}".format(loglike))
        print("-" * 10)
        row = str(loglike) + " " + strvalues
        self.fgenetic.write(row + "\n")
        return loglike

############# for emcee
    def logPosterior(self, theta):
        """
        The natural logarithm of the joint posterior.
        
        Args:
            theta (tuple): a sample containing individual parameter values
            data (list): the set of data/observations
            sigma (float): the standard deviation of the data points
            x (list): the abscissa values at which the data/model is defined
        """
    
        lp = self.logPrior(theta) # get the prior
        
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
# ### pool from multiprocessing 

    def mppool(self):
        if self.nproc <= 0:
            ncores = mp.cpu_count()
            print("--"*10 )
            nprocess = ncores//2
            print("Using  {} processors of {}.".format(nprocess, ncores))
            pool = mp.Pool(processes=nprocess)
        elif self.nproc == 1:
            print("Using 1 processor.")
            print("--"*10 )
            nprocess = None
            pool = None
        else:
            nprocess = self.nproc
            print("Using {} processors.".format(nprocess))
            print("--"*10 )
            pool = mp.Pool(processes=nprocess)

        return pool, nprocess