import sys
sys.path = ["Analizers", "Cosmo", "pybambi", "py"] + sys.path

from MaxLikeAnalyzer import MaxLikeAnalyzer
from MCMCAnalyzer import MCMCAnalyzer
from Parameter import Parameter
from RunBase import ParseModel, ParseDataset
from scipy.special import ndtri
#from SimpleGenetic import SimpleGenetic

import dynesty
#import emcee
import multiprocessing as mp
import numpy as np
import time

class DriverMC():
    """
        This class is the manager and wrapper between all 
        the analyzers and the perninent functions.
    """
    def __init__(self, file):
        #nsigma is the default value for sigma in gaussian priors
        self.nsigma  = 4.
        self.inifile = file
        self.iniReader()   
        self.T, self.L = self.TLinit()
        #if you like to change the bounds, you need edit ParamDefs.py file.
        self.bounds, self.means    = self.priorValues()
        self.dims, self.paramsList = self.getDims()
        #self.n_sigma = self.dims * self.nsigma
        self.outputname = self.model + "_" + self.prefact + \
            "_" + self.datasets + "_" + self.samplername + "_" + \
            "[" + time.strftime("%H:%M:%S") + "]_"
#quitar el tiempo y talvez samplername

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
            
            elif self.samplername in ['nested']:
                print('\nUsing dynesty library')
                self.nlivepoints   = int(config['nested']['nlivepoints'])
                self.accuracy      = float(config['nested']['accuracy'])
                self.priortype     = config['nested']['priortype']
                self.nestedType    = config['nested']['nestedType']
                self.neuralNetwork = config['nested']['neuralNetwork']
                self.dynamic = config['nested']['dynamic'] 
                
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
        weights = None
        self.outputname += str(self.nsamp)
      
        
        M = MCMCAnalyzer(self.L, self.chainsdir + "/" + self.outputname,\
            skip=self.skip, nsamp=self.nsamp, temp = self.temp, chain_num=self.chainno)
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
        self.outputname += self.nestedType+'_nested_'+str(self.nlivepoints)

        ncores = mp.cpu_count()
        print("Number of processors: ", ncores)
        
        if ncores > 1: 
            nprocess = ncores//2
        else:
            nprocess = 1
        
        pool = mp.Pool(processes=nprocess)

        if self.dynamic == 'yes':
            print("Using dynamic nested sampling...")

            sampler = dynesty.DynamicNestedSampler(self.logLike, self.priorTransform, \
                self.dims, bound = self.nestedType,\
                        pool = pool, queue_size = nprocess)
            
            sampler.run_nested(nlive_init=self.nlivepoints, dlogz_init=0.05, nlive_batch=100,\
                            maxiter_init=10000, maxiter_batch=1000, maxbatch=10)

        else:
            sampler = dynesty.NestedSampler(self.logLike, self.priorTransform, self.dims,
                        bound=self.nestedType, sample = 'rwalk', nlive = self.nlivepoints,\
                        pool = pool, queue_size = ncores)
            sampler.run_nested(dlogz=self.accuracy)
        M = sampler.results   

        return ['nested', M, 'nested : ' + self.nestedType, 'dynamic : ' +\
                 self.dynamic, 'ANN : ' + self.neuralNetwork]


    def BambiRunner(self):
        from bambi import run_pyBAMBI
        self.outputname += self.nestedType + '_nested+ANN_'+str(self.nlivepoints)

        
        M = run_pyBAMBI(self.logLike, self.priorTransform, self.dims,\
            eff = self.accuracy, nestedType = self.nestedType, dynamic = self.dynamic,\
            nlive = self.nlivepoints, learner = 'keras')
        
        return ['bambi', M, 'nested : ' + self.nestedType, 'dynamic : ' +\
                 self.dynamic, 'ANN : ' + self.neuralNetwork]


###################################Optimizers##############################################

    def MaxLikeRunner(self):
        self.outputname += 'optimization' 
        print("Using MaxLikeAnalyzer")
        self.T.printFreeParameters()
        res, nloglike = MaxLikeAnalyzer(self.L).result()
        return ['MaxLikeAnalyzer', None, "\nParameters : ", res, "\nOptimal nloglike : ", nloglike ]

        #self.T.printParameters(A.params)
       


    def geneticRunner(self):
        self.outputname += 'optimization' 
        print("Usinge Simple Genetic Algorithm")
        M = SimpleGenetic(self.logLike2, self.dims, self.bounds)
        
        return ['genetic', M ]



###############################Post-processing############################################
    def postprocess(self):
        from PostProcessing import PostProcessing 
        pp = PostProcessing(self.result, self.paramsList , self.outputname,\
             chainsdir=self.chainsdir)        
        
        if self.samplername == 'nested':
            pp.paramFiles(self.T, self.L)
            pp.saveNestedChain()


#Do it later -- plots and stats analis

        #if self.samplername in ['mcmc', 'nested']:
        #    stats = pp.getdistAnalyzer(cov = True)
        #    pp.writeSummary(self.ttime, stats)
        #else:
        #    pp.writeSummary(self.ttime)

        

    def plotter(self):
        from PlotterMC import PlotterMC
        plot = PlotterMC(self.dims, chainsdir = self.chainsdir, chainsfile = self.outputname)

#############################Testing functions#################################################
    
    def logposterior(self, theta):
        lp = self.priorTransform(theta)
        if not np.isfinite(lp).any():
            return -np.inf

        return lp + self.logLike(theta)

    def emceeRunner(self):
        #pos = soln.x + 1e-4 * np.random.randn(32, 3)
        nwalkers = 100
        ndim = self.dims
        p0 = np.random.rand(nwalkers, ndim)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.logposterior)
        sampler.run_mcmc(p0, 5000) 


    def logLike2(self, *v):
        values = []
        print("Entrando al logLike2")
        for i, element in enumerate(v):
            print("Entrando al for del loglike", i)
            values.append(element)
            print("element", element)
        print("saliendo for loglike")
        #values = np.array(values)
        listPars = self.instantiatePars(values)
        self.T.updateParams(listPars)
        self.L.setTheory(self.T)
        if (self.L.name()=="Composite"):
            cloglikes=self.L.compositeLogLikes_wprior()
            loglike=cloglikes.sum()
        else:
            loglike = self.L.loglike_wprior()
        print("loglike", loglike)
        return loglike