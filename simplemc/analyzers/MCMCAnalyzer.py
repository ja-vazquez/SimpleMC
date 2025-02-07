

#TODO print the Accepted # for each rank


from simplemc.cosmo.Derivedparam import AllDerived
import scipy.linalg as la
import scipy as sp
import copy
import random
import sys
import numpy as np
from numpy import loadtxt


class MCMCAnalyzer:
    """
    MCMC sampler (Metropolis-Hastings).
    This is the MCMC module. It spits out chains that are compatible with CosmoM
    it calculates cov matrix during burn-in.
    optional temperature makes it sample at a higher temperature but note that
    this guy, as opposed to cosmomc, reweights the weights on the fly.

    Parameters
    -----------
    like : Likelihood object
        Object of a Likelihood class.

    outfile : str
         Output file.

    skip : int
        Burn-in.

    nsamp : int
        Number of mcmc samples.

    temp : float
        Temperature

    cov : numpy.array
        Covariance matrix. Default: None.

    addDerived : bool
        In order to ad derived parameters such as age of the universe and Omega_{Lambda}.

    GRstop : float
        Gelman-Rubin criteria.
    """
    def __init__(self, like, outfile, skip=5000, nsamp=100000, temp=1.0,
                 cov=None, chain_num=None, addDerived=False, GRstop=0.01, checkGR=500):
        try:
            from mpi4py import MPI
            self.comm = MPI.COMM_WORLD
            name = MPI.Get_processor_name()
            self.chain_num = self.comm.rank+1  #chain_num
            print("Hello, World! "
                  "I am process %d of %d on %s" %
                  (self.comm.rank, self.comm.size, name))
        except:
            self.chain_num = 1
            print("Running only 1 chain without MPI.")

        self.like      = like
        self.outfile   = outfile
        self.nsamp     = nsamp
        self.skip      = skip
        self.temp      = float(temp)  # temperature
        self.cpars     = like.freeParameters()
        self.N         = len(self.cpars)
        self.derived   = addDerived
        self.checkgr   = checkGR

        minvals, maxvals = [], []
        for lb, hb in [p.bounds for p in self.cpars]:
            minvals.append(lb)
            maxvals.append(hb)
        self.minvals = np.array(minvals)
        self.maxvals = np.array(maxvals)
        print("Bounds:", self.minvals, self.maxvals)

        if (like.name() == "Composite"):
            self.sublikenames = like.compositeNames()
            self.composite = True
        else:
            self.composite = False

        if (cov == None):
            # make initial cov matrix from diagonal "errors"
            errs = [0.01*p.error**2 for p in self.cpars]
            self.init_pcov(sp.diag(errs))
        else:
            self.init_pcov(cov)

        if self.derived: self.AD = AllDerived()
        self.GRcondition = GRstop

        self.RunChain()


    def init_pcov(self, mat):
        self.chol = la.cholesky(mat)


    def RunChain(self):
        self.openFiles()
        self.cloglike, self.cloglikes = self.getLikes()
        # set up logofs based on the first log like which should be
        # the same for all chains. Better than nothing.
        #self.logofs=self.cloglike
        # Actually, above doesn't seem to work very well.
        # Instead, use zero, as our likelihoods never became very large
        self.logofs = 0
        # current weight
        self.cw     = 0
        # current counter
        self.co     = 0
        # mean for burin
        self.swx    = 0
        self.meanx  = sp.zeros(self.N)
        self.meanxx = sp.zeros((self.N, self.N))
        # max loglike
        self.maxloglike = -1e30
        # are we done
        self.done   = False
        #converge
        self.lpars   = []
        #Uses the last percentage of the chain
        self.percen  = 0.4

        self.gr = None

        print("Starting chain...")

        while not (self.done):
            ppars, numout = self.GetProposal()
            self.cw += numout  ## things hitting outside the prior are formally rejected samples
            self.like.updateParams(ppars)
            ploglike, ploglikes = self.getLikes()
            if (sp.isnan(ploglike)):
                print("\nSomething bad has happened, nan in loglike, assuming zero log")
                ploglike = -1e50
            # print cloglike, ploglike, [p.value for p in like.freeParameters()], [p.value for p in self.cpars]
            if (ploglike > self.cloglike):
                accept = True
            else:
                accept = (np.exp((ploglike-self.cloglike)/self.temp)
                          > random.uniform(0., 1.))

            # print [p.value for p in ppars], accept, ploglike
            # stop
            if (accept):
                self.ProcessAccepted(ppars, ploglike, ploglikes)
                # for i, item in enumerate(ppars):
                #     print("\nPPARS", i, item.value, ploglike, ploglikes)
            else:
                self.cw += 1

            print("Accepted: {:d} | loglike: {:3.4f} | "
                  "Gelman-Rubin: {}".format(self.co, self.cloglike, self.gr), end='\r')
            sys.stdout.flush()
            if (self.co >0 and self.co % self.checkgr == 0):
                try:
                    chains = self.comm.gather(self.lpars, root=0)
                    if self.comm.rank ==0:
                        self.gr = self.GRDRun(chains)
                        #print('Gelman-Rubin R-1:', self.gr)
                        if (sp.all(self.gr < self.GRcondition)):
                            condition = 1
                            self.closeFiles()
                        else:
                            condition = 0
                    else:
                            condition = None
                    recvmsg = self.comm.bcast(condition, root=0)
                    if recvmsg ==1:
                        print('\n---- Gelman-Rubin achived ---- ')
                        self.closeFiles()
                        return True
                except:
                # Without mpi4py installed
                    self.gr = self.GRDRun(self.lpars)
                    if (sp.all(self.gr < self.GRcondition)):
                        print('\n---- Gelman-Rubin achived ---- ')
                        self.closeFiles()
                        return True

    def GRDRun(self, chains):
        """
        This is a implementation of the Gelman Rubin diagnostic.
        If the number of chains is 1, then this method divides it in two
        and does the diagnostic for convergence.

        Parameters
        ----------
        chains : list
            List with the chains to perform the GR-diagnostic.

        Returns
        -------
        result : float
            Gelman-Rubin diagnostic.

        """
        mean_chain = []
        var_chain  = []

        if len(chains) == 1:
            lchain = len(chains[0])//2
            chains = [chains[0][:lchain], chains[0][lchain:]]
        else:
            clen = [len(chain) for chain in chains]
            if len(set(clen)) == 1:
                lchain = clen[0]
            else:
                #print('take same # steps', clen)
                lchain = min(clen)

        try:
            for chain in chains:
                mean_chain.append(sp.mean(chain[-lchain:], axis=0))
                var_chain.append(sp.var(chain[-lchain:], axis=0))
        except:
            return 1

        M = sp.mean(mean_chain, axis=0)
        W = sp.mean(var_chain,  axis=0)

        B= sum([(b-M)**2 for b in mean_chain])
        B = lchain/(len(chains)- 1.)*B
        R = (1. - 1./lchain)*W +  B/lchain

        result = np.array(sp.absolute(1- np.sqrt(R/W)))
        return result



    def openFiles(self):
        """
        Open the files to save the samples and maxlike.
        Also add the Dervided Parameters if addDerived option is True.

        """
        outfile = self.outfile
        formstr = '%g ' + '%g '*(self.N+1)
        if self.derived:
            formstr += '%g '*(len(self.AD.list))

        if (self.composite):
            formstr += '%g '*(len(self.sublikenames)+1)
        formstr += '\n'

        if (self.chain_num == None):
            self.cfname  = outfile + ".txt"
            mlfname = outfile + ".maxlike"
        else:
            self.cfname  = outfile + "_%i.txt" % (self.chain_num)
            mlfname = outfile + "_%i.maxlike" % (self.chain_num)

        self.fout    = open(self.cfname, 'w')
        self.mlfout  = open(mlfname, 'w')
        self.formstr = formstr


    def closeFiles(self):
        chain = loadtxt(self.cfname)
        self.weights = chain[:, 0]
        self.samples = chain[:, 2:self.N+2]
        self.loglikes = chain[:, 1]
        self.mlfout.close()



    def getLikes(self):
        """
        Get loglikelihoods values from the used data.
        """
        if (self.composite):
            cloglikes = self.like.compositeLogLikes_wprior()
            cloglike  = cloglikes.sum()
        else:
            cloglikes = []
            cloglike  = self.like.loglike_wprior()
        return cloglike, cloglikes


    def GetProposal(self):
        """
        Generation of proposal point in mcmc.

        """
        vec = sp.zeros(self.N)
        numreject = 0
        while True:
            ppars = copy.deepcopy(self.cpars)
            step  = self.draw_pcov()
            #print ('step #', [p.value for p in  ppars])
            for i, p in enumerate(ppars):
                p.value += step[i]
                vec[i]   = p.value

            if all(vec > self.minvals) and all(vec < self.maxvals):
                return ppars, numreject
            numreject += 1



    def draw_pcov(self):
        a = np.array([random.gauss(0., 1,) for _ in range(self.N)])
        return np.dot(a, self.chol)



    def ProcessAccepted(self, ppars, ploglike, ploglikes):
        self.co += 1
        # if (self.co % 100 == 0): #JAV 1000
        #     print("Accepted samples", self.co, self.cw)
        vec = [p.value for p in self.cpars]

        #for convergence
        self.lpars.append(vec)
        if self.co % 10 ==0:
            del self.lpars[:int((1-self.percen)*10)]


        if (self.co > self.skip):
            # weight rescaled
            wers = self.cw*np.exp((self.cloglike-self.logofs)
                               * (self.temp-1.0)/self.temp)

            tmp = [wers, -self.cloglike] + vec
            if self.derived:
                tmp += [pd.value for pd in self.AD.listDerived(self.like)]

            if (self.composite):
                outstr = self.formstr % tuple(tmp + self.cloglikes.tolist())
            else:
                outstr = self.formstr % tuple(tmp)
            self.fout.write(outstr)
            # Flush file on regular basis
            if (self.co % 100 == 0): #JAV 1000
                self.fout.flush()

            if (self.cloglike > self.maxloglike):
                self.maxloglike = self.cloglike
                #print("New maxloglike", self.maxloglike)
                self.mlfout.seek(0)
                self.mlfout.write(outstr)
                self.mlfout.flush()
                self.maxloglike

            if self.co > self.nsamp:
                print('Number of steps achieved')
                self.done = True
                self.closeFiles()

        elif (self.co < self.skip):
            self.swx += self.cw
            v = np.array(vec)
            self.meanx  += v*self.cw
            self.meanxx += sp.outer(v, v)*self.cw
            if (self.cw > 30):
                print("\nStill burning in, weight too large")
                self.chol *= 0.9
                # print(self.cw)
        else:  # co==skip
            self.meanx  /= self.swx
            self.meanxx /= self.swx
            self.meanxx -= sp.outer(self.meanx, self.meanx)
            print("\nRe-initializing covariance matrix after burn-in")
            print(self.meanxx)
            print()
            # for i, p in enumerate(self.cpars):
            #     print("{}: {} +/- {}".format(p.name, p.value, np.sqrt(self.meanxx[i, i])))

            self.init_pcov(self.meanxx)

        self.cw    = 1
        self.cpars = ppars
        self.cloglike = ploglike
        if self.composite:
            self.cloglikes = ploglikes

    def get_results(self):
        return {'samples': self.samples, 'weights': self.weights, 'loglikes': self.loglikes,
                'gr_diagnostic': self.gr}
