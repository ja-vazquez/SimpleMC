import sys,os
sys.path=["py","../py"]+sys.path
import pandas as pd
import numpy as np 
import os

class GelmanRubinDiagnostic:
    """Gelman Rubin Diagnostics class.
    
    counter is the current number of stpes in the MCMC method. """

    def __init__ (self, df, freeParams):
        self.df = df
        self.counter = len(df)
        #self.df = np.array(df).reshape(counter, 2*freeParams+2)
        data = np.array(df)
        print(self.counter)
        self.freeParams = freeParams
        print(np.shape(data))

        burnin = int(self.counter*0.1)
        midle = int((self.counter-burnin)/2)

        burninarray = np.array([x for x in range(burnin,self.counter)])
        print(np.shape(burninarray))
       
        #even = np.array([x for x in range(burnin,self.counter) if x%2==0])
        #odd = np.array([x for x in range(burnin,self.counter) if x%2==1])
        firsthalf = np.array([x for x in range(burnin, self.counter//2)])
        secondhalf = np.array([x for x in range(self.counter//2, self.counter)])

        columns = np.array([x for x in range(1,self.freeParams+2)])
       
        #data2 = data[burninarray[:,None],columns]
        
        self.chain1 = np.ndarray.tolist(data[firsthalf[:,None],columns])
        self.chain2 = np.ndarray.tolist(data[secondhalf[:,None],columns])
        
        #one = np.array([x for x in range(burnin,self.counter) if x%5==4])
        #two = np.array([x for x in range(burnin,self.counter) if x%5==3])
        #three = np.array([x for x in range(burnin,self.counter) if x%5==2])
        #four = np.array([x for x in range(burnin,self.counter) if x%5==1])
        #five = np.array([x for x in range(burnin,self.counter) if x%5==0])
        #self.chain1 = np.ndarray.tolist(data[one[:,None],columns])
        #self.chain2 = np.ndarray.tolist(data[two[:,None],columns])
        #self.chain3 = np.ndarray.tolist(data[three[:,None],columns])
        #self.chain4 = np.ndarray.tolist(data[four[:,None],columns])
        #self.chain5 = np.ndarray.tolist(data[five[:,None],columns])




    def GRDRun(self):
        """This is a implementation of the Gelman Rubin diagnostic"""
        chains = [self.chain1, self.chain2]
        chains = np.array(chains)
        mean_chain = []
        var_chain  = []
        #n_params puede ser sust por self.N
        n_params = self.freeParams
        n_steps = int(self.counter/2)
        chain_num = 2
        #bi is the burn-in, right?
        bi = int(n_steps*0.1)
        for chain in chains:
            mean_chain.append(np.mean(chain, axis=0))
            var_chain.append(np.var(chain, axis=0))
                
        M = np.mean(mean_chain, axis=0)
        W = np.mean(var_chain,  axis=0)
        B = 0.

        for i in np.arange(chain_num):
            B += (mean_chain[i] - M)**2

        B = n_steps/(chain_num - 1.)*B
        R = (1. - 1./n_steps)*W +  B/n_steps
      
        result = np.array(np.abs(np.sqrt(R/W) -1))
        print ('Gelman-Rubin Diagnostic:', result)

        return result


