##
## Chain Iterator,
##
## Iterate over public chains and get a theory at each point in the chain
## 

from cosmich import *
from RunBase import *

class ChainIterator:
    def __init__(self, directory, model, pre, data):
        self.chain=cosmochain("%s/%s_%s_%s"%(directory,model,pre,data),'auto')
        #print self.chain.parcol['mnu']
        #stop
        self.T=ParseModel(model)
        self.L=ParseDataset(data)
        if pre=="pre":
            T.setVaryPrefactor()
        self.pars=self.T.freeParameters()
        self.L.setTheory(self.T)
        self.pindices=[self.chain.parcol[name] for name in [par.name for par in self.pars]]
        self.N=self.chain.N


    def update_theory(self,i):
        for j in range(len(self.pars)):
            self.pars[j].setValue(self.chain.chain[i,self.pindices[j]])
            #print self.pars[j].name, self.chain.chain[i,self.pindices[j]]
        self.T.updateParams(self.pars)

    def theory(self,i):
        self.update_theory(i)
        return self.T

    def pvalue (self, i, name):
        return self.chain.chain[i,self.chain.parcol[name]]

    def weight(self,i):
        return self.chain.chain[i,0]

    def chi2(self,i):
        return self.chain.chain[i,1]*2
