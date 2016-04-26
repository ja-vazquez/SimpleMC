#
# This thing takes likelihoods and cooks one
# which is a composite of many.
#
from BaseLikelihood import *
from scipy import *

class CompositeLikelihood (BaseLikelihood):
    def __init__(self, llist=[]):
        BaseLikelihood.__init__(self,"Composite")
        self.llist_ = llist

    def setTheory(self,theory):
        BaseLikelihood.setTheory(self,theory)
        for l in self.llist_:
            l.setTheory(theory)
            assert(self.theory_ is l.theory_)

    def addLikelihood (self,like):
        self.llist_.append(like)
        if hasattr(self, 'theory_'):
            like.setTheory(self.theory_)
            assert(self.theory_ is like.theory_)

    def addLikelihoods (self,likes):
        for like in likes:
            self.addLikelihood(like)

    def compositeNames(self):
        return [like.name() for like in self.llist_]

    def compositeLogLikes(self):
        return array([like.loglike() for like in self.llist_])

    def compositeLogLikes_wprior(self):
        return array([like.loglike() for like in self.llist_]+[self.theory_.prior_loglike()])

    def loglike(self):
        likes=[like.loglike() for like in self.llist_]
        return array(likes).sum()
    
    ## The base si good enough as 
    # Everybody is holding the same reference
    #def updateParams(self,params):
        #ok=True
        #for l in self.llist_:
        #    ok=ok and l.updateParams(params)
        #return ok
            
            



