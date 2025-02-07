

from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
import scipy as sp
import numpy as np

class CompositeLikelihood(BaseLikelihood):
    def __init__(self, llist=[]):
        """
        This thing takes likelihoods and cooks one
        which is a composite of many.

        Parameters
        ----------
        llist

        Returns
        -------

        """
        BaseLikelihood.__init__(self, "Composite")
        self.llist_ = llist


    def setTheory(self, theory):
        BaseLikelihood.setTheory(self, theory)
        for l in self.llist_:
            l.setTheory(theory)
            assert(self.theory_ is l.theory_)


    def addLikelihood(self, like):
        self.llist_.append(like)
        if hasattr(self, 'theory_'):
            like.setTheory(self.theory_)
            assert(self.theory_ is like.theory_)


    def addLikelihoods(self, likes):
        for like in likes:
            self.addLikelihood(like)


    def compositeNames(self):
        return [like.name() for like in self.llist_]


    def compositeLogLikes(self):
        return np.array([like.loglike() for like in self.llist_])


    def compositeLogLikes_wprior(self):
        return np.array([like.loglike() for like in self.llist_] + [self.theory_.prior_loglike()])


    def loglike(self):
        likes = [like.loglike() for like in self.llist_]
        return np.array(likes).sum()

    # The base is good enough as
    # Everybody is holding the same reference
    # def updateParams(self,params):
        # ok=True
        # for l in self.llist_:
        #    ok=ok and l.updateParams(params)
        # return ok
