#
# This is the base likelihood object. Base for all other likelihoods.
#


class BaseLikelihood:

    def __init__(self, name):
        self.name_ = name

    def name(self):
        return self.name_

    def setTheory(self, theory):
        self.theory_ = theory

    def theory(self):
        return self.theory_

    def freeParameters(self):
        return self.theory_.freeParameters()

    def updateParams(self, params):
        return self.theory_.updateParams(params)

    def loglike(self):
        return 0.0

    def theory_loglike_prior():
        return self.theory_.prior_loglike()

    def loglike_wprior(self):
        return self.loglike()+self.theory_.prior_loglike()
