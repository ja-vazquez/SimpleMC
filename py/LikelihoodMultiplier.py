# This is a helper class that takes a likelihood
# and just multiplies loglike with some value to simulater
# better or shittier data.
##
# You have to give python a credit for how effortlessly
# one can hack classes!!
##


def LikelihoodMultiplier(LikeInstance, mult):
    class MultipliedLikelihood(LikeInstance.__class__):
        def __init__(self, LikeInstance, mult):
            # copy  instance into my class
            self.__dict__.update(LikeInstance.__dict__)
            # save type since we need to call parent type
            self.LikeType = LikeInstance.__class__
            self.mult = mult

        def loglike(self):
            return self.LikeType.loglike(self)*self.mult

    return MultipliedLikelihood(LikeInstance, mult)
