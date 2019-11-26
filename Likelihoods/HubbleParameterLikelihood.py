#
# Distance ladder hubble likelihoods.
#


from BaseLikelihood import *
from scipy import *
import scipy.linalg as la


class HubbleParameterLikelihood (BaseLikelihood):
    def __init__(self, mean=0.70, sigma=0.02):
        BaseLikelihood.__init__(self, "Hubble_%g_pm_%g" % (mean, sigma))
        self.mean = mean
        self.pr = 1/(2*sigma**2)

    def loglike(self):
        return -self.pr*(self.mean-self.theory_.h)**2


class RiessH0(HubbleParameterLikelihood):
    def __init__(self):
        HubbleParameterLikelihood.__init__(self, 0.738, 0.024)
