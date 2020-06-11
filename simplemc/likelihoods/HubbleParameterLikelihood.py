


from simplemc.likelihoods.BaseLikelihood import BaseLikelihood


class HubbleParameterLikelihood (BaseLikelihood):
    def __init__(self, mean=0.70, sigma=0.02):
        """
         Distance ladder hubble likelihoods.
        Parameters
        ----------
        mean
        sigma

        Returns
        -------

        """
        BaseLikelihood.__init__(self, "Hubble_%g_pm_%g" % (mean, sigma))
        self.mean = mean
        self.pr = 1/(2*sigma**2)

    def loglike(self):
        return -self.pr*(self.mean-self.theory_.h)**2


class RiessH0(HubbleParameterLikelihood):
    def __init__(self):
        HubbleParameterLikelihood.__init__(self, 0.738, 0.024)
