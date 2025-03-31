
from simplemc.likelihoods.BaseLikelihood import BaseLikelihood


class BBNLikelihood (BaseLikelihood):
    def __init__(self, mean=0.02218, sigma=0.00055):
        """
        A.-K. Burns, T.M.P. Tait and M. Valli, PRyMordial: the first
        three minutes, within and beyond the standard model,
        European Physical Journal C 84 (2024) 86 [2307.07061].
        .
        Parameters
        ----------
        mean
        sigma

        Returns
        -------

        """
        BaseLikelihood.__init__(self, "BBN_%g_pm_%g" % (mean, sigma))
        self.mean = mean
        self.pr = 1/(2*sigma**2)

    def loglike(self):
        return -self.pr*(self.mean-self.theory_.Obh2)**2


class BBN(BBNLikelihood):
    def __init__(self):
        BBNLikelihood.__init__(self, 0.02218, 0.00055)
