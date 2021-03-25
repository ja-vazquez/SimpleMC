
from simplemc.likelihoods.BaseLikelihood import BaseLikelihood
from simplemc.likelihoods.JLA_SN import SN_likelihood


class SNLikelihood(BaseLikelihood):
    def __init__ (self, name, filename):
        """
        This module calculates likelihood for the full JLA SN.
        Parameters
        ----------
        ninterp

        Returns
        -------

        """

        self.name_= name
        
        # JLA with alpha, beta parameters passed in, fairly fast (one matrix inversion)
        self.like = SN_likelihood(filename,  marginalize=False)

        # JLA marginalized over alpha, beta, e.g. for use in importance sampling with no nuisance parameters.
        # Quite fast as inverses precomputed. Note normalization is not same as for alpha, beta varying.
        #like = SN_likelihood(filename,  marginalize=True)

        # as above, but very slow (but lower memory) using non-precomputed inverses (and non-threaded in python)
        #like = SN_likelihood(filename, precompute_covmats=False, marginalize=True)


    def loglike(self):
        zs = self.like.get_redshifts()
        angular_distance = [self.theory_.AD_z(z) for z in zs]
        chi2 = self.like.loglike(angular_distance, {'alpha': 0.1325237, 'beta': 2.959805}) * 2
        return -chi2/2


class JLASN_Full(SNLikelihood):
    def __init__(self):
        SNLikelihood.__init__(self, "JLASN", "simplemc/data/jla_cosmo_v2/jla.dataset")
