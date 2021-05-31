


class BaseLikelihood:
    """
    This is the base likelihood object. Base for all other likelihoods.
    """

    def __init__(self, name):
        self.name_ = name

    def name(self):
        """
        Name of the dataset

        :return: name
        """
        return self.name_

    def setTheory(self, theory):
        """
        Define a theoretical model to use in the likelihood

        Parameters
        ------------
        theory : object
            Instance of BaseCosmology class.
            For example, theory = LCDMCosmology()


        :return: theory
        """
        self.theory_ = theory

    def theory(self):
        """
        :return: theory
        """
        return self.theory_

    def freeParameters(self):
        """
        :return: free parameters of the theory
        """
        return self.theory_.freeParameters()

    def updateParams(self, params):
        """
        Update values of the model parameters

        Parameters
        ----------
        params : list
            List of instance of Parameter class.

       :return: list of updated parameters
        """
        return self.theory_.updateParams(params)

    def loglike(self):
        return 0.0

    def theory_loglike_prior(self):
        return self.theory_.prior_loglike()

    def loglike_wprior(self):
        return self.loglike()+self.theory_.prior_loglike()
