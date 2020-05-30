"""Base predictor class.

Author: Will Handley (wh260@cam.ac.uk)
Date: November 2018
"""
import numpy


class Predictor(object):
    """Base predictor class.

    This takes in a training set params -> logL, and aims to construct a
    mapping between them.

    Parameters
    ----------
    params:
        `numpy.array` of physical parameters to train on
        shape (ntrain, ndims)

    logL:
        `numpy.array` of loglikelihoods to learn
        shape (ntrain,)

    """
#IGV: ntrain is 80% (split) of multinest sampling points as in arXiv:1110.2997
    def __init__(self, params, logL, split=0.8):
        """Construct predictor from training data."""
        params = numpy.array(params)
        logL = numpy.array(logL)

        if len(params) != len(logL):
            raise ValueError("input and target must be the same length")
        elif params.ndim != 2:
            raise ValueError("input must be two-dimensional")
        elif logL.ndim != 1:
            raise ValueError("target must be one-dimensional")

        nparams = len(params)
        randomize = numpy.random.permutation(nparams)
        params = params[randomize]
        logL = logL[randomize]

        self._maxLogL = numpy.max(logL)
        self._minLogL = numpy.min(logL)
        ntrain = int(split*nparams)
        indx = [ntrain]
        self.params_training, self.params_testing = numpy.split(params, indx)
        self.logL_training, self.logL_testing = numpy.split(logL, indx)

    def __call__(self, x):
        """Calculate proxy loglikelihood.

        Parameters
        ----------
        x:
            `numpy.array` of physical parameters to predict

        Returns
        -------
        proxy loglikelihood value(s)

        """
        err = "Predictor: You need to implement a call function"
        raise NotImplementedError(err)

    def uncertainty(self):
        """Uncertainty value for the trained model."""
        err = "Predictor: You need to implement an uncertainty function"
        raise NotImplementedError(err)

    def valid(self, loglikelihood):
        """Check validity of proxy.

        Checks to see if the supplied log likelihood value is within the
        current range of likelihoods, including the uncertainty

        Parameters
        ----------
        loglikelihood:
            Value of the log likelihood that needs checking

        """
        inRange = True
        if loglikelihood > self._maxLogL + self.uncertainty() \
                or loglikelihood < self._minLogL - self.uncertainty():
            inRange = False
        return inRange
