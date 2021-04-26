"""BAMBI management object.

Author: Pat Scott (p.scott@imperial.ac.uk)
Date: Feb 2019

Modified for SimpleMC and dynesty use by I Gomez-Vargas (igomezv0701@alumno.ipn.mx)
Date: June 2020
"""

from simplemc.analyzers.pybambi.kerasnet import KerasNetInterpolation
import sys

try:
    import tensorflow as tf
except:
    import warnings
    warnings.warn("Please install tensorflow library if you want to use neural networks")


class BambiManager(object):
    """Does all the talking for BAMBI.

    Takes a new set of training data from the dumper and trains (or retrains) a
    neural net, and assesses whether or not it can be used for a given
    parameter combination.

    """

    def __init__(self, loglikelihood, proxy_tolerance,
                 failure_tolerance, updInt, split=0.8, numNeurons=200, epochs=300,
                 model=None, savedmodelpath=None, it_to_start_net=None, dlogz_start=None):
        """Construct bambi object."""

        self.proxy_tolerance = proxy_tolerance
        self._loglikelihood = loglikelihood
        self._proxy_tolerance = proxy_tolerance
        self._failure_tolerance = failure_tolerance
        self._proxy_trained = False
        self.old_learners = []
        self.split = split
        self.numNeurons = numNeurons
        self.epochs = epochs
        self.model = model
        self.savedmodelpath = savedmodelpath
        self.it_to_start_net = it_to_start_net
        self.updInt = updInt
        self.dlogz_start = dlogz_start

    def make_learner(self, params, loglikes):
        """Construct a Predictor."""
        return KerasNetInterpolation(params, loglikes,
                                     split=self.split, numNeurons=self.numNeurons,
                                     epochs=self.epochs, model=self.model,
                                     savedmodelpath=self.savedmodelpath)

    def dumper(self, params, live_loglks=None, dlogz=1e4, it=0):
        """It sends datasets of physical points and likelihoods to neural net"""
        # params is a dictionary if dynesty is running in parallel
        # or the numpy array of the live_params if is running without parallelism.
        if live_loglks is None:
            # if live_loglks is None -> dynesty is running in parallel.
            live_loglks = params['live_logl']
            live_params = params['live_v']
            counter = params['it'] - 3
            dlogz = params['dlogz']
        else:
            # dynesty isn't running in parallel
            live_params = params
            counter = it - 3

        if self.dlogz_start:
            if dlogz <= self.dlogz_start and counter < self.it_to_start_net:
                self.it_to_start_net = counter

        mod = counter % self.updInt
        if mod == 0 and self._proxy_trained is False and counter >= self.it_to_start_net:
            self.train_new_learner(live_params, live_loglks)

        return self._proxy_trained

    def loglikelihood(self, params):
        """Bambi Proxy wrapper for original loglikelihood."""
        # Short circuit to the full likelihood if proxy not yet fully trained
        if not self._proxy_trained:
            return self._loglikelihood(params)

        # Call the learner
        candidate_loglikelihood = self._current_learner(params)

        # If the learner can be trusted, use its estimate,
        # otherwise use the original like and update the failure status
        if self._current_learner.valid(candidate_loglikelihood):
            return candidate_loglikelihood
        else:
            self._rolling_failure_fraction = (1.0 + (self.updInt - 1.0) *
                                              self._rolling_failure_fraction
                                              ) / self.updInt
            if self._rolling_failure_fraction > self._failure_tolerance:
                self._proxy_trained = False
            return self._loglikelihood(params)

    def train_new_learner(self, params, loglikes):
        """Train a new Predictor."""
        try:
            self.old_learners.append(self._current_learner)
        except AttributeError:
            pass

        self._current_learner = self.make_learner(params, loglikes)
        sigma = self._current_learner.uncertainty()
        print("\nCurrent uncertainty in network log-likelihood predictions: {}".format(sigma))

        if sigma < self._proxy_tolerance:
            self._proxy_trained = True
            self._rolling_failure_fraction = 0.0