"""Keras neural net predictor.

This implements a Keras Sequential model (a deep MLP)

Author: Martin White (martin.white@adelaide.edu.au)
Date: December 2018

"""
import numpy
from .base import Predictor
from keras.models import Sequential
from keras.layers import Dense
from keras.callbacks import EarlyStopping


class KerasNetInterpolation(Predictor):
    """Keras neural net interpolation.

    Returns the loglikelihood from a Keras neural net-based interpolator

    Trains a basic 3-layer neural network with 200 neurons per layer.

    Parameters
    ----------
    params:
        `numpy.array of` physical parameters to train on
        shape (ntrain, ndims)

    logL:
        `numpy.array` of loglikelihoods to learn
        shape (ntrain,)

    """
    # IGV: ntrain is 80% (split) of multinest sampling points as in arXiv:1110.2997
    def __init__(self, params, logL, split=0.8, numNeurons=200, epochs=300, model=None):
        """Construct predictor from training data."""
        super(KerasNetInterpolation, self).__init__(params, logL, split)
        if model is None:
            self.numNeurons = numNeurons
            self.model = self._default_architecture()
        else:
            self.model = model

        callbacks = [EarlyStopping(monitor='val_loss', mode='min',
                                   min_delta=0.001,
                                   patience=10,
                                   restore_best_weights=True)]

        self.history = self.model.fit(self.params_training,
                                      self.logL_training,
                                      validation_data=(self.params_testing,
                                                       self.logL_testing),
                                      epochs=epochs,
                                      callbacks=callbacks)

    def _default_architecture(self):
        # Create model
        model = Sequential()

        # Get number of input parameters
        # Note: if params contains extra quantities (ndim+others),
        # we need to change this
        n_cols = self.params_training.shape[1]

        # Add model layers, note choice of activation function (relu)
        # We will use 3 hidden layers and an output layer
        # Note: in a Dense layer, all nodes in the previous later connect
        # to the nodes in the current layer

        model.add(Dense(self.numNeurons, activation='relu', input_shape=(n_cols,)))
        model.add(Dense(self.numNeurons, activation='relu'))
        model.add(Dense(self.numNeurons, activation='relu'))
        model.add(Dense(1))

        # Now compile the model
        # Need to choose training optimiser, and the loss function
        model.compile(optimizer='adam', loss='mean_squared_error')

        return model

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
        x_ = numpy.atleast_2d(x)
        y = self.model.predict(x_)
        return float(numpy.squeeze(y))

    def uncertainty(self):
        """Uncertainty value for the trained keras model."""
        test_loss = numpy.sqrt(self.history.history['val_loss'])

        return numpy.squeeze(test_loss.min())
