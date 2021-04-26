"""Keras neural net predictor.

This implements a Keras Sequential model (a deep MLP)

Author: Martin White (martin.white@adelaide.edu.au)
Date: December 2018

Modified for SimpleMC use by I Gomez-Vargas (igomezv0701@alumno.ipn.mx)
Date: June 2020
"""
import numpy
import sys

try:
    from sklearn.preprocessing import StandardScaler
except:
    import warnings
    warnings.warn("Please install sklearn library if you want to use neural networks in pybambi")

try:
   import tensorflow as tf
except:
    import warnings
    warnings.warn("Please install tensorflow library if you want to use neural networks")

class KerasNetInterpolation:
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
    def __init__(self, params, logL, split=0.8, numNeurons=300, epochs=100, model=None,
                 savedmodelpath=None):

        params = numpy.array(params)
        logL = numpy.array(logL)
        batch_size = 8

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
        ntrain = int(split * nparams)
        indx = [ntrain]
        self.params_training, self.params_testing = numpy.split(params, indx)
        self.logL_training, self.logL_testing = numpy.split(logL, indx)

        # create scaler
        scaler = StandardScaler()
        # fit scaler on data
        scaler.fit(self.params_training)
        # apply transform
        self.params_training = scaler.transform(self.params_training)
        self.params_testing = scaler.transform(self.params_testing)

        if savedmodelpath is not None:
            self.model = tf.keras.models.load_model(savedmodelpath)
        elif model is None:
            self.numNeurons = numNeurons
            self.model = self._default_architecture()
        else:
            self.model = model

        callbacks = [tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min',
                                   min_delta=0.0,
                                   patience=5,
                                   restore_best_weights=True),
                     tf.keras.callbacks.ReduceLROnPlateau(patience=2)]

        self.history = self.model.fit(self.params_training,
                                      self.logL_training,
                                      validation_data=(self.params_testing,
                                                       self.logL_testing),
                                      epochs=epochs, batch_size=batch_size,
                                      callbacks=callbacks)

    def _default_architecture(self):

        n_cols = self.params_training.shape[1]

        model = tf.keras.models.Sequential([
            tf.keras.layers.Dense(self.numNeurons, activation='relu', input_shape=(n_cols,)),
            tf.keras.layers.Dense(self.numNeurons, activation='relu'),
            tf.keras.layers.Dense(self.numNeurons, activation='relu'),
            tf.keras.layers.Dense(1)
        ])
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

        if loglikelihood > self._maxLogL + 0.1 \
                           or loglikelihood < self._minLogL - 0.1:
            inRange = False
        return inRange
