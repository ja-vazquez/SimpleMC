"""Keras neural net predictor.

This implements a Keras Sequential model (a deep MLP)

Author: Martin White (martin.white@adelaide.edu.au)
Date: December 2018

"""
import numpy
from base import Predictor
#from neuralnetworks.base import Predictor
from keras.models import Sequential
from keras.layers import Dense
from keras.callbacks import EarlyStopping
#IGV: Se agrega la siguiente
from keras.layers.normalization import BatchNormalization

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
    #IGV: ntrain is 80% (split) of multinest sampling points as in arXiv:1110.2997
    def __init__(self, params, logL, split=0.8, model=None):
        """Construct predictor from training data."""
        super(KerasNetInterpolation, self).__init__(params, logL, split)

        if model is None:
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
                                      epochs=300,
                                      callbacks=callbacks)

####################### This is the original _default_architecture. Uncomment for recovery

    def _default_architecture(self):
        # Create model
        model = Sequential()

        # Number of neurons in each hidden layer, could make this configurable?
        numNeurons = 200

        # Get number of input parameters
        # Note: if params contains extra quantities (ndim+others),
        # we need to change this
        n_cols = self.params_training.shape[1]

        # Add model layers, note choice of activation function (relu)
        # We will use 3 hidden layers and an output layer
        # Note: in a Dense layer, all nodes in the previous later connect
        # to the nodes in the current layer

        model.add(Dense(numNeurons, activation='relu', input_shape=(n_cols,)))
        model.add(Dense(numNeurons, activation='relu'))
        model.add(Dense(numNeurons, activation='relu'))
        model.add(Dense(1))

        # Now compile the model
        # Need to choose training optimiser, and the loss function
        model.compile(optimizer='adam', loss='mean_squared_error')

        return model
 
####################### This is the end of the original _default_architecture

############IGV: begin modifications of _default_architecture
#    def _default_architecture(self):
        # Create model
#        model = Sequential()

        # Number of neurons in each hidden layer, could make this configurable?
#        numNeurons = 200

        # Get number of input parameters
        # Note: if params contains extra quantities (ndim+others),
        # we need to change this
#        n_cols = self.params_training.shape[1]

        # Add model layers, note choice of activation function (relu)
        # We will use 3 hidden layers and an output layer
        # Note: in a Dense layer, all nodes in the previous later connect
        # to the nodes in the current layer

#        model.add(Dense(numNeurons, input_shape=(n_cols,), init=uniform))
        #model.add(BatchNormalization())
        #model.add(Activation('relu'))
        #model.add(Dropout(0.5))

        #IGV: we can think of this chunk as the hidden layer    
        #model.add(Dense(numNeurons, init='uniform'))
        #model.add(BatchNormalization())
        #model.add(Activation('relu'))
        #model.add(Dropout(0.5))

        #model.add(Dense(numNeurons, init='uniform'))     
        #model.add(BatchNormalization())
        #model.add(Activation('relu'))
        #model.add(Dropout(0.5))

        #model.add(Dense(1))


        #IGV: The next 4 lines are the original layers
#        model.add(Dense(numNeurons, activation='relu', input_shape=(n_cols,)))
#        model.add(Dense(numNeurons, activation='relu'))
#        model.add(Dense(numNeurons, activation='relu'))
#        model.add(Dense(1))

        # Now compile the model
        # Need to choose training optimiser, and the loss function
        #IGV:
        #sgd = SGD(lr=0.1, decay=1e-6, momentum=0.9, nesterov=True)
        #model.compile(loss='mean_squared_error', optimizer=sgd)

        # Original:
#        model.compile(optimizer='adam', loss='mean_squared_error')

#        return model

############IGV: end modifications of _default_architecture



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
