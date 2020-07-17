"""Torch neural net predictor.

This implements a Torch Sequential model

Author: I Gomez-Vargas (igomezv0701@alumno.ipn.mx)
Date: July 2020
"""
import numpy as np
import sys

try:
    import torch
    from torch.autograd import Variable
    from torch import optim
except:
    print("You need to install torch if want to use torch neural net. \n"
          "Keras can be used instead.")

class TorchNetInterpolation:
    """Torch neural net interpolation.

    Returns the loglikelihood from a Torch neural net-based interpolator

    Trains a basic 3-layer neural network with 200 neurons per layer.

    Parameters
    ----------
    params:
        `np.array of` physical parameters to train on
        shape (ntrain, ndims)

    logL:
        `np.array` of loglikelihoods to learn
        shape (ntrain,)

    """
    # IGV: ntrain is 80% (split) of multinest sampling points as in arXiv:1110.2997
    def __init__(self, params, logL, split=0.8, numNeurons=200, epochs=300, model=None,
                 savedmodelpath=None):

        params = np.array(params)
        logL = np.array(logL)

        if len(params) != len(logL):
            raise ValueError("input and target must be the same length")
        elif params.ndim != 2:
            raise ValueError("input must be two-dimensional")
        elif logL.ndim != 1:
            raise ValueError("target must be one-dimensional")

        nparams = len(params)
        randomize = np.random.permutation(nparams)
        params = params[randomize]
        logL = logL[randomize]

        self._maxLogL = np.max(logL)
        self._minLogL = np.min(logL)
        self._lastLogL = logL[-1]
        ntrain = int(split * nparams)
        indx = [ntrain]
        self.params_training, self.params_testing = np.split(params, indx)
        self.logL_training, self.logL_testing = np.split(logL, indx)

        self.numNeurons = numNeurons
        self.model = self._default_architecture()

        # define a loss function
        self.loss = torch.nn.MSELoss(reduction='mean')
        # define an optimizer
        optimizer = optim.Adam(self.model.parameters())

        Xtrain = torch.from_numpy(self.params_training).float()
        Ytrain = torch.from_numpy(self.logL_training).float()
        Xtest = torch.from_numpy(self.logL_testing).float()

        epochs = 15
        batch_size = 32
        n_batches = nparams // batch_size

        costs = []
        test_accuracies = []
        for i in range(epochs):
            cost = 0.
            for j in range(n_batches):
                Xbatch = Xtrain[j * batch_size:(j + 1) * batch_size]
                Ybatch = Ytrain[j * batch_size:(j + 1) * batch_size]
                cost += self.train(self.model, self.loss, optimizer, Xbatch, Ybatch)

            # Ypred = self.predict(self.model, Xtest)
            # acc = np.mean(self.logL_testing == Ypred)
            # print("Epoch: %d, cost: %f, acc: %.2f" % (i, cost / n_batches, acc))
            #
            # # for plotting
            # costs.append(cost / n_batches)
            # test_accuracies.append(acc)

    def _default_architecture(self):

        n_cols = self.params_training.shape[1]

        # the model will be a sequence of layers
        model = torch.nn.Sequential()

        model.add_module("dense1", torch.nn.Linear(n_cols, self.numNeurons))
        model.add_module("relu1", torch.nn.ReLU())
        model.add_module("dense2", torch.nn.Linear(self.numNeurons, self.numNeurons))
        model.add_module("relu2", torch.nn.ReLU())
        model.add_module("dense3", torch.nn.Linear(self.numNeurons, 1))

        return model

    def train(self, model, loss, optimizer, inputs, labels):
        inputs = Variable(inputs, requires_grad=False)
        labels = Variable(labels, requires_grad=False)
        # Reset gradient
        optimizer.zero_grad()
        # Forward
        logits = model.forward(inputs)
        self.output = loss.forward(logits, labels)
        # Backward
        self.output.backward()
        # Update parameters
        optimizer.step()

        return self.output.item()

    def predict(self, model, inputs):
        inputs = Variable(inputs, requires_grad=False)
        logits = model.forward(inputs)
        return logits.data.numpy().argmax(axis=1)

    def __call__(self, x):
        """Calculate proxy loglikelihood.

        Parameters
        ----------
        x:
            `np.array` of physical parameters to predict

        Returns
        -------
        proxy loglikelihood value(s)

        """
        x_ = np.atleast_2d(x)
        x_ = torch.from_numpy(x_).float()
        x_.retain_grad()
        y = self.predict(self.model, x_)
        return float(np.squeeze(y))

    def uncertainty(self):
        """Uncertainty value for the trained keras model."""
        test_loss = np.sqrt(np.array(self.output.item()))
        #np.sqrt(self.history.history['val_loss'])

        return np.squeeze(test_loss.min())

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
