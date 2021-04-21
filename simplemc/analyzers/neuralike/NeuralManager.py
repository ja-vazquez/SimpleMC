from .GridLikes import GridLikes
from .NeuralNet import NeuralNet
import os


class NeuralManager:
    def __init__(self, loglikelihood, pars_bounds, rootname, ndivsgrid=5, hidden_layers_neurons=None,
                 epochs=100, plot=True, **kwargs):
        if hidden_layers_neurons is None:
            hidden_layers_neurons = [100, 100, 100]
        self.loglikelihood = loglikelihood
        self.pars_bounds = pars_bounds
        self.ndivsgrid = ndivsgrid
        self.epochs = epochs
        self.plot = plot
        self.grid_path = 'simplemc/analyzers/neuralike/neural_models/{}'.format(rootname)
        self.model_path = 'simplemc/analyzers/neuralike/neural_models/{}.h5'.format(rootname)
        self.fig_path = 'simplemc/analyzers/neuralike/neural_models/{}.png'.format(rootname)

        self.learning_rate = kwargs.pop('learning_rate', 5e-4)
        self.batch_size = kwargs.pop('batch_size', 32)
        self.early_tol = kwargs.pop('early_tol', 100)
        self.psplit = kwargs.pop('psplit', 0.8)
        self.topology = [len(pars_bounds)] + hidden_layers_neurons + [1]

        if not self.modelChecker():
            self.training()
        self.neural_model = self.load()


    def training(self):
        grid = GridLikes(self.loglikelihood, self.pars_bounds, ndivs=self.ndivsgrid, files_path=self.grid_path)
        samples, likes = grid.make_dataset()
        neural_model = NeuralNet(X=samples, Y=likes, topology=self.topology, epochs=self.epochs,
                                 batch_size=self.batch_size, learrning_rate=self.learning_rate)

        neural_model.train()
        neural_model.save_model('{}'.format(self.model_path))
        if self.plot:
            neural_model.plot(save=True, figname='{}'.format(self.fig_path))

        return True

    def load(self):
        neural_model = NeuralNet(load=True, model_path=self.model_path)
        return neural_model

    def modelChecker(self):
        """
        This method checks if the name of the model.h5 exists, if it already exists use it, otherwise train a
        new neural net in order to generate a new model.
        """
        if os.path.isfile('{}'.format(self.model_path)):
            return True
        else:
            return False

    def loglikelihood(self, params):
        return self.neural_model.predict(params)



