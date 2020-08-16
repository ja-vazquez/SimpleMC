import numpy as np

class ToyModels:
    def __init__(self, model):
        """
        This class contains some toy models to test nested samplers

        Parameters
        ----------
        model : str
            {'egg', 'ring', 'gaussian', 'square', 'himmel'}
        """
        # self.bounds contains x and y bounds.
        if model == 'egg':
            self.bounds = [[0., 1.], [0., 1.]]
            self.bounds_z = [-100., 300.]
            self.loglike = self.eggLoglike

        elif model in ['ring', 'gaussian', 'himmel', 'square']:
            self.bounds = [[-5., 5.], [-5., 5.]]
            if model == 'ring':
                self.bounds_z = [0., 10.]
            else:
                self.bounds_z = [0., 1.]

    def eggLoglike(self, cube):
        tmax = 5.0 * np.pi
        t = 2.0 * tmax * cube - tmax
        return (2.0 + np.cos(t[0] / 2.0) * np.cos(t[1] / 2.0)) ** 5.0

    def himmelLoglike(self, cube):
        return -(cube[0] ** 2 + cube[1] - parameterlist[0]) ** 2.0 - (cube[0] + cube[1] ** 2 - parameterlist[1]) ** 2

    def gaussLoglike(self, x):
        return -((x[0]) ** 2 + (x[1]) ** 2 / 2.0 - 1.0 * x[0] * x[1]) / 2.0

    def ringLoglike(self, x):
        r2 = x[0] ** 2 + x[1] ** 2
        return -(r2 - 4.0) ** 2 / (2 * 0.5 ** 2)

    def squareLoglike(self, x):
        sq = 1.
        if abs(x[0]) > 5 or abs(x[1]) > 5:
            sq = 0.
        return sq

    def priorTransform(self, theta):
        priors = []
        for c, bound in enumerate(self.bounds):
            priors.append(theta[c] * (bound[1] - bound[0]) + bound[0])
        return np.array(priors)










