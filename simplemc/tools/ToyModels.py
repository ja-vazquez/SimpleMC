#!/usr/bin/env python
import numpy as np
import time
from scipy.special import ndtri

class ToyModels:
    def __init__(self, model):
        """
        model_bounds are lists [infx, supx, infy, supy, infz, supz]
        """
        self.ring_bounds  = [-5., 5., -5., 5., 0., 10.]
        self.gauss_bonds = [-5., 5., -5., 5., 0., 1.]
        self.egg_bounds = [0., 1., 0., 1., -100., 300.]
        self.himmel_bounds = [-5., 5., -5., 5., 0., 1.]
        self.square_bounds = [-5., 5., -5., 5., 0., 1.]

    def eggLoglike(self, cube):
        tmax = 5.0 * np.pi
        t = 2.0 * tmax * cube - tmax
        return (2.0 + np.cos(t[0] / 2.0) * np.cos(t[1] / 2.0)) ** 5.0

    def himmelLoglike(self, cube):
        return -(cube[0] ** 2 + cube[1] - 11) ** 2.0 - (cube[0] + cube[1] ** 2 - 7) ** 2

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

    def genericPriorTransform(self, cube):
        return cube * (suppr - infpr) + infpr











