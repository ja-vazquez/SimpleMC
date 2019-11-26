#!/usr/bin/env python
import MCMCAnalyzer
from BaseLikelihood import BaseLikelihood
from Parameter import Parameter
import sys
import os
# print sys.path
sys.path = ["py", "../py"]+sys.path


class testLike(BaseLikelihood):

    def __init__(self):
        BaseLikelihood.__init__(self, "testlike")
        self.x = 0
        self.y = 0

    def freeParameters(self):
        return [Parameter("x", 0.1, err=1, bounds=(0, 2),), Parameter("y", 0, err=1)]

    def updateParams(self, params):
        for p in params:
            if p.name == "x":
                self.x = p.value
            if p.name == "y":
                self.y = p.value

    def loglike_wprior(self):
        x = self.x
        y = self.y
        return (-(x**2+y**2-2*x*y)/2.)


L = testLike()
MCMCAnalyzer.MCMCAnalyzer(L, "output_test")
