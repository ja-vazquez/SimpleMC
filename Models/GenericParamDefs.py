"""
This module is for arbitrary parameters of an arbitrary model. 
"""

from Parameter import *

## Parameters are name, value, variation, bounds, LaTeX name, "\Omega_{b}h^2")
a_par = Parameter("a", 0., 0.5, (-10., 10.), "a" )
b_par = Parameter("b", 0., 0.5, (-10., 10.), "b")