#!/usr/bin/env python
from RunBase import *
from scipy.optimize import minimize as minf
A = DR11LyaAuto()
B = DR11LyaCross()


def fu(x):

    print(x, A.loglike_aperp_apar(x[0], x[1]), end=' ')
    print(B.loglike_aperp_apar(x[0], x[1]), end=' ')
    toret = -A.loglike_aperp_apar(x[0], x[1])
    toret += -B.loglike_aperp_apar(x[0], x[1])
    print(toret)
    return toret


print(minf(fu, [1., 1.], method='SLSQP', bounds=[(0.85, 1.15), (0.85, 1.15)]))
