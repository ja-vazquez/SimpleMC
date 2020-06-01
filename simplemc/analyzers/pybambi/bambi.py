"""Driving routine for pyBAMBI.

Author: Will Handley (wh260@cam.ac.uk)
Date: November 2018

Modified for SimpleMC use by I Gomez-Vargas (2020)
"""
import os
from .pybambimanager import BambiManager


def loglike_thumper(loglikelihood, nDims, **kwargs):
    """loglike_thumper.

    Parameters
    ----------
    sampler: str
        Choice of sampler. Options: `['mnest', 'snest']`.
        Default `'multinest'`.

    nlive: int
        Number of live points.
        Default `nDims*25`

    root: str
        root of filename.
        Default `'chains/<nested_sampler>'`

    num_repeats: int
        number of repeats for polychord.
        Default `nDims*5`

    eff: float
        efficiency for multinest.
        Default `0.5**nDims`

    learner: object
        information indicating what learning algorithm to use for approximating
        the likelihood. Can be the string `'keras'`, or a `keras.models.Model`
        Default `'keras'`

    ntrain: int
        Number of training points to use
        Default `nlive`

    proxy_tolerance: float
        Required accuracy of proxy.
        Default `0.01`


    """
    # Process kwargs
    sampler = kwargs.pop('sampler', 'mnest')
    nlive = kwargs.pop('nlive', nDims*25)
    root = kwargs.pop('root', os.path.join('chains', sampler, '/'))
    num_repeats = kwargs.pop('num_repeats', nDims*5)
    eff = kwargs.pop('eff', 0.5**nDims)
    learner = kwargs.pop('learner', 'keras')
    proxy_tolerance = kwargs.pop('proxy_tolerance', 0.1)
    failure_tolerance = kwargs.pop('failure_tolerance', 0.5)
    ntrain = kwargs.pop('ntrain', nlive)
    if kwargs:
        raise TypeError('Unexpected **kwargs: %r' % kwargs)

    # Set up the global manager of the BAMBI session.
    thumper = BambiManager(loglikelihood, learner, proxy_tolerance,
                           failure_tolerance, ntrain)

    return thumper.loglikelihood()

    


