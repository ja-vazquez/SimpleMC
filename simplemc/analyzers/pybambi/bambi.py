"""Driving routine for pyBAMBI.

Author: Will Handley (wh260@cam.ac.uk)
Date: November 2018

Modified for SimpleMC use by I Gomez-Vargas (igomezv0701@alumno.ipn.mx)
Date: June 2020
"""
import os
from .pybambimanager import BambiManager
from simplemc.analyzers.dynesty import dynesty
import multiprocessing as mp

nprocess = 2

pool = mp.Pool(processes=nprocess)

def loglike_thumper(loglikelihood, prior, nDims, **kwargs):
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
    nlive = kwargs.pop('nlive', nDims*25)
    learner = kwargs.pop('learner', 'keras')
    #proxy_tolerance = kwargs.pop('proxy_tolerance', 0.01)
    proxy_tolerance = kwargs.pop('proxy_tolerance', 10.0)
    failure_tolerance = kwargs.pop('failure_tolerance', 0.5)
    # ntrain = kwargs.pop('ntrain', nlive)
    ntrain = kwargs.pop('ntrain', 100)
    split = kwargs.pop('split', 0.8)
    numNeurons = kwargs.pop('numNeurons', 200)
    epochs = kwargs.pop('epochs', 0.8)
    model = kwargs.pop('model', None)
    savedmodelpath = kwargs.pop('savedmodelpath', None)
    simpleLike = kwargs.pop('simpleLike', None)

    if kwargs:
        raise TypeError('Unexpected **kwargs: %r' % kwargs)
    # Set up the global manager of the BAMBI session.
    thumper = BambiManager(loglikelihood, learner, proxy_tolerance,
                           failure_tolerance, ntrain, split=split,
                           numNeurons=numNeurons, epochs=epochs, model=model,
                           savedmodelpath=savedmodelpath)

    s = dynesty.NestedSampler(thumper.loglikelihood, prior, nDims,
                  bound='multi', sample='unif', nlive=nlive)

    s.run_nested(dlogz=0.01, simpleLike=simpleLike, dumper=thumper.dumper)



    


