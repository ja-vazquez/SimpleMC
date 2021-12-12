SimpleMC
==========

A code for cosmological parameter estimation where only
expansion history matters. `Simplemc` contains several tools 
to perform parameter estimation and maximum likelihood estimation. 

This code is not intendent as a replacement of CosmoMC. It simply
reimplements the homogeneus
part of cosmomc code. This allows one to
fit BAO data, SNeIa, cosmic chronometers and fsigma8 among others. 

Any crazy model can be implemented without much hassle and very quickly. 
One just needs to define
how Hubble parameter varies with redshift. It might also be useful for
pedagogical reasons.

Temporary documentation in: https://igomezv.github.io/SimpleMC/ 

You can install simplemc on your computer using:

     $ git clone https://github.com/ja-vazquez/SimpleMC 
          
     $ cd SimpleMC
     
     $ pip3 install -e .

then you can delete the cloned repo because you must have `simplemc` installed locally. 

Other way to install `simplemc` (without clonning) is:

     $ pip3 install -e git+https://github.com/ja-vazquez/SimpleMC#egg=simplemc


In the `requirements.txt` file are the basic libraries to run `simplemc`, but some functions such as graphics or neural networks may not be available. To get the full requirements use `requirements_full.txt`. You can manually install these dependencies with `pip3 install -r requirements_full.txt`.



Quick start:
------------

A quick start with the code, try:

`python3 test.py` 


You can set the `baseConfig.ini` with the analyzer, model and datasets that
you need. 


Directory structure:
--------------------
simplemc (root) - actual executables - you should read them like scripts

analyzers - analysis code (samplers and optimizers)

data - data files 

cosmo - Cosmoligical issues

likelihoods - Several likelihoods for several datasets

models - Cosmological models based on LCDM 

tools - various tools we use for plotting, etc. 



How does it work:
-----------------

It was supposed to be minimal, allowing any one to add models, without
being spaghetti.

There are three basic kind of objects:

- Cosmology objects that define models. These are based on
BaseCosmology.py 

- Likelihood objects that calculate Likelihood given some theory
predictions from Cosmoloy. These are based on BaseLikelihood. You tell
them the theory you want to use with setTheory and then you can ask
them about loglikelihood(). They can also tell you which Parameters
they have with freeParameters. These are based on Parameter.py and are
really simple -- name defines the parameter.

- Analyzers -- these take likelihoods and do something with it (see below).

To see how data is used, review Likelihood folder. 

To see how models are created, have a look, for example at
LCDMCosmology.py and oLCDMCosmology.py.

Parameters and priors are defined in paramDefs.py

LCDMCosmology inherits BaseCosmology and defines three functions:

- freeParameters - here we announced what are the parameters in this
model

- updateParameters - here we process a request to update
parameters. Here we store parameters, but must also update rs of the
BaseCosmogy if it is affected. At the moment we keep rs constant.

- HSSquared_a - this is the relative Hubble parameter squared, 
              in simple ternms H(a)**2/H0**2

BaseCosmology deals with integrating the above to various distance
moduli. It also handles the prefactors of the form c/(H0*rs). We can
either have h=H0/100 as a parameter and then believe that rs is
whatever it is, or have the entire phenomenological prefactor
c/(H0*rs) as a free parameter. One can switch this on as demonstrated
in Run/TestRun.py with T.setVaryPrefactor().

RadiationAndNeutrinos.py deals with, well, radiation and neutrino components.
NuDensity.py calculates the annoying neutrino integral.


Analyzers:
----------

So, one can then create a composite Likelihood by saying e.g.

L=CompositeLikelihood([
    DR11LOWZ(),
    DR11CMASS(),
    DR11LyaAuto(),
    DR11LyaCross()
    ])

and a then connect it to a theory, like 
L.setTheory(oLCDMCosmology())

This is now a nice package:
L.freeParameters() will tell you which parameters you can jiggle.
L.updateParams() allows you to update them.
L.loglike() returns the total log likelihood.

Different Analyzers can use this interface to  do different things
with this likelihood, without knowing anything about the underlying
theory. There are two attached:

- MaxLikeAnalyzer - finds the maximum and then uses second derivative to
get errors, only that that doesn't quite work and I think this is due
to linear interpolation in the chi2 tables

- MCMCAnalyzer - a vanilla MCMC samples compatible with CosmoMC format

- Nested Sampling - through a modified version of the dynesty library. 
It calculates Bayesian evidence and allows model comparison.

- Genetic algorithms - finds the maximum with the genetic algorithms of the deap library.

- emcee - emcee library.

- mcevidence - approximate Bayesian evidence to emcee and MCMC. 



