april
=====

A simple MCMC code for cosmological parameter estimation where only
expansion history matters. Written by Anže Slosar and Jose Vazquez.

This is the github version, based off internal BOSS svn. It should
work, but feel free to complain to me (Anže Slosar).

This code is not intendent as a replacement of CosmoMC. It simply
reimplements the homogeneus part of cosmomc code. This allows one to
fit BAO data, etc. in the spirit of arxiv:xxxx for any bat-shit crazy
model without much hassle and very quickly. One just needs to define
how Hubble parameter varies with redshift. It might also be useful for
pedagogical reasons.


Directory structure:
--------------------

py - analysis code  
Run - actual executables - you should read them like scripts  
data - data files  
chains - empty directory to store chains  
example - a simple example to see how the code works  
tools - various tools we use for plotting, etc.  
attick - various hacks and tests  

Quick start:
------------

A quick start with the code, try:

`examples/TestRun.py`  
and then, after a few minutes  
`examples/TestPlot.py`  


But this is now how we normally run the code. To fit for LCDM using
BOSS BAO without assuming you know what rd is, say something

`Run/driver.py pre LCDM BBAO`

To fit for oLCDM using Lyman-alpha BAO and LyaCross BAO +Planck say 

`Run/driver.py phy oLCDM LBAO+Planck`

`Run/driver.py` will print all possible options


How does it work:
-----------------

It was supposed to be minimal, allowing any one to add models, without
being spaghetti.

There are three basic kind of objects:

Cosmology objects that define models. These are based on
BaseCosmology.py 

Likelihood objects that calculate Likelihood given some theory
predictions from Cosmoloy. These are based on BaseLikelihood. You tell
them the theory you want to use with setTheory and then you can ask
them about loglikelihood(). They can also tell you which Parameters
they have with freeParameters. These are based on Parameter.py and are
really simple -- name defines the parameter (like white socks define
an american).

Analyzers -- these take likelihoods and do something with it (see below)

To see how data is used, start with BOSSLikelihoods.py This defines
objects that actually correspond to real likelihoods.  These are then
chucked into one single CompositeLikelihood as seen in Run/TestRun.py.

To see how models are created, have a look, for example at
LCDMCosmology.py and oLCDMCosmology.py.

Parameters and priors are defined in ParamDefs.py

LCDMCosmology inherits BaseCosmology and defines three functions:

freeParameters - here we announced what are the parameters in this
model

updateParameters - here we process a request to update
parameters. Here we store parameters, but must also update rs of the
BaseCosmogy if it is affected. At the moment we keep rs constant.

HSSquared_a - this is the relative Hubble parameter squared, 
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

MaxLikeAnalyzer - finds the maximum and then uses second derivative to
get errors, only that that doesn't quite work and I think this is due
to linear interpolation in the chi2 tables

MCMCAnalyzer - a vanilla MCMC samples compatible with CosmoMC format

tools/cosmich.py is an archaic ugly script that I use for plotting.




