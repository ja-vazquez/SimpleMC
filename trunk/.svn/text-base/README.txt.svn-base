Welcome!

Directory structure:
--------------------
py - analysis code
Run - actual executables - you should read them like scripts
data - data files
chains - empty directory to store chains

Getting existing chains:
------------------------

Run Pack/unpack.py

You can check some existing plotting routines (they are really quick 
hacks, but to get the idea) in Run/plot_nice.py.

Best fits can be found in maxlike files. For example:

anze@basquiat:~/work/April$ cat chains_140401_155316M/owaCDM_phy_BBAO+SN+CMBP.maxlike
5 22.0976 0.300376 0.0226451 0.679125 -0.89428 -0.510333 -0.00522912 -0.0546475 -1.59601 -1.48502 -2.28881 -0.184969 -16.4881 0

This means that the best fit point has weight 5 in the chain and -loglike
of 22.0976. The rest are actual parameters, whose meaning you can get:  

anze@basquiat:~/work/April$ cat chains_140401_155316M/owaCDM_phy_BBAO+SN+CMBP.params
Om
Obh2
h
w
wa
Ok
DR11LOWZ_like
DR11CMASS_like
DR11LyaAuto_like
DR11LyaCross_like
S_Planck_like
BetouleSN_like
theory_prior



Quick start:
------------

To fit for LCDM using BOSS BAO without assuming you know what rd is, say
something

Run/driver.py pre LCDM BBAO


To fit for oLCDM using Lyman-alpha BAO and LyaCross BAO +Planck say 

Run/driver.py phy oLCDM LBAO+CMBP

Run/driver.py will print all possible options


How does it work:
-----------------

It was supposed to be minimal, allowing any one to add models, without
being spaghetti.

To see how data is used, start with BOSSLikelihoods.py
This defines objects that actually correspond to real likelihoods.
These are then chucked into one single CompositeLikelihood as seen
in Run/TestRun.py

To see how models are created, have a look at LCDMCosmology.py and
oLCDMCosmology.py.

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

Analyzers:
----------

So, one can then create a composite Likelihood by saying e.g.

L=CompositeLikelihood([
    DR11LOWZ(),
    DR11CMASS(),
    DR11LyaAuto(),
    DR11LyaCross()
    ])

and then connect it to a theory, like 
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
MCMCAnalyzer - rudimentary MCMC - could be improved, but really kinda
works for the simple problems we have here... Output is compatible
with cosmoMC output (in the sense that the format is  weight -logLike
parameter values). Run/cosmich is an archaic ugly script that I use
for plotting. Run/TestPlot uses that, but this is just a quick hack at
the moment.


