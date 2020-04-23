==================
Introduction
==================

Requirements
-------------

This code runs both in Python 2x and 3x. However, we highly recommend Python 3x.

You need the following scientific modules:

.. code-block:: bash
   
   sudo pip install numpy matplotlib scipy dynesty


In addition, MCEvidence is necessary to estimate the bayesian evidence in the Metropolis-Hastings sampler:

.. code-block:: bash
   
   pip install git+https://github.com/yabebalFantaye/MCEvidence

For use Artificial Neural Networks with multinest and ellipsoidal sampling (as in pyBAMBI), you need to install:

.. code-block:: bash
   
   pip install tensorflow keras


If you want to use the full options to plot:

.. code-block:: bash
   
   pip install corner getdist

.. note:: All in one copy-paste line: 

   .. code-block:: bash
   
      pip install numpy matplotlib scipy nestle tensorflow keras corner getdist git+https://github.com/yabebalFantaye/MCEvidence



Quick Start
------------

In this section we show a basic use of SuperMC. 

1) First you need to create an *ini file* as follows:

.. code-block:: none

        [custom]
	chainsdir = chains
	model = LCDM
	prefact = py
	datasets = SN+Planck 
	;sampler can be {mcmc, nested, emcee}
	;or analyzers {MaxLike, genetic}
	sampler = mcmc

	[mcmc]
	nsamp = 5000 
	skip = 0 
	chainno = 1
	GRcriteria = 0.05
	temp = 2

	;it is recommended around nlivepoints=50*ndim
	[nested]
	;dynamic and neuralNetwork can be yes/no
	dynamic = no
	neuralNetwork = no
	;type: {'none','single','multi', 'balls', 'cubes'}
	nestedType = cubes	
	nlivepoints = 1024 
	accuracy = 0.5
	;u for flat(uniform) or g for gaussian prior
	priortype = u


	;data options: HD, BBAO, GBAO, GBAO_no6dF, CMASS, LBAO, LaBAO, LxBAO, MGS, Planck, WMAP, PlRd, WRd, PlDa, PlRdx10, CMBW, SN, SNx10, UnionSN, RiessH0, 6dFGS, generic

	;model options: LCDM, LCDMasslessnu, nuLCDM, NeffLCDM, noradLCDM, nuoLCDM, nuwLCDM, oLCDM, wCDM, waCDM, owCDM, owaCDM, JordiCDM, WeirdCDM, TLight, StepCDM, Spline, PolyCDM, fPolyCDM, Decay, Decay01, Decay05, EarlyDE, EarlyDE_rd_DE, SlowRDE, nled, generic

	;sampler options: mh, snest, mnest, bambi, sbambi, emcee, MaxLike
	;mcmc -> metropolis-hastings
	;none -> Prior mass without bounds
	;single -> Ellipsoidal nested sampling
	;multi -> multinest
	;balls ->  balls centered on each live point
	;cube -> cubes centered on each live point
	;MaxLikeAnalyzer -> Maximum Likelihood Analyzer
      ;genetic -> simple genetic algorithm for optimization


.. note::

   Considerations:
  
   * *prefact* and *nsamp* are only for Metropolis-Hastings.

   * *nlivepoints* and *accuracy* are only for nested sampling.

   * *sampler* options are:
   
      * mcmc : Metropolis-Hastings.
      * nested : Nested Sampling

   * *sampler* can be one *optimizer* of the following:
      
      * MaxLikeAnalyzer : from scipy.optimize.minimize
      * genetic : a Simple Genetic Algorithm
      

   * *skip* is burnin. 
  
   * For *priortype* u is uniform prior and g gaussian prior. At this time, only nested sampling accept both of them.
   
   * *chainsdir* is the directory where the chains in a text file and the plots will be saved.

2) Then in *driver.py*, please put the path of the *ini file*:

.. code-block:: python
   
   from DriverMC import DriverMC
   
   fileConfig = "/home/isidro/SuperMC_/baseConfig.ini"
   D = DriverMC(fileConfig)



3) For last, you can run in the *SuperMC* directory:

.. code-block:: bash
   
   $ python3 Run/driver.py

You can see your outputs in the *chains* directory.

* See the `plots <plotters.html>`_ .


General flow
-------------

.. figure:: /img/SuperMCDiagram.png

Samplers
---------

SuperMC contains the nest samplers:

   * Metropolis-Hastings.
      It is implemented in the **MCMCAnalyzer class**.

   * Ellipsoidal Nested Sampling.
      `Mukherjee, P., Parkinson, D., & Liddle, A. R. (2006). A nested sampling algorithm for cosmological model selection. The Astrophysical Journal Letters, 638(2), L51.  <https://iopscience.iop.org/article/10.1086/501068/metal>`_.

   * MULTINEST
      `Feroz, F., Hobson, M. P., & Bridges, M. (2009). MultiNest: an efficient and robust Bayesian inference tool for cosmology and particle physics. Monthly Notices of the Royal Astronomical Society, 398(4), 1601-1614. <https://academic.oup.com/mnras/article/398/4/1601/981502>`_.

   * BAMBI = MULTINEST + Artificial Neural Networks.
      `Graff, P., Feroz, F., Hobson, M. P., & Lasenby, A. (2012). BAMBI: blind accelerated multimodal Bayesian inference. Monthly Notices of the Royal Astronomical Society, 421(1), 169-180. <https://academic.oup.com/mnras/article/421/1/169/989639>`_.

   * Ellipsoidal Nested Sampling with Artificial Neural Networks


Also, for a previous quickly test, SuperMC have **MaxAnalyzer** that analyze the most probably values of the parameters.


Sampler comparison
-------------------

.. note:: 

   To verify the consistency of the parameter estimation among the different samplers available, we have made the following graph.

.. figure:: /img/samplersTriangle.png

     We estimates the posteriors of the parameters of the owaCDM model (dark energy with timedependent equation-of-state in a model of unknown curvature) using Supernovae type Ia, Cosmic Chronometers (Hubble Distance) and BAO .




