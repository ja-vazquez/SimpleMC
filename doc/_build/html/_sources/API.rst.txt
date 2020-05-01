===
API
===

This page details the methods and classes provided by the `SuperMC`.

**Run modules**
===================


.. driver
.. -------
.. .. py:class:: Run.driver


DriverMC
----------------
.. py:class:: Run.DriverMC


    

PostProcessing
---------------
.. automodule:: Run.PostProcessing
    :members:
    :undoc-members:
    :show-inheritance:

RunBase
----------------
.. py:class:: Run.RunBase
    


wqdriver
--------------

.. automodule:: Run.wqdriver
    :members:
    :undoc-members:
    :show-inheritance:



Analizers
=========


MCMCAnalyzer
---------------
.. automodule:: py.MCMCAnalyzer
    :members:
    :undoc-members:
    :show-inheritance:


CosmoMCImportanceSampler
-------------------------
.. automodule:: py.CosmoMCImportanceSampler
    :members:
    :undoc-members:
    :show-inheritance:


GelmanRubinDiagnostic
-----------------------
.. automodule:: py.GelmanRubinDiagnostic
    :members:
    :undoc-members:
    :show-inheritance:

PlotterMC
----------------
.. automodule:: py.PlotterMC
    :members:
    :undoc-members:
    :show-inheritance:


MaxLikeAnalyzer
-----------------
.. automodule:: Analizers.MaxLikeAnalyzer
    :members:
    :undoc-members:
    :show-inheritance:

Individuo
----------------
.. automodule:: Analizers.Individuo
    :members:
    :undoc-members:
    :private-members:
    :show-inheritance:

Poblacion
----------------
.. automodule:: Analizers.Poblacion
    :members:
    :undoc-members:
    :private-members:
    :show-inheritance:

GeneticMedel
----------------
.. automodule:: Analizers.GeneticMedel
    :members:
    :undoc-members:
    :private-members:
    :show-inheritance:


**Models**
===========

BaseCosmology class
---------------------

This class is the general structure of the models. It have the following methods:
   
.. note::

   Base Cosmology class doesn't know about your parameterization of the equation of state or densities or anything. However, it does know about Hubble's constant at z=0 OR the prefactor c/(H0*rd) which should be fit for in the case of "rd agnostic" fits. That is why you should let it declare those parameterd based on its settings However, to get the angular diameter distance you need to pass it its Curvature parameter (Omega k basically), so you need to update it. 
  

Models.LCDMCosmology
---------------------
.. py:class:: Models.LCDMCosmology
.. function:: LCDMCosmology.setNoObh2prior()


.. note::

   This is LCDM cosmology. It is used as a base class for most other cosmologies, mostly because it treats Neutrinos and Radiation hassle.


**Likelihoods**
================

.. note::

   * **BaseLikelihood** is the basic structure of the likelihoods objects.
   * **CompositeLikelihood** and **LikelihoodMultiplier** allow to create more complex probabilities from the combination of two or more.

Likelihoods.BaseLikelihood
---------------------------

.. automodule:: Likelihoods.BaseLikelihood
    :members:
    :undoc-members:
    :private-members:
    :show-inheritance:

CompositeLikelihood
--------------------
.. automodule:: Likelihoods.CompositeLikelihood
    :members:
    :undoc-members:
    :show-inheritance:

LikelihoodMultiplier
---------------------
.. automodule:: Likelihoods.LikelihoodMultiplier
    :members:
    :undoc-members:
    :show-inheritance:


Examples:
----------

CompressedSNLikelihood
------------------------
.. automodule:: Likelihoods.CompressedSNLikelihood
    :members:
    :undoc-members:
    :show-inheritance:

CompressedGenericLikelihood
-------------------------------

.. automodule:: Likelihoods.CompressedGenericLikelihood
    :members:
    :undoc-members:
    :show-inheritance:


.. Likelihoods.BAOLikelihoods
.. ----------------------------

.. .. automodule:: Likelihoods.BAOLikelihoods
..     :members:
..     :undoc-members:
..     :show-inheritance:


Others modules with Likelihoods
----------------------------------

* BAOLikelihoods
* GaussBAODVLikelihood
* TabulatedBAOLikelihood
* TabulatedBAODVLikelihood
* HubbleParameterLikelihood
* SimpleCMB
* WangWangCMB
* CompressedHDLikelihood




**Cosmo modules**
===================


BaseCosmology
--------------
.. automodule:: Cosmo.BaseCosmology
    :members:
    :undoc-members:
    :show-inheritance:



CosmoApprox
----------------
.. automodule:: Cosmo.CosmoApprox
    :members:
    :undoc-members:
    :show-inheritance:


.. NuDensity
.. ---------------
.. .. automodule:: Cosmo.NuDensity
..     :members:
..     :undoc-members:
..     :show-inheritance:

Parameter
----------------
.. automodule:: Cosmo.Parameter
    :members:
    :undoc-members:
    :show-inheritance:


ParamDefs
--------------

.. automodule:: Cosmo.ParamDefs
    :members:
    :undoc-members:
    :show-inheritance:

.. RadiationAndNeutrinos
.. ----------------------

.. .. automodule:: Cosmo.RadiationAndNeutrinos
..     :members:
..     :undoc-members:
..     :show-inheritance:


**pybambi modules**
===================


bambi
-------
.. automodule:: pybambi.bambi
    :members:
    :undoc-members:
    :show-inheritance:



base
----------------
.. automodule:: pybambi.base
    :members:
    :undoc-members:
    :show-inheritance:


pybambimanager
---------------
.. automodule:: pybambi.pybambimanager
    :members:
    :undoc-members:
    :show-inheritance:

kerasnet
----------------
.. automodule:: pybambi.kerasnet
    :members:
    :undoc-members:
    :show-inheritance:











