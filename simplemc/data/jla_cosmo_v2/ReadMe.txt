	       __________________________________________

		SDSS-II/SNLS3 JOINT LIGHT-CURVE ANALYSIS

			      Marc Betoule
	       __________________________________________


Table of Contents
_________________

1 Installation of the cosmomc plugin
2 SALT2 model
3 Error propagation
.. 3.1 Error decomposition
.. 3.2 SALT2 light-curve model uncertainties
4 Recalibrated light-curves in native SALT2 format
5 Presentation material


This page contains links to data associated with the SDSS-II/SNLS3 Joint
Light-Curve Analysis ([Betoule et al. 2014], submitted to A&A).

Light-curve fit parameters and associated uncertainties for the JLA
SDSS-II/SNLS3 sample are released in the form of a plugin for the
[cosmomc] software. We also deliver intermediate products and
presentation material.

Download the archive: [jla_cosmo_v2.tgz].


[Betoule et al. 2014] http://arxiv.org/abs/1401.4064

[cosmomc] http://cosmologist.info/cosmomc/

[jla_cosmo_v2.tgz] ./jla_cosmo_v2.tgz


1 Installation of the cosmomc plugin
====================================

  Untar the archive and copy the following files to the corresponding
  cosmomc directory:

  - source/supernovae_JLA.f90
  - data/jla.dataset
  - data/JLA.paramnames
  - data/jla_*

  Add the JLA likelihood in source/supernovae.f90:

  ,----
  | use JLA
  | ...
  | call JLALikelihood_Add(LikeList, Ini)
  `----

  Modify the "SUPERNOVAE =" line of the cosmomc source/Makefile line to
  compile the JLA likelihood:

  ,----
  | SUPERNOVAE = supernovae_Union2.o supernovae_SNLS.o supernovae_JLA.o
  `----

  Depending on your cosmomc installation, you may obtain the same result
  by extracting directly the archive in your cosmomc directory.

  Light-curve parameters are given in file `data/jla_lcparams.txt'.  The
  columns are:
  ,----
  | name: name of the SN
  | zcmb: CMB frame redshift (including peculiar velocity corrections for
  |       nearby supernova based on the models of M.J. Hudson)
  | zhel: Heliocentric redshift (note both zcmb and zhel are needed
  |       to compute the luminosity distance)
  | dz: redshift error (no longer used by the plugin)
  | mb: B band peak magnitude
  | dmb: Error in mb
  | x1: SALT2 shape parameter
  | dx1: Error in shape parameter
  | colour: Colour parameter
  | dcolour: Error in colour
  | 3rdvar: In these files, the log_10 host stellar mass
  | d3rdvar: Error in 3rdvar
  | cov_m_s: The covariance between mb and x1
  | cov_m_c: The covariance between mb and colour
  | cov_s_c: The covariance between x1 and colour
  | set: A number indicating which sample this SN belongs to, with
  |    0 - SNLS, 1 - lowz, 2 - Riess HST 3 - SDSS
  `----


2 SALT2 model
=============

  Version 2.4 of the SALT2 code, is available from the [SALT2 home
  page].  The corresponding data, including the model trained on the JLA
  training sample, and the description of instruments involved in the
  cosmology sample can be downloaded here: [salt2-4_data.tgz]
  (md5:af829c7f177c6037adcbe0a89360464a).

  With this model, the fit of the lightcurves of SNLS 03D4ag (see
  Recalibrated light-curves in native SALT2 format) should return:
  ,----
  | BEGIN_OF_FITPARAMS Salt2Model
  | DayMax 52830.9312582 0.0955792582608 
  | Redshift 0.285 0 F 
  | Color -0.0851965244364 0.0234763445429 
  | X0 5.65786629794e-05 1.52744999131e-06 
  | X1 0.937399343815 0.104657850633 
  | RestFrameMag_0_B 21.2581394711 0.0292572971214 
  | CovColorColor 0.000551138753097 -1 
  | CovColorDayMax -0.00012040087613 -1 
  | CovColorRestFrameMag_0_B 0.000524478290388 -1 
  | CovColorX0 -2.73917683701e-08 -1 
  | CovColorX1 -0.000243030505277 -1 
  | CovDayMaxDayMax 0.00913539460969 -1 
  | CovDayMaxRestFrameMag_0_B 0.000131399790871 -1 
  | CovDayMaxX0 -6.85109541277e-09 -1 
  | CovDayMaxX1 0.00179671536042 -1 
  | CovRestFrameMag_0_BX0 -4.46890539999e-08 -1 
  | CovRestFrameMag_0_BX1 0.00100854829146 -1 
  | CovX0X0 2.33310347595e-12 -1 
  | CovX0X1 -5.26149351948e-08 -1 
  | CovX1X1 0.0109532656992 -1 
  | END_OF_FITPARAMS Salt2Model
  `----
  Note: This value of msb is different from the value published in
  Betoule et al. 2014 by -0.0012 because of malmquist bias correction.


  [SALT2 home page] http://supernovae.in2p3.fr/salt/

  [salt2-4_data.tgz] ./salt2-4_data.tgz


3 Error propagation
===================

3.1 Error decomposition
~~~~~~~~~~~~~~~~~~~~~~~

  We deliver separate covariances matrices following the decomposition
  of the error on light-curve parameters proposed in [Betoule et
  al. 2014], Eq. 11. The ordering of the \(3 \times N_{\rm SN} = 2256\)
  vector of light-curve parameters is: \[ \vec \eta = (m^\star_1,
  {X_1}_1, C_1, \cdots, m^\star_{N_{\rm SN}}, {X_1}_{N_{\rm SN}},
  C_{N_{\rm SN}}).  \] The covariance matrix of this vector can be
  assembled as: \[ C_{\vec \eta} = C_{\rm stat} + C_{\rm cal} + C_{\rm
  model} + C_{\rm bias} + C_{\rm host} + C_{\rm dust} + C_{\rm pecvel} +
  C_{\rm nonIa} \] with:
  - \(C_{\rm stat} \): statistical uncertainty (including the
    statistical uncertainty of the SALT2 Model, but not including the
    $\sigma_{\rm coh}$ and $\sigma_{lens}$ of Eq. 13)
  - \(C_{\rm cal} \): calibration uncertainty
  - \(C_{\rm model} \): systematic uncertainty on the model (See [Mosher
    et al. 2014])
  - \(C_{\rm bias} \): uncertainty on the bias correction
  - \(C_{\rm host} \): uncertainty on the fonctionnal form of the
    luminosity-host mass relation.
  - \(C_{\rm dust} \): uncertainty on the Milky Way dust column density
  - \(C_{\rm pecvel}\): uncertainty on the peculiar velocity correction
    (Systematic only: does not include the \(\sigma_z\) term of Eq. 13.)
  - \(C_{\rm nonIa} \): potential contamination by non-Ia.

  The \(3N_{\rm SN}\times 3N_{\rm SN} = 2256\times2256\) matrices are
  provided in [covmat_v1.tgz].


  [Betoule et al. 2014] http://arxiv.org/abs/1401.4064

  [Mosher et al. 2014] http://arxiv.org/abs/1401.4065

  [covmat_v1.tgz] ./covmat_v1.tgz


3.2 SALT2 light-curve model uncertainties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  We plan to release shortly material to enable easy propogation of
  calibration and statistical uncertainties affecting the SALT2 model to
  the fit of other SNe samples.


4 Recalibrated light-curves in native SALT2 format
==================================================

  Light curves of SN includes in the cosmology sample:
  [jla_light_curves.tgz] (md5:14f402041003661d284ac18dbcc7ec24).

  It includes:
  - the recalibrated SNLS3 light curve (lc-??D??.list)
  - the recalibrated SDSS-II light-curves (lc-SDSS*.list, subsample of
      [http://sdssdp62.fnal.gov/sdsssn/DataRelease/index.html])
  - Nearby and HST light-curves from the C11 compilation (lc-sn*.list
    and lc-name.list, available from
    [https://tspace.library.utoronto.ca/snls/])

  Six SNe (sn2005hc, sn2005ir, sn2006nz, sn2006oa, sn2006ob, sn2006on)
  are in common with SDSS and either CfAIII or CSP. The light-curve file
  contains a merge of all available photometry.


  [jla_light_curves.tgz] ./jla_light_curves.tgz


5 Presentation material
=======================

  [./hd.png]

  [./althd.png]

  [./contour_omegam_w.png]

  [./contour_w_wa.png]

  [./diff_mu.png]
