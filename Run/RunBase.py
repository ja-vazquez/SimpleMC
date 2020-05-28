# coding=utf-8
# Add paths, we want to be able to run in either root or Run/

import sys
sys.path = ["py", "../py", "Models", "Cosmo", "Likelihoods"] + sys.path

#TODO -- Include models used in several papers
#TODO -- Add Planck 15
#TODO -- Add DR12 Galaxies
#TODO -- Add Compress Pantheon

# Cosmologies already included
from LCDMCosmology import LCDMCosmology
from oLCDMCosmology import oLCDMCosmology
from wCDMCosmology import wCDMCosmology
from owa0CDMCosmology import owa0CDMCosmology
from PolyCDMCosmology import PolyCDMCosmology
from JordiCDMCosmology import JordiCDMCosmology
from WeirdCDMCosmology import WeirdCDMCosmology
from TiredLightDecorator import TiredLightDecorator
from SplineLCDMCosmology import SplineLCDMCosmology
from DecayLCDMCosmology import DecayLCDMCosmology
from StepCDMCosmology import StepCDMCosmology
from EarlyDECosmology import EarlyDECosmology
from SlowRDECosmology import SlowRDECosmology
from BinnedWCosmology import BinnedWCosmology
from PhiCDMCosmology import PhiCosmology
#from STCDMCosmology import STCDMCosmology

#from DGPCDMCosmology import DGPCDMCosmology

#Generic model
from GenericCosmology import GenericCosmology
from GenericPantheon import GenericPantheon

# Composite Likelihood
from CompositeLikelihood import CompositeLikelihood

# Likelihood Multiplier
from LikelihoodMultiplier import LikelihoodMultiplier

# Likelihood modules
from BAOLikelihoods import DR11LOWZ, DR11CMASS, DR14LyaAuto, DR14LyaCross, \
        SixdFGS, SDSSMGS, DR11LyaAuto, DR11LyaCross, eBOSS, DR12Consensus
from SimpleCMB import PlanckLikelihood, PlanckLikelihood_15, WMAP9Likelihood
from CompressedSNLikelihood    import BetouleSN, UnionSN, BinnedPantheon
from PantheonSNLikelihood      import PantheonSNLikelihood
from HubbleParameterLikelihood import RiessH0
from CompressedHDLikelihood    import HubbleDiagram
from Compressedfs8Likelihood import fs8Diagram

from GenericLikelihood import StraightLine
from GenericPantheonSNLikelihood import GenericPantheonSNLikelihood


#Importance Sampling
from CosmoMCImportanceSampler import *


# String parser Aux routines
model_list = "LCDOM, LCDMasslessnu, nuLCDM, NeffLCDM, noradLCDM, nuoLCDM, nuwLCDM, oLCDM, wCDM, waCDM, owCDM,"\
    "owaCDM, JordiCDM, WeirdCDM, TLight, StepCDM, Spline, PolyCDM, fPolyCDM, Decay, Decay01, Decay05,"\
    "EarlyDE, EarlyDE_rd_DE, SlowRDE"


def ParseModel(model):
    """ 
    Parameters
    -----------
    model:
         name of the model, i.e. LCDM

    Returns
    -----------
    object - info/calculations based on this model: i.e. d_L, d_A, d_H

    """

    if model == "LCDM":
        T = LCDMCosmology()
    elif model == "LCDMmasslessnu":
        T = LCDMCosmology(mnu=0)
    elif model == "nuLCDM":
        T = LCDMCosmology()
        T.setVaryMnu()
    elif model == "NeffLCDM":
        LCDMCosmology.rd_approx = "CuestaNeff"
        T = LCDMCosmology()
        T.setVaryNnu()
    elif model == "noradLCDM":
        T = LCDMCosmology(disable_radiation=True)
    elif model == "oLCDM":
        T = oLCDMCosmology()
    elif model == "nuoLCDM":
        T = oLCDMCosmology()
        T.setVaryMnu()
    elif model == "wCDM":
        T = wCDMCosmology()
    elif model == "nuwCDM":
        T = wCDMCosmology()
        T.setVaryMnu()
    elif model == "waCDM":
        T = owa0CDMCosmology(varyOk=False)
    elif model == "owCDM":
        T = owa0CDMCosmology(varywa=False)
    elif model == "owaCDM":
        T = owa0CDMCosmology()
    elif model == "JordiCDM":
        T = JordiCDMCosmology()
    elif model == "WeirdCDM":
        T = WeirdCDMCosmology()
    elif model == "TLight":
        T = TiredLightDecorator(PolyCDMCosmology())
    elif model == "StepCDM":
        T = StepCDMCosmology()
    elif model == "Spline":
        T = SplineLCDMCosmology()
    elif model == "DecayFrac":
        T = DecayLCDMCosmology()
    elif model == "Decay":
        T = DecayLCDMCosmology(varyxfrac=False, xfrac=1.0)
    elif model == "Decay01":
        T = DecayLCDMCosmology(varyxfrac=False, xfrac=0.1)
    elif model == "Decay05":
        T = DecayLCDMCosmology(varyxfrac=False, xfrac=0.5)
    elif model == "PolyCDM":
        T = PolyCDMCosmology()
    elif model == "fPolyCDM":
        T = PolyCDMCosmology(polyvary=['Om1', 'Om2'])
    elif model == "PolyOk": ## polycdm for OK
        T = PolyCDMCosmology(Ok_prior=10.)
    elif model == "PolyOkc": ## polycdm sans Om2 term to couple two
        T = PolyCDMCosmology(polyvary=['Om1', 'Ok'], Ok_prior=10.)
    elif model == "PolyOkf": ## polycdm sans Om2 term to couple two
        T = PolyCDMCosmology(polyvary=['Om1', 'Om2'])
    elif model == "EarlyDE":
        T = EarlyDECosmology(varyw=False, userd_DE=False)
    elif model == "EarlyDE_rd_DE":
        T = EarlyDECosmology(varyw=False)
    elif model == "SlowRDE":
        T = SlowRDECosmology(varyOk=False)
    elif model == "Binned":
        T = BinnedWCosmology()
    elif model == 'wDM':
        T = wCDMCosmology()
    elif model == 'sline':
        T = GenericCosmology()
    elif model == 'GPantheon':
        T = GenericPantheon()
    #elif model == 'DGP':
    #    T = DGPCDMCosmology()
    elif model == "Phi_exp_p0":
        T = PhiCosmology(mu=0, alpha=1, varybeta=True)
    elif model == "Phi_pow_test_i":
        T = PhiCosmology(beta=0, varymu=True, varyilam=True)
    elif model == "Phi_exp_pow2":
        T = PhiCosmology(mu=0, alpha=2, varybeta=True, varyilam=True)
    elif model == "Phi_pow_exp":
        T = PhiCosmology(alpha=1, varybeta=True, varymu=True, varyilam=True)
    elif model == "Phi_exp_pow_a":
        T = PhiCosmology(mu=0, varybeta=True, varyalpha=True, varyilam=True)
    elif model == "Phi_pow2_exp_pow2":
        T = PhiCosmology(mu=2, alpha=2, varybeta=True, varyilam=True)
    elif model == "Phi_cosh":
        T = PhiCosmology(beta=0, mu=-1, varyalpha=True, varyilam=True)
    elif model == "Phi_cosh_1":
        T = PhiCosmology(beta=-1, mu=-1, varyalpha=True, varyilam=True)
    elif model == "Phi_cos_1":
        T = PhiCosmology(beta=1, mu=-1, varyalpha=True, varyilam=True)
    #elif model == 'ST':
    #    T = STCDMCosmology()
    else:
        print("Cannot recognize model", model)
        sys.exit(1)

    return T


data_list = "BBAO, GBAO, GBAO_no6dF, CMASS, LBAO, LaBAO, LxBAO, MGS, Planck, WMAP, PlRd, WRd, PlDa, PlRdx10,"\
    "CMBW, SN, SNx10, UnionSN, RiessH0, 6dFGS"


def ParseDataset(datasets):
    """ 
    Parameters
    -----------
    datasets:
         name of datasets, i.e. BBAO

    Returns
    -----------
    object - likelihood

    """
    dlist = datasets.split('+')
    L = CompositeLikelihood([])
    for name in dlist:
        if name == 'BBAO':
            L.addLikelihoods([
                DR11LOWZ(),
                DR11CMASS(),
                DR11LyaAuto(),
                DR11LyaCross(),
                SixdFGS(),
                SDSSMGS()
            ])
        elif name == 'GBAO11':
            L.addLikelihoods([
                DR11LOWZ(),
                DR11CMASS(),
                SixdFGS(),
                SDSSMGS()
            ])
        elif name == 'CBAO':
            L.addLikelihoods([
                DR12Consensus(),
                DR14LyaAuto(),
                DR14LyaCross(),
                SixdFGS(),
                SDSSMGS(),
                eBOSS()
            ])
        elif name == 'GBAOx10':
            L.addLikelihoods([
                LikelihoodMultiplier(DR11LOWZ(),  100.0),
                LikelihoodMultiplier(DR11CMASS(), 100.0),
                LikelihoodMultiplier(SixdFGS(),   100.0)
            ])
        elif name == 'GBAO_no6dF':
            L.addLikelihoods([
                DR11LOWZ(),
                DR11CMASS()
            ])
        elif name == 'CMASS':
            L.addLikelihoods([
                DR11CMASS()
            ])
        elif name == 'LBAO':
            L.addLikelihoods([
                DR14LyaAuto(),
                DR14LyaCross()
            ])
        elif name == 'LBAO11':
            L.addLikelihoods([
                DR11LyaAuto(),
                DR11LyaCross()
            ])
        elif name == 'LaBAO':
            L.addLikelihoods([
                DR14LyaAuto(),
            ])
        elif name == 'LxBAO':
            L.addLikelihoods([
                DR14LyaCross(),
            ])
        elif name == "MGS":
            L.addLikelihood(SDSSMGS())
        elif name == '6dFGS':
            L.addLikelihood(SixdFGS())
        elif name == 'eBOSS':
            L.addLikelihood(eBOSS())
        elif name == 'Planck':
            L.addLikelihood(PlanckLikelihood())
        elif name == 'Planck_15':
            L.addLikelihood(PlanckLikelihood_15())
        elif name == 'WMAP':
            L.addLikelihood(WMAP9Likelihood())
        elif name == 'PlRd':
            L.addLikelihood(PlanckLikelihood(kill_Da=True))
        elif name == 'WRd':
            L.addLikelihood(WMAP9Likelihood(kill_Da=True))
        elif name == 'PlDa':
            L.addLikelihood(PlanckLikelihood(kill_rd=True))
        elif name == 'PlRdx10':
            L.addLikelihood(LikelihoodMultiplier(
                PlanckLikelihood(kill_Da=True), 100.0))
        elif name == 'CMBW':
            L.addLikelihood(WMAP9Likelihood())
        elif name == 'Pantheon':
            L.addLikelihood(PantheonSNLikelihood())
        elif name == 'BPantheon_15':
            L.addLikelihood(BinnedPantheon())
        elif name == 'SN':
            L.addLikelihood(BetouleSN())
        elif name == 'SNx10':
            L.addLikelihood(LikelihoodMultiplier(BetouleSN(), 100.0))
        elif name == 'UnionSN':
            L.addLikelihood(UnionSN())
        elif name == 'RiessH0':
            L.addLikelihood(RiessH0())
        elif name == 'HD':
            L.addLikelihood(HubbleDiagram())
        elif name == 'fs8':
            L.addLikelihood(fs8Diagram())
        elif name == 'dline':
            L.addLikelihood(StraightLine())
        elif name == 'CPantheon_15':
            L.addLikelihood(GenericPantheonSNLikelihood())
        else:
            print("Cannot parse data, unrecognizable part:", name)
            sys.exit(1)

    return L
