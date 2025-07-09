# coding=utf-8
import sys

# Cosmologies already included
from .models import LCDMCosmology
from .models import oLCDMCosmology
from .models import wCDMCosmology
from .models import owa0CDMCosmology
from .models import PolyCDMCosmology
from .models import JordiCDMCosmology
from .models import WeirdCDMCosmology
from .models import TiredLightDecorator
from .models import DecayLCDMCosmology
from .models import EarlyDECosmology
from .models import SlowRDECosmology
from .models import DGPCDMCosmology
from .models import AnisotropicCosmology
from .models import GraduatedCosmology
from .models import QuintomCosmology
from .models import RotationCurves
from .models import ChaplyginCosmology
from .models import HolographicCosmology

#Non-parametric functions
from .models import ReconstructedCosmology
from .models import BinnedWCosmology
from .models import CompressPantheon
from .models import TanhCosmology

#Generic model
from .models.SimpleModel import SimpleModel, SimpleCosmoModel

# Composite Likelihood
from .likelihoods.CompositeLikelihood import CompositeLikelihood

# Likelihood Multiplier
from .likelihoods.LikelihoodMultiplier import LikelihoodMultiplier

# Likelihood modules
from .likelihoods.BAOLikelihoods import DR11LOWZ, DR11CMASS, DR14LyaAuto, DR14LyaCross, \
                                        SixdFGS, SDSSMGS, DR11LyaAuto, DR11LyaCross, eBOSS, \
                                        DR12Consensus, DR16BAO, DESIBAO, DESIDR2BAO
from .likelihoods.SimpleCMBLikelihood import PLK, PLK15, PLK18, WMAP9
from .likelihoods.CompressedSNLikelihood import BetouleSN, UnionSN
from .likelihoods.SNLikelihood import JLASN_Full
from .likelihoods.PantheonSNLikelihood import PantheonSN, BinnedPantheon
from .likelihoods.PantheonPlusSNLikelihood import PantheonPlus
from .likelihoods.UNION3Likelihood import UNION3
from .likelihoods.DESY5Likelihood import DESY5
from .likelihoods.CompressedHDLikelihood import HubbleDiagram, HD23
from .likelihoods.Compressedfs8Likelihood import fs8Diagram
from .likelihoods.HubbleParameterLikelihood import RiessH0, RiessH0_21
from .likelihoods.BBNLikelihood import BBN
from .likelihoods.StrongLensingLikelihood import StrongLensing
from .likelihoods.SimpleLikelihood import GenericLikelihood, StraightLine
from .likelihoods.RotationCurvesLikelihood import RotationCurvesLike

#Importance Sampling
#from .CosmoMCImportanceSampler import *


# String parser Aux routines
model_list = "LCDOM, LCDMasslessnu, nuLCDM, NeffLCDM, noradLCDM, nuoLCDM, nuwLCDM, oLCDM, wCDM, waCDM, owCDM,"\
    "owaCDM, JordiCDM, WeirdCDM, TLight, Spline, PolyCDM, fPolyCDM, Decay, Decay01, Decay05,"\
    "EarlyDE, EarlyDE_rd_DE, SlowRDE"


def ParseModel(model, **kwargs):
    """ 
    Parameters
    -----------
    model:
         name of the model, i.e. LCDM

    Returns
    -----------
    object - info/calculations based on this model: i.e. d_L, d_A, d_H

    """
    custom_parameters = kwargs.pop('custom_parameters', None)
    custom_function = kwargs.pop('custom_function', None)

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
    elif model == "NumuLCDM":
        LCDMCosmology.rd_approx = "CuestaNeff"
        T = LCDMCosmology()
        T.setVaryNnu()
        T.setVaryMnu()
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
    elif model == "Recon":
        T = ReconstructedCosmology()
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
    elif model == "Anisotropic":
        T = AnisotropicCosmology(varybd=False)
        LCDMCosmology.rd_approx = "CuestaNeff"
        T.setVaryNnu()
        T.setVaryMnu()
    elif model == "Binned":
        T = BinnedWCosmology()
    elif model == "Tanh":
        T = TanhCosmology()
    elif model == 'CPantheon':
        T = CompressPantheon()
    elif model == 'DGP':
        T = DGPCDMCosmology()
    elif model == 'Grad_Ok':
        T = GraduatedCosmology(varyOk=True)
    elif model == 'Quintess':
        T = QuintomCosmology(vary_mquin=True)
    elif model == 'Phantom':
        T = QuintomCosmology(vary_mphan=True)
    elif model == 'Quintom':
        T = QuintomCosmology(vary_mquin=True, vary_mphan=True)
    elif model == 'Quintom_couple':
        T = QuintomCosmology(vary_mquin=True, vary_mphan=True, vary_coupling=True)
    elif model == "Rotation":
        T = RotationCurves()
    elif model == "Chaplygin_s":
        T = ChaplyginCosmology(usesigmoid=True)
    elif model == "HDE_node1":
        T = HolographicCosmology(nodes=1)
    elif model == 'simple':
        T = SimpleModel(custom_parameters, custom_function)
    elif model == 'simple_cosmo':
        T = SimpleCosmoModel(custom_parameters, RHSquared=custom_function)
    else:
        print("Cannot recognize model", model)
        sys.exit(1)

    return T


data_list = "BBAO, GBAO, GBAO_no6dF, CMASS, LBAO, LaBAO, LxBAO, MGS, Planck, WMAP, PlRd, WRd, PlDa, PlRdx10,"\
    "CMBW, SN, SNx10, UnionSN, RiessH0, 6dFGS, PantheonPlus, DR16BAO, HD23, DESIBAO, DESIDR2BAO"


def ParseDataset(datasets, **kwargs):
    """ 
    Parameters
    -----------
    datasets:
         name of datasets, i.e. BBAO

    Returns
    -----------
    object - likelihood

    """
    path_to_data = kwargs.pop('path_to_data', None)
    path_to_cov = kwargs.pop('path_to_cov', None)
    fn = kwargs.pop('fn', 'generic')

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
        elif name == 'DR16BAO':
            L.addLikelihood(DR16BAO())
        elif name == 'DESI':
            L.addLikelihood(DESIBAO())
        elif name == 'DESIDR2':
            L.addLikelihood(DESIDR2BAO())  
        elif name == 'PLK':
            L.addLikelihood(PLK())
        elif name == 'PLK15':
            L.addLikelihood(PLK15())
        elif name == 'PLK18':
            L.addLikelihood(PLK18())
        elif name == 'PLKW':
            from .likelihoods.WangWangCMB import PlanckLikelihood
            L.addLikelihood(PlanckLikelihood())
        elif name == 'WMAP9':
            L.addLikelihood(WMAP9())
        elif name == 'PlRd':
            L.addLikelihood(PLK(kill_Da=True))
        elif name == 'WRd':
            L.addLikelihood(WMAP9(kill_Da=True))
        elif name == 'PlDa':
            L.addLikelihood(PLK(kill_rd=True))
        elif name == 'PlRdx10':
            L.addLikelihood(LikelihoodMultiplier(
                PLK(kill_Da=True), 100.0))
        elif name == 'Pantheon':
            L.addLikelihood(PantheonSN())
        elif name == 'BPantheon':
            L.addLikelihood(BinnedPantheon())
        elif name == 'PantheonPlus':
            L.addLikelihood(PantheonPlus())
        elif name == 'JLA':
            L.addLikelihood(JLASN_Full())
        elif name == 'Union3':
            L.addLikelihood(UNION3())
        elif name == 'DESY5':
            L.addLikelihood(DESY5())
        elif name == 'SN':
            L.addLikelihood(BetouleSN())
        elif name == 'SNx10':
            L.addLikelihood(LikelihoodMultiplier(BetouleSN(), 100.0))
        elif name == 'UnionSN':
            L.addLikelihood(UnionSN())
        elif name == 'RiessH0':
            L.addLikelihood(RiessH0())
        elif name == 'RiessH0_21':
            L.addLikelihood(RiessH0_21())
        elif name == 'HD':
            L.addLikelihood(HubbleDiagram())
        elif name == 'HD23':
            L.addLikelihood(HD23())
        elif name == 'BBN':
            L.addLikelihood(BBN())
        elif name == 'fs8':
            L.addLikelihood(fs8Diagram())
        elif name == 'dline':
            L.addLikelihood(StraightLine())
        #elif name == 'CPantheon_15':
        #    L.addLikelihood(PantheonLikelihood())
        elif name == 'RC':
            L.addLikelihood(RotationCurvesLike())
        elif name == "SL":
            L.addLikelihood(StrongLensing())
        elif name == 'generic':
            L.addLikelihood(GenericLikelihood(path_to_data=path_to_data,
                                              path_to_cov=path_to_cov,
                                              fn=fn))
        else:
            print("Cannot parse data, unrecognizable part:", name)
            sys.exit(1)

    return L
