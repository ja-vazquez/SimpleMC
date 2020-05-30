__version__ = '1.9.0'
__author__  = 'JA Vazquez, I Gomez-Vargas, A Slosar'

import sys
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

if sys.version_info < (3, 0):
    logger.critical("You must use Python 3!")
    exit(1)

from simplemc.DriverMC import DriverMC
from simplemc.PostProcessing import PostProcessing

from simplemc.analyzers.MaxLikeAnalyzer import MaxLikeAnalyzer
from simplemc.analyzers import SimpleGenetic

from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.models.oLCDMCosmology import oLCDMCosmology
from simplemc.models.wCDMCosmology import wCDMCosmology
from simplemc.models.owa0CDMCosmology import owa0CDMCosmology
from simplemc.models.PolyCDMCosmology import PolyCDMCosmology
from simplemc.models.JordiCDMCosmology import JordiCDMCosmology
from simplemc.models.WeirdCDMCosmology import WeirdCDMCosmology
from simplemc.models.TiredLightDecorator import TiredLightDecorator
from simplemc.models.SplineLCDMCosmology import SplineLCDMCosmology
from simplemc.models.DecayLCDMCosmology import DecayLCDMCosmology
from simplemc.models.StepCDMCosmology import StepCDMCosmology
from simplemc.models.EarlyDECosmology import EarlyDECosmology
from simplemc.models.SlowRDECosmology import SlowRDECosmology
from simplemc.models.BinnedWCosmology import BinnedWCosmology
from simplemc.models.SureshCosmology import SureshCosmology
from simplemc.models.PhiCDMCosmology import PhiCDMCosmology

# Generic model
from simplemc.models.GenericCosmology import GenericCosmology
from simplemc.models.GenericPantheon import GenericPantheon

# Composite Likelihood
from simplemc.likelihoods.CompositeLikelihood import CompositeLikelihood

# Likelihood Multiplier
from simplemc.likelihoods.LikelihoodMultiplier import LikelihoodMultiplier

# Likelihood modules
from simplemc.likelihoods.BAOLikelihoods import DR11LOWZ, DR11CMASS, DR14LyaAuto, DR14LyaCross, \
    SixdFGS, SDSSMGS, DR11LyaAuto, DR11LyaCross, eBOSS, DR12Consensus
from simplemc.likelihoods.SimpleCMB import PlanckLikelihood, PlanckLikelihood_15, WMAP9Likelihood
from simplemc.likelihoods.CompressedSNLikelihood import BetouleSN, UnionSN, BinnedPantheon
from simplemc.likelihoods.PantheonSNLikelihood import PantheonSNLikelihood
from simplemc.likelihoods.HubbleParameterLikelihood import RiessH0
from simplemc.likelihoods.CompressedHDLikelihood import HubbleDiagram

from simplemc.likelihoods.GenericLikelihood import StraightLine
from simplemc.likelihoods.GenericPantheonSNLikelihood import GenericPantheonSNLikelihood
