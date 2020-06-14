__version__ = '1.9.0'
__author__  = 'JA Vazquez, I Gomez-Vargas, A Slosar'

import sys
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

if sys.version_info < (3, 0):
    logger.critical("You must use Python 3!")
    exit(1)

from .PostProcessing import PostProcessing
from .CosmoCalc import CosmoCalc
from .runbase import ParseDataset, ParseModel
from . import cosmo
from . import analyzers
from . import likelihoods
from . import models
from . import tools
from . import attick


__version__ = "1.9.0"