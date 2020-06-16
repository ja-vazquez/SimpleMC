__version__ = '1.9.0'
__author__  = 'JA Vazquez, I Gomez-Vargas, A Slosar'

import sys
import logging
import datetime

log = False

if log:
    date = datetime.datetime.now()
    date = date.strftime("%Y-%b-%d, %A %I:%M:%S")
    logging.basicConfig(filename="simplemc_{}.log".format(date),
                                filemode='a',
                                format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                                datefmt='%H:%M:%S',
                                level=logging.INFO)
else:
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
from . import plots
