__version__ = '1.9.0'
__author__  = 'JA Vazquez, I Gomez-Vargas, A Slosar'

import sys
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

if sys.version_info < (3, 0):
    logger.critical("You must use Python 3!")
    exit(1)