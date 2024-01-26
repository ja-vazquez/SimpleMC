
import sys
import logging
import datetime
from pathlib import Path

cdir = str(Path(__file__).parent.resolve())

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