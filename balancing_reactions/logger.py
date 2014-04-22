#!/usr/bin/env python

import logging
import time


#Set up logging
logfn = "rxn_retrieve_and_balance_{0}.log".format(time.strftime('%Y%m%d_%H%M'))
logger = logging.getLogger(logfn)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(logfn)
formatter = logging.Formatter('%(name)s: %(levelname)s: %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)
logger.info('Script is running, logging events to {0}'.format(logfn))
