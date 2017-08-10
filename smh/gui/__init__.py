#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
from numpy import RankWarning
from warnings import simplefilter

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

# Do something with the GUI logger. Put to file? Show warnings? Save to sesion?

# Ignore warnings in the GUI, but not in the main code.
simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)
