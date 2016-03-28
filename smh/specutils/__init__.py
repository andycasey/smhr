#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Spectroscopy-related utilities. """

import logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))

#handler.setFormatter(logging.Formatter(
#    "%(asctime)s [%(levelname)-8s] (%(name)s/%(lineno)d): %(message)s"))

logger.addHandler(handler)

#simplefilter("ignore", RankWarning)
#simplefilter("ignore", RuntimeWarning)

import motions
from spectrum import Spectrum1D
from rv import *
from utils import *
