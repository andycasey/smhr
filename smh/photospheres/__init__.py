#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model photospheres """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging

from .photosphere import Photosphere
from .abundances import asplund_2009 as solar_abundance
from .castelli_kurucz import Interpolator as ck_interp
from .marcs import Interpolator as marcs_interp
from .stagger import Interpolator as stagger_interp
from . import utils

logger = logging.getLogger("oracle")

def interpolator(kind="castelli/kurucz", **kwargs):

    logger.debug("Initialising {0} photosphere interpolator".format(kind.upper()))
    kind = kind.lower()
    if kind == "castelli/kurucz":
        return ck_interp(**kwargs)

    elif kind == "marcs":
        return marcs_interp(**kwargs)

    elif kind[:9] == "stagger-o":
        return stagger_interp("stagger-2013-optical.pkl", **kwargs)

    elif kind[:9] == "stagger-m":
        return stagger_interp("stagger-2013-mass-density.pkl", **kwargs)

    elif kind[:9] == "stagger-r":
        return stagger_interp("stagger-2013-rosseland.pkl", **kwargs)

    elif kind[:8] == "stagger-" and kind[8] in ("h", "g", "z"):
        return stagger_interp("stagger-2013-height.pkl", **kwargs)

    else:
        raise ValueError("'{}' model photospheres not recognised".format(kind))
    