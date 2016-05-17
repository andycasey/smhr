#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model photospheres """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
from pkg_resources import resource_stream


from .castelli_kurucz import Interpolator as ck_interp
from .marcs import Interpolator as marcs_interp
from .stagger import Interpolator as stagger_interp
from . import utils

logger = logging.getLogger(__name__)

def interpolator(kind="castelli/kurucz", **kwargs):

    logger.debug("Initialising {0} photosphere interpolator"\
        .format(kind.upper()))
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


# Detect which photospheres are available.
_photospheres = [
    ("MARCS 1D (2011)",
        "marcs",  "marcs-2011-standard.pkl"),
    ("Castelli & Kurucz 1D (2004)",
        "castelli/kurucz", "castelli-kurucz-2004.pkl"),
    ("Stagger <3D> (2013; optical)",
        "stagger-o", "stagger-2013-optical.pkl"),
    #("Stagger <3D> (2013; column mass density)",
    #    "stagger-m", "stagger-2013-mass-density.pkl"),
    #("Stagger <3D> (2013; height)",
    #    "stagger-h", "stagger-2013-height.pkl"),
]
available = []
for description, kind, basename in _photospheres:

    # Try to open the photosphere grid.
    try:
        with resource_stream(__name__, basename) as fp:
            None

    except ValueError:
        continue

    else:
        available.append((description, kind, basename))