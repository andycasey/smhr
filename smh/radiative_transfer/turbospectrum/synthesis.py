#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Synthesis functionality using the Turbospectrum radiative transfer code. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import numpy as np
import yaml
from pkg_resources import resource_stream

logger = logging.getLogger(__name__)

# Load defaults.
with resource_stream(__name__, "defaults.yaml") as fp:
    _turbospectrum_defaults = yaml.load(fp)



def synthesize(photosphere, transitions, abundances=None, isotopes=None,
    verbose=False, twd=None, **kwargs):
    """
    Sythesize a stellar spectrum given the model photosphere and transitions
    provided. 

    """

    raise a
    # Get regions to synthesise.


    # Update abundances

    # Calculate opacities

    # Synthesise



