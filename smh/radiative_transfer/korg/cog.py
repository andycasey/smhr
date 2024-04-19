#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Functionality to derive abundances from a transition's curve-of-growth using the
Korg radiative transfer code.
This doesn't actually work yet, because Korg assumes everything is on the linear part of the COG.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from juliacall import Main as jl
jl.seval("using Korg")
Korg = jl.Korg

import logging
import numpy as np
import re
import yaml
from pkg_resources import resource_stream

from . import utils
from .utils import RTError
from smh.utils import element_to_species
from smh import linelists

logger = logging.getLogger(__name__)


def abundance_cog(photosphere, transitions, full_output=False, verbose=False,
    twd=None, **kwargs):
    """
    Calculate atomic line abundances by interpolating the measured 
    equivalent width from the curve-of-growth. 
    This doesn't actually work yet, because Korg assumes everything is on the linear part of the COG.
    
    # :param photosphere:
    #     A formatted photosphere.

    # :param transitions:
    #     A list of atomic transitions with measured equivalent widths.

    # :param verbose: [optional]
    #     Specify verbose flags to MOOG. This is primarily used for debugging.

    """
    
    # Create a temporary directory.
    path = utils.twd_path(twd=twd,**kwargs)

    # Write out the photosphere.
    model_in, lines_in = path("model.in"), path("lines.in")
    ## TODO
    photosphere.write(model_in, format="marcs")
    ## TODO
    transitions.write(lines_in, format="moog")
    
    atm = Korg.read_photosphere(model_in)
    MH = TODO
    linelist = Korg.read_linelist(lines_in, format="moog")
    # The first argument specifies the default [X/H] abundance, the second lets you specify the default [alpha/H] abundance, and the third lets you specify individual abundances.
    A_X = Korg.format_A_X(MH) # the need for this should go away
    measured_EWs = transition["equivalent_width"]
    logeps = Korg.Fit.ews_to_abundances(atm, linelist, A_X, measured_EWs)
    # return logeps
    raise NotImplementedError
