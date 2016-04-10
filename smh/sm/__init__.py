#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Models for fitting spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)


from .base import *
from .profile import *
from .synthesis import *



# Inputs common to EW and synthesis:
# - transitions
# - parent session

# Properties common to EW and synthesis:
# - transitions
# - parent session
# - metadata
# - fitting_method (be as part of metadata?)
# - parameters (named, and values if they are optimized/fixed)
# - bounds (bounds on parameters)
# - fixed (whether parameters are fixed or not)

# Methods common to EW and syntheses:
# - __call__ function to produce spectra given the model parameters
# - abundance property (internally store an abundance and the params it was calculated on? then only re-calculate from session model if those params
# have been updated??) --> probably not because that will TRIGGER, but do something clever


# EW-only:
# fit(spectrum, wl_range, *fitting_args)
# --> return optimised fitting parameters


# Synthesis-only:
# fit(spectrum, wl_range, teff, logg, which_atmospheres, etc, *fitting_args)
# --> returns optimized fitting parameters including abundances, etc
