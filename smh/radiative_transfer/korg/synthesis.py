#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Synthesis functionality with the Korg radiative transfer code. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)


import logging
import numpy as np
import yaml
from pkg_resources import resource_stream

from . import utils

logger = logging.getLogger(__name__)

def synthesize(photosphere, transitions, abundances=None, isotopes=None,
    verbose=False, twd=None, **kwargs):
    """
    Sythesize a stellar spectrum given the model photosphere and list of
    transitions provided with Korg.
    
    :param photosphere:
        A formatted photosphere.

    :param transitions:
        A list of atomic transitions with measured equivalent widths.

    :param abundances: [optional]
        The non-scaled-Solar abundances to use in the synthesis.

    :param isotopes: [optional]
        The isotopic fractions to use in the synthesis.

    :param verbose: [optional]
        Specify verbose flags to Korg. This is primarily used for debugging.

    """

    # Create a temporary directory and write out the photoshere and transitions.
    path = utils.twd_path(twd=twd,**kwargs)
    model_in, lines_in = path("model.in"), path("lines.in")
    ## TODO: can't currently write a correct MARCS format
    photosphere.write(model_in, format="marcs")
    transitions.write(lines_in, format="moog")
    wavemin, wavemax = transitions[0]-1.0,transitions[-1]+1.0
    atm = Korg.read_model_atmosphere(model_in)
    ## TODO: this doesn't support broadening parameters or dissociation energies
    ## TODO: Korg will automatically apply isotopic ratios when it reads, need to update this?
    linelist = Korg.read_linelist(lines_in, format="moog")
    # read_linelist(filename; format="vald", isotopic_abundances=Korg.isotopic_abundances)

    
    # Abundances.
    #A_X = Korg.format_A_X(-0.5, Dict("C" => -0.25))
    #solution = synthesize(atm, linelist, A_X, 5000, 5100)
    
    #synthesize(atm, linelist, A_X, λ_start, λ_stop, [λ_step=0.01]; kwargs... )
    #synthesize(atm, linelist, A_X, wavelength_ranges; kwargs... )
    
    #sol = Korg.synthesize(...)
    # sol.flux
    # sol.wavelengths
    # sol.cntm
    #flux: the output spectrum
    #cntm: the continuum at each wavelength
    #alpha: the linear absorption coefficient at each wavelength and atmospheric layer a Matrix of size (layers x wavelengths)
    #number_densities: A dictionary mapping Species to vectors of number densities at each atmospheric layer
    #electron_number_density: the electron number density at each atmospheric layer
    #wavelengths: The vacuum wavelengths (in Å) over which the synthesis was performed. If air_wavelengths=true this will not be the same as the input wavelengths.
    #subspectra: A vector of ranges which can be used to index into flux to extract the spectrum for each range provided in wavelength_ranges. If you use the standard λ_start, λ_stop, λ_step arguments, this will be a vector containing only one range.

    spectra = []
    all_abundances = []
    for abundances in all_abundances:
        out = Korg.synthesize(atm, linelist, abundances, wavemin, wavemax)
        meta = {}
        spectra.append((out.wavelengths, out.flux/out.cntm, meta))
    
    #return spectra
    raise NotImplementedError
