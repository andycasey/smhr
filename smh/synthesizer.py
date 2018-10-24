#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
An object to quickly and easily generate synthetic spectra
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Synthesizer"]

import numpy as np
import matplotlib.pyplot as plt
import logging

from .linelists import LineList
from . import (photospheres, radiative_transfer, specutils, isoutils, utils)
from .spectral_models import SpectralSynthesisModel
from smh.photospheres.abundances import asplund_2009 as solar_composition

logger = logging.getLogger(__name__)

class Synthesizer(object):
    """
    An object with convenience functions for creating synthetic spectra
    """
    def __init__(self,twd=None,isotopes=None,R=None,verbose=False):
        from .radiative_transfer import moog as rt
        self.rt = rt
        self.interpolator = photospheres.interpolator()

        # Make a working directory
        if twd is None:
            twd = utils.mkdtemp()
        self.twd = twd
        if verbose:
            logger.info("Synthesizer working directory: {}".format(twd))
        
        if isotopes is None:
            # Ba is r-process isotopes
            self.isotopes = {"Eu":{151:0.48,153:0.52},
                             "Ba":{134:0.0,135:0.370,136:0.0,137:0.351,138:0.279},
                             "H-C":{112:0.8,113:0.2},
                             "C-C":{1212:0.8,1213:0.2},
                             "C-N":{1214:0.8,1314:0.2}}
        
        try: self.R = float(R)
        except: self.R = None
        
        self.verbose = verbose

    def run_synth(self, atmosphere, linelist, abundances, smooth_fwhm=None, isotopes=None, verbose=False,
                  **kwargs):
        """
        Calculate a single synthetic spectrum

        TODO allow easy changing of moog parameters (e.g. damping)
        """
        verbose = verbose or self.verbose
        if isotopes is None:
            isotopes = self.isotopes
            if verbose:
                logger.info("Using default isotopes")
                logger.info(isotopes)
        # First assume you have a set of stellar parameters for "atmosphere"
        if not isinstance(atmosphere, photospheres.photosphere.Photosphere):
            if len(atmosphere) == 4:
                Teff, logg, vt, MH = atmosphere
                alpha = 0.4
            elif len(atmosphere) == 5:
                Teff, logg, vt, MH, alpha = atmosphere
            else:
                raise ValueError(atmosphere)
            atmosphere = self.make_photosphere(Teff, logg, vt, MH, alpha)
            
        # Run synthesis, create a spectrum with huge S/N
        out = self.rt.synthesize(atmosphere, linelist,
                                 abundances=abundances,
                                 isotopes=isotopes, verbose=verbose,
                                 twd=self.twd, **kwargs)
        if len(out) > 1:
            raise ValueError("Currently don't allow multiple abundances")
        out = out[0]
        spec = specutils.Spectrum1D(out[0],out[1],1e10*np.ones(len(out[0])))
        
        # Smooth if requested
        if self.R is not None or smooth_fwhm is not None:
            if self.R is not None:
                wlmid = np.median(linelist["wavelength"])
                smooth_fwhm = wlmid/self.R
            if verbose:
                logger.info("Smoothing to FWHM={:.3f}A".format(smooth_fwhm))
            spec = spec.gaussian_smooth(smooth_fwhm)
        return spec

    def make_photosphere(self, Teff, logg, vt, MH, alpha):
        """ Make a stellar photosphere """
        photosphere = self.interpolator(Teff, logg, MH, alpha)
        photosphere.meta["stellar_parameters"]["microturbulence"] = vt
        return photosphere

    def plot(self):
        """ Do something for making easy plots """
        raise NotImplementedError

    def run_abfind(self, atmosphere, linelist, verbose=False,
                   **kwargs):
        return self.rt.abundance_cog(atmosphere, linelist, verbose=verbose, twd=self.twd,
                                     **kwargs)
    
