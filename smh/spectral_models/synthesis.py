
#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model for fitting synthetic spectra to spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["SpectralSynthesisModel"]

import logging
from six import string_types

from .base import BaseSpectralModel

logger = logging.getLogger(__name__)


class SpectralSynthesisModel(BaseSpectralModel):
    """
    A class to fit spectral data by synthesis, given a photospheric structure.
    """
    
    def __init__(self, transitions, session, elements, *args, **kwargs):
        super(SpectralSynthesisModel, self).__init__(
            transitions, sessions, **kwargs)


        # Initialize metadata with default fitting values.
        self.metadata.update({
            "mask": [],
            "window": 2, 
            "continuum_order": 1,
            "velocity_tolerance": None,
            "smoothing_kernel": True
        })

        self._verify_transitions(elements)

        # Inputs about which RT code to use, stellar params + RT inputs, etc.


        # Set the model parameter names based on the current metadata.
        self._update_parameter_names()


        return None


    def _verify_transitions(self, elements):
        """
        Verify that the atomic or molecular transitions associated with this
        class are valid.
        """

        # Check format first.
        super(ProfileFittingModel, self)._verify_transitions()
        if isinstance(elements, string_types):
            elements = [elements]

        self.metadata["elements"] = elements

        # Check each element is real and that it exists in the transition table.
        # TODO:

        return True


    def _update_parameter_names(self):
        """ Update the model parameter names based on the current metadata. """

        bounds = {}

        parameter_names = []

        # Abundances of different elements to fit simultaneously.
        parameter_names.extend(
            ["log_eps({0})".format(e) for e in self.metadata["elements"]])

        # Radial velocity?
        vt = abs(self.metadata["velocity_tolerance"] or 0)
        if vt > 0:
            parameter_names.append("vrad")
            bounds["vrad"] = [-vt, +vt]

        # Continuum coefficients?
        parameter_names += ["c{0}".format(i) \
            for i in range(self.metadata["continuum_order"] + 1)]

        # Gaussian smoothing kernel?
        if self.metadata["smoothing_kernel"]:
            # TODO: Better init of this
            bounds["sigma_smooth"] = (0, 1)
            parameter_names.append("sigma_smooth")

        self._parameter_bounds = bounds
        self._parameter_names = parameter_names
        return True


    def fit(self, spectrum=None, **kwargs):
        """
        Fit a synthesised model spectrum to the observed spectrum.

        :param spectrum: [optional]
            The observed spectrum to fit the synthesis spectral model. If None
            is given, this will default to the normalized rest-frame spectrum in
            the parent session.
        """

        spectrum = spectrum or self.session.normalized_spectrum

        

        raise NotImplementedError



    def __call__(self, dispersion, *parameters):
        """
        Generate data at the dispersion points, given the parameters.

        :param dispersion:
            An array of dispersion points to calculate the data for.

        :param parameters:
            Keyword arguments of the model parameters and their values.
        """



        # Parse parameters.

        # Call synthesis via radiative transfer by accessing parent session
        # stellar atmosphere.

        # apply convolution, velocity + continuum parameters.

        # Interpolate onto dispersion points.

        raise NotImplementedError



    def abundances(self):
        """
        Calculate the abundances (model parameters) given the current stellar
        parameters in the parent session.
        """

        self.fit(self.session.normalized_spectrum)

        # Parse the abundances into the right format.

        

# Fit synthesised region to spectrum.

# HOOK in to any RT.


