#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model for fitting an absorption profile to spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["ProfileFittingModel"]

import logging
import numpy as np
from scipy.special import wofz

from base import BaseSpectralModel

logger = logging.getLogger(__name__)



def gaussian(x, *parameters):
    position, sigma, amplitude = parameters
    return amplitude * np.exp(-(x - position)**2 / (2.0 * sigma**2))

def lorentzian(x, *parameters):
    position, width, amplitude = parameters
    return (amplitude/np.pi) * (width/((x - position)**2 + width**2))

def voigt(x, *parameters):
    
    try:
        n = len(x)
    except TypeError:
        n = 1

    position, fwhm, amplitude, shape = parameters
    
    profile = 1 / wofz(np.zeros((n)) + 1j * np.sqrt(np.log(2.0)) * shape).real
    profile *= amplitude * wofz(2*np.sqrt(np.log(2.0)) * (x - position)/fwhm \
        + 1j * np.sqrt(np.log(2.0))*shape).real
    return profile



class ProfileFittingModel(BaseSpectralModel):
    """
    A class to fit spectral data with an absorption profile and continuum.
    """

    _profiles = {
        "gaussian": (gaussian, ("mean", "sigma", "amplitude")),
        "lorentzian": (lorentzian, ("mean", "width", "amplitude")),
        "voigt": (voigt, ("mean", "fwhm", "amplitude", "shape"))
    }

    def __init__(self, *args, **kwargs):
        super(ProfileFittingModel, self).__init__(*args, **kwargs)

        # Initialize metadata with default fitting values.
        self.metadata.update({
            "profile": "gaussian",
            "local_continuum": False,
            "central_weighting": True,
            "window": 10,
            "continuum_order": 2,
            "detection_sigma": 1.0,
            "max_iterations": 4
        })

        # Set the model parameter names based on the current metadata.
        self._update_parameters()
        return None


    @property
    def parameter_names(self):
        """ Return the model parameter names. """
        return self._parameter_names


    def _update_parameters(self):
        """ Update the model parameter names based on the current metadata. """

        # Profile parameter names.
        func, parameter_names = self._profiles[self.metadata["profile"]]
        parameter_names = list(parameter_names)

        # Continuum coefficients.
        if self.metadata["local_continuum"]:
            parameter_names += ["c{0}".format(i) \
                for i in range(self.metadata["continuum_order"])]

        self._parameter_names = parameter_names
        return True


    def fit(self, spectrum, **kwargs):
        """
        Fit an asborption profile to the transition in the spectrum.

        :param spectrum:
            The observed spectrum to fit the profile transition model.
        """

        # Update internal metadata with any input parameters.
        # Ignore additional parameters because other BaseSpectralModels will
        # have different input arguments.
        for key in set(self.metadata).intersection(kwargs):
            self.metadata[key] = kwargs[key]


        # What model parameters are in the fitting process?
        # In synth these would be abundances, etc. Here they are profile/cont
        # parameters.


        # Optimize the model parameters.

        # Store and return the model parameters.




    def __call__(self, dispersion, *parameters):
        """
        Generate data at the dispersion points, given the parameters.

        :param dispersion:
            An array of dispersion points to calculate the data for.

        :param parameters:
            Keyword arguments of the model parameters and their values.
        """

        function, profile_parameters = self._profiles[self.metadata["profile"]]

        N = len(profile_parameters)
        y = 1.0 - function(dispersion, *parameters[:N])

        # Continuum.
        if parameters[N:]:
            y *= np.polyval(parameters[N:], dispersion)
        
        if len(parameters) != len(self.parameter_names):
            logger.warn("Number of parameters given "
                        "did not match the expected number ({0} != {1})".format(
                            len(parameters), len(self.parameter_names)))

        return y



if __name__ == "__main__":


    import smh
    a = smh.Session([
        "/Users/arc/codes/smh/hd44007red_multi.fits",
        "/Users/arc/codes/smh/hd44007blue_multi.fits",
    ])


    foo = ProfileFittingModel(None, a)


