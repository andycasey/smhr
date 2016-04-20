
#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model for fitting synthetic spectra to spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["SpectralSynthesisModel"]

import logging
import numpy as np
from six import string_types

from .base import BaseSpectralModel
from smh import utils
from smh.photospheres.abundances import asplund_2009 as solar_composition
# HACK TODO REVISIT SOLAR SCALE: Read from session defaults?


logger = logging.getLogger(__name__)


class SpectralSynthesisModel(BaseSpectralModel):
    """
    A class to fit spectral data by synthesis, given a photospheric structure.
    """
    
    def __init__(self, transitions, session, elements, *args, **kwargs):
        super(SpectralSynthesisModel, self).__init__(
            transitions, session, **kwargs)


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

        :param elements:
            The element(s) (string or list-type of strings) that will be
            measured in this model.
        """

        # Formatting checks.
        super(SpectralSynthesisModel, self)._verify_transitions()

        # Format the elements and then check that all are real.
        if isinstance(elements, string_types):
            elements = [elements]

        elements = [str(element).title() for element in elements]
        for element in elements:
            # Is the element real?
            if element not in utils.periodic_table:
                raise ValueError("element '{}' is not valid".format(element))

            # Is it in the transition list?
            if  element not in self.transitions["elem1"] \
            and element not in self.transitions["elem2"]:
                raise ValueError(
                    "element '{0}' does not appear in the transition list"\
                        .format(element))

        self.metadata["elements"] = elements

        return True


    def _initial_guess(self):
        """
        Return an initial (uninformed) guess about the model parameters.
        """

        # Potential parameters:
        # elemental abundances, continuum coefficients, smoothing kernel,
        # velocity offset

        # The elemental abundances are in log_epsilon format.
        # We will assume a scaled-solar initial value based on the stellar [M/H]

        p0 = []
        for parameter in self.parameter_names:

            if parameter in ("sigma_smooth", "vrad"):
                p0.append(0)

            elif parameter.startswith("log_eps"):
                # Elemental abundance.
                element = parameter.split("(")[1].rstrip(")")

                # Assume scaled-solar composition.
                p0.append(solar_composition(element) + \
                self.session.metadata["stellar_parameters"]["metallicity"])

            elif parameter.startswith("c"):
                # Continuum coefficient.
                p0.append((0, 1)[parameter == "c0"])

            else:
                raise ParallelUniverse("this should never happen")

        return np.array(p0)


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

        # Check the observed spectrum for validity.
        spectrum = spectrum or self.session.normalized_spectrum
        self._verify_spectrum(spectrum)

        # Update internal metadata with any input parameters.
        # Ignore additional parameters because other BaseSpectralModels will
        # have different input arguments.
        for key in set(self.metadata).intersection(kwargs):
            self.metadata[key] = kwargs[key]

        # Update the parameter names in case they have been updated due to the
        # input kwargs.
        self._update_parameter_names()
        
        # Get a bad initial guess.
        p0 = self._initial_guess(spectrum, **kwargs)

        # Build a mask based on the window fitting range, and prepare the data.
        mask = self.mask(spectrum)
        x, y = spectrum.dispersion[mask], spectrum.flux[mask]
        yerr, absolute_sigma = ((1.0/spectrum.ivar[mask])**0.5, True)
        if not np.all(np.isfinite(yerr)):
            yerr, absolute_sigma = (np.ones_like(x), False)

        # Fit the data!
        try:
            p_opt, p_cov = op.curve_fit(self.fitting_function, xdata=x, ydata=y,
                sigma=yerr, p0=p0, absolute_sigma=absolute_sigma)

        except:
            logger.exception("Exception raised in fitting {}".format(self))
            raise


        raise NotImplementedError



    def __call__(self, dispersion, *parameters):
        """
        Generate data at the dispersion points, given the parameters.

        :param dispersion:
            An array of dispersion points to calculate the data for.

        :param parameters:
            Keyword arguments of the model parameters and their values.
        """


        # Parse the elemental abundances, because they need to be passed to
        # the synthesis routine.
        abundances = {}
        names =  self.parameter_names
        for name, parameter in zip(names, parameters):
            if name.startswith("log_eps"):
                element = name.split("(")[1].rstrip(")")
                abundances[element] = parameter
            else: break # The log_eps abundances are always first.

        # Produce a synthetic spectrum.
        synth_dispersion, intensities = self.session.rt.synthesize(
            self.session.stellar_photosphere, self.transitions,
            abundances=abundances) # TODO: Other RT kwargs......

        # Continuum.
        O = self.metadata["continuum_order"]
        if 0 > O:
            continuum = 1
        else:
            continuum = np.polyval(
                [parameters[names.index("c{}".format(i))] for i in range(O)],
                synth_dispersion)

        model = intensities * continuum

        # Smoothing.
        try:
            index = names.index("sigma_smooth")
        except IndexError:
            None
        else:
            model = gaussian_filter(model, abs(parameters[index]))

        # Interpolate the model spectrum onto the requested dispersion points.
        return np.interp(dispersion, synth_dispersion, model, left=1, right=1)



    def abundances(self):
        """
        Calculate the abundances (model parameters) given the current stellar
        parameters in the parent session.
        """

        self.fit(self.session.normalized_spectrum)

        # Parse the abundances into the right format.

        

# Fit synthesised region to spectrum.

# HOOK in to any RT.


