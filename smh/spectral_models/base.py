#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Models for fitting spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["BaseSpectralModel"]

import numpy as np
from smh.session import BaseSession


class BaseSpectralModel(object):

    def __init__(self, session, transition_hashes, **kwargs):
        """
        Initialize a base class for modelling spectra.

        :param session:
            The session that this spectral model will be associated with.

        :param transition_hashes:
            The hashes of transitions in the parent session that will be
            associated with this model.
        """

        if not isinstance(session, BaseSession):
            raise TypeError("session must be a sub-class from BaseSession")

        if "line_list" not in session.metadata \
        or len(session.metadata["line_list"]) == 0:
            raise ValueError("session does not have a line list")

        # Validate the transition_hashes
        transition_hashes = list(transition_hashes)
        for transition_hash in transition_hashes:
            if transition_hash not in session.metadata["line_list"]["hash"]:
                raise ValueError(
                    "transition hash '{}' not found in parent session".format(
                        transition_hash))

        self._session = session
        self._transition_hashes = transition_hashes

        indices = self._index_transitions()
        transitions = self._session.metadata["line_list"][indices]

        self.metadata = {
            "use_for_stellar_composition_inference": True,
            "use_for_stellar_parameter_inference": (
                "Fe I" in transitions["element"] or
                "Fe II" in transitions["element"])
        }

        # Create a _repr_wavelength property.
        if len(transitions) == 1:
            self._repr_wavelength \
                = "{0:.1f}".format(transitions["wavelength"][0])
        else:
            self._repr_wavelength \
                = "~{0:.0f}".format(np.mean(transitions["wavelength"]))

        return None


    @property
    def is_acceptable(self):
        """ Return whether this spectral model is acceptable. """
        return self.metadata.get("is_acceptable", False)


    @property
    def use_for_stellar_parameter_inference(self):
        """
        Return whether this spectral model should be used during the
        determination of stellar parameters.
        """
        return self.metadata["use_for_stellar_parameter_inference"]


    @property
    def use_for_stellar_composition_inference(self):
        """
        Return whether this spectral model should be used for the determination
        of stellar composition.
        """
        return self.metadata["use_for_stellar_composition_inference"]


    @property
    def transitions(self):
        """ Return the transitions associateed with this class. """

        try:
            indices = self._transition_indices
        except AttributeError:
            indices = self._index_transitions()
        else:
            # Check hashes.
            actual_hashes = self._session.metadata["line_list"][indices]["hash"]
            if not np.all(actual_hashes == self._transition_hashes):
                indices = self._index_transitions()

        return self._session.metadata["line_list"][indices]

    def _index_transitions(self):
        """
        Index the transitions to the parent session.
        """

        indices = np.zeros(len(self._transition_hashes), dtype=int)
        for i, each in enumerate(self._transition_hashes):
            index = np.where(
                self._session.metadata["line_list"]["hash"] == each)[0]
            assert len(index) == 1
            indices[i] = index

        self._transition_indices = indices
        return indices



    @property
    def elements(self):
        """ Return the elements to be measured from this class. """
        return self.metadata["elements"]

    @property
    def session(self):
        """ Return the parent session that this model is associated with. """
        return self._session


    @property
    def parameters(self):
        """
        Return the model parameters.
        This must be implemented by the sub-classes.
        """
        raise NotImplementedError(
            "parameters property must be implemented by the sub-classes")


    @property
    def parameter_bounds(self):
        """ Return the fitting limits on the parameters. """
        return self._parameter_bounds


    @property
    def parameter_names(self):
        """ Return the model parameter names. """
        return self._parameter_names


    def __call__(self, dispersion, *args, **kwargs):
        """ The data-generating function. """
        raise NotImplementedError(
            "the data-generating function must be implemented by sub-classes")


    def _verify_transitions(self):
        """
        Verify that the transitions provided are valid.
        """
        # TODO
        return True


    def _verify_spectrum(self, spectrum):
        """
        Check that the spectrum provided is valid and has data in the wavelength
        range that we are interested.

        :param spectrum:
            The observed rest-frame normalized spectrum.
        """

        spectrum = spectrum or self.session.normalized_spectrum

        # Check the transition is in the spectrum range.
        wavelength = self.transitions["wavelength"]
        try:
            wavelength = wavelength[0]
        except IndexError:
            None
        if wavelength + 1 > spectrum.dispersion[-1] \
        or wavelength - 1 < spectrum.dispersion[0]:
            raise ValueError(
                "the observed spectrum contains no data over the wavelength "
                "range we require")

        return spectrum


    def mask(self, spectrum):
        """
        Return a pixel mask based on the metadata and existing mask information
        available.

        :param spectrum:
            A spectrum to generate a mask for.
        """

        window = abs(self.metadata["window"])
        wavelengths = self.transitions["wavelength"]
        try:
            lower_wavelength = wavelengths[0]
            upper_wavelength = wavelengths[-1]
        except IndexError:
            # Single row.
            lower_wavelength, upper_wavelength = (wavelengths, wavelengths)

        mask = (spectrum.dispersion >= lower_wavelength - window) \
             * (spectrum.dispersion <= upper_wavelength + window)

        # Any masked ranges specified in the metadata?
        for start, end in self.metadata["mask"]:
            mask *= ~((spectrum.dispersion >= start) \
                    * (spectrum.dispersion <= end))
        return mask


    def fitting_function(self, dispersion, *parameters):
        """
        Generate data at the dispersion points, given the parameters, but
        respect the boundaries specified on model parameters.

        :param dispersion:
            An array of dispersion points to calculate the data for.

        :param parameters:
            Keyword arguments of the model parameters and their values.
        """

        for parameter_name, (lower, upper) in self.parameter_bounds.items():
            value = parameters[self.parameter_names.index(parameter_name)]
            if not (upper >= value and value >= lower):
                return np.nan * np.ones_like(dispersion)

        return self.__call__(dispersion, *parameters)


    def _fill_masked_arrays(self, spectrum, x, *y):
        """
        Detect masked regions and fill masked regions in y-axis arrays with
        NaNs.

        :param spectrum:
            The spectrum used in the fit.

        :param x:
            The x values that were used in the fit.

        :param *y:
            The y-axis arrays to fill.
        """

        indices = spectrum.dispersion.searchsorted(x)
        x_actual = spectrum.dispersion[indices[0]:1 + indices[-1]]

        filled_arrays = [x_actual]
        for yi in y:
            yi_actual = np.nan * np.ones_like(x_actual)
            if len(yi_actual.shape) == 2:
                yi_actual[:, indices - indices[0]] = yi
            else:
                yi_actual[indices - indices[0]] = yi
            filled_arrays.append(yi_actual)

        return tuple(filled_arrays)

