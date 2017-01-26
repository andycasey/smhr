#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Models for fitting spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["BaseSpectralModel"]

import numpy as np

from .quality_constraints import constraints

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

        # Here we will have to just assume that the user knows what they are
        # doing, otherwise to import BaseSession implies that we will (probably)
        # have a recursive import loop.
        #if not isinstance(session, BaseSession):
        #    raise TypeError("session must be a sub-class from BaseSession")

        if len(session.metadata.get("line_list", [])) == 0:
            raise ValueError("session does not have a line list")

        # Validate the transition_hashes
        transition_hashes = list(transition_hashes)
        for transition_hash in transition_hashes:
            if transition_hash not in session.metadata["line_list"]["hash"]:
                raise ValueError(
                    "transition hash '{}' not found in parent session"\
                    .format(transition_hash))

        self._session = session
        self._transition_hashes = transition_hashes

        # Link the .transitions attribute to the parent session.
        self.index_transitions()
        
        self.metadata = {
            "use_for_stellar_composition_inference": True,
            "use_for_stellar_parameter_inference": (
                "Fe I" in self.transitions["element"] or
                "Fe II" in self.transitions["element"])
        }

        # Create a _repr_wavelength property.
        if len(self.transitions) == 1:
            self._repr_wavelength \
                = "{0:.1f}".format(self.transitions["wavelength"][0])
        else:
            self._repr_wavelength \
                = "~{0:.0f}".format(np.mean(self.transitions["wavelength"]))
        
        return None


    @property
    def wavelength(self):
        """
        Return a (sometimes approximate) wavelength for where this spectral line
        occurs.
        """

        wavelength = np.mean(self.transitions["wavelength"])
        return int(wavelength) if len(self.transitions) > 1 else wavelength


    @property
    def is_acceptable(self):
        """ Return whether this spectral model is acceptable. """
        return self.metadata.get("is_acceptable", False)


    @is_acceptable.setter
    def is_acceptable(self, decision):
        """
        Mark the spectral model as acceptable or unacceptable.

        :param decision:
            A boolean flag.
        """
        decision = bool(decision)
        if not decision or (decision and "fitted_result" in self.metadata):
            self.metadata["is_acceptable"] = bool(decision)
        return None


    @property
    def use_for_stellar_parameter_inference(self):
        """
        Return whether this spectral model should be used during the
        determination of stellar parameters.
        """
        return self.metadata["use_for_stellar_parameter_inference"]


    @use_for_stellar_parameter_inference.setter
    def use_for_stellar_parameter_inference(self, decision):
        """
        Mark whether this spectral model should be used when inferring stellar
        parameters.

        :param decision:
            A boolean flag.
        """
        self.metadata["use_for_stellar_parameter_inference"] = bool(decision)
        return None


    @property
    def use_for_stellar_composition_inference(self):
        """
        Return whether this spectral model should be used for the determination
        of stellar composition.
        """
        return self.metadata["use_for_stellar_composition_inference"]



    @use_for_stellar_composition_inference.setter
    def use_for_stellar_composition_inference(self, decision):
        """
        Mark whether this spectral model should be used when inferring the 
        stellar composition.

        :param decision:
            A boolean flag.
        """
        self.metadata["use_for_stellar_composition_inference"] = bool(decision)
        return None


    def apply_quality_constraints(self, quality_constraints):
        """
        Apply quality constraints to the this spectral model. If the model does
        not meet the quality constraints it will be marked as unacceptable.

        :param quality_constraints:
            A dictionary containing constraint names as keys and a 2-length
            tuple with the (lower, upper) bounds specified as values.

        :returns:
            Whether this model met the specified quality constraints.
        """

        is_ok = constraints(self, quality_constraints)
        self.is_acceptable = is_ok
        return is_ok


    def meets_quality_constraints(self, quality_constraints):
        """
        Check whether this spectral model meets specified quality constraints.

        :param quality_constraints:
            A dictionary containing constraint names as keys and a 2-length
            tuple with the (lower, upper) bounds specified as values.

        :returns:
            Whether this model met the specified quality constraints.
        """
        return constraints(self, quality_constraints)


    @property
    def meets_quality_constraints_in_parent_session(self):
        """
        Returns whether this spectral model meets the quality constraints set
        in the parent session.
        """

        return self.meets_quality_constraints(
            self._session.setting(("spectral_model_quality_constraints", )) or {})


    @property
    def transitions(self):
        """ Return the transitions associateed with this class. """

        # This is left as a property to prevent users from arbitrarily setting
        # the .transitions attribute.
        return self._transitions


    def index_transitions(self):
        """
        Index the transitions to the parent session.
        
        This step is very slow for large linelists.
        """

        if len(self._transition_hashes) < 50:
        #if True:
            ## Brute force loop for small N
            indices = np.zeros(len(self._transition_hashes), dtype=int)
            for i, each in enumerate(self._transition_hashes):
                index = np.where(
                    self._session.metadata["line_list"]["hash"] == each)[0]
                #assert len(index) == 1, len(index)
                if len(index) != 1: 
                    #print("WARNING: hash {} appears {} times in session linelist".format(each,len(index)))
                    index = index[0]
                indices[i] = index
        else:
            ## Use searchsorted binary search, speeds up by ~1.8x
            #iisort = np.argsort(self._session.metadata["line_list"]["hash"])
            iisort = self._session.metadata["line_list_argsort_hashes"]
            sorted = np.searchsorted(self._session.metadata["line_list"]["hash"],
                                     self._transition_hashes,
                                     sorter=iisort)
            indices = iisort[sorted]
            ## Check for correctness, i.e. all hashes are actually in linelist
            assert np.all(self._transition_hashes == self._session.metadata["line_list"]["hash"][indices])

        self._transition_indices = indices
        self._transitions = self._session.metadata["line_list"][indices]

        return indices


    @property
    def elements(self):
        """ Return the elements to be measured from this class. """
        return self.metadata["elements"]


    @property
    def species(self):
        """ Return the species to be measured from this class. """
        return self.metadata["species"]


    @property
    def session(self):
        """ Return the parent session that this model is associated with. """
        return self._session

    @property
    def abundances(self):
        """ Return abundances if fit, else None """
        try:
            return self.metadata["fitted_result"][2]["abundances"]
        except KeyError:
            return None
        

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


    def __getstate__(self):
        """ Return a serializable state of this spectral model. """

        state = {
            "type": self.__class__.__name__,
            "transition_hashes": self._transition_hashes,
            "metadata": self.metadata
        }
        return state


    def __setstate__(self, state):
        """ Disallow the state to be instantiated from a serialised object. """
        return None


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

        # HACK
        if "antimask_flag" not in self.metadata:
            self.metadata["antimask_flag"] = False
        if self.metadata["antimask_flag"]:
            antimask = np.ones_like(spectrum.dispersion,dtype=bool)
            for start, end in self.metadata["mask"]:
                antimask *= ~((spectrum.dispersion >= start) \
                            * (spectrum.dispersion <= end))
            return ~antimask

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

