#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Fit Balmer lines. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from scipy import (interpolate, ndimage)

import utils

__all__ = ["BalmerLineModel"]


class BalmerLineModel(object):

    _minimum_required_spectra = 10

    def __init__(self, paths):
        """
        Initiate a Balmer line model to fit spectra.

        :param paths:
            A list of disk paths containing pre-computed spectra for a Balmer
            line profile.
        """

        if self._minimum_required_spectra > len(paths):
            raise ValueError("insufficient number of spectra provided")

        # The grid has different spectrum points at each set of stellar
        # parameters. But the synthesized region is the same: Extra points are
        # sampled to account for fast-changing flux. That's good for physics
        # but bad for interpolation. Instead we'll use the most sampled spectrum
        # as a common dispersion map.

        spectra = [utils.parse_spectrum(path) for path in paths]

        index = np.argmax([spectrum.size for spectrum in spectra])
        self.wavelengths = spectra[index][:, 0]

        # Skip over empty spectra and resample the others.
        
        fluxes = []
        usable_paths = []
        for i, (path, spectrum) in enumerate(zip(paths, spectra)):
            if spectrum.size > 0:
                usable_paths.append(path)

                # Interpolate the spectrum onto the common dispersion map.
                fluxes.append(np.interp(self.wavelengths, spectrum[:, 0], 
                    spectrum[:, 1], left=None, right=None))

        self.fluxes = np.array(fluxes)
        self.paths = usable_paths

        # Parse the stellar parameters from the usable paths.
        self.stellar_parameters = np.array(
            [utils.parse_stellar_parameters(path) for path in self.paths])

        # Create an interpolator.
        self._interpolator = interpolate.LinearNDInterpolator(
            self.stellar_parameters, self.fluxes, rescale=True)

        self.metadata = {
            "mask": [],
            "redshift": False,
            "smoothing": False,
            "continuum_order": -1
        }
        self._update_parameters()

        return None


    @property
    def parameter_names(self):
        """ Return the model parameters. """
        return self._parameter_names


    @property
    def parameter_bounds(self):
        return self._parameter_bounds


    def _update_parameters(self):
        """ Update the model parameter names based on the current metadata. """

        parameters = ["effective_temperature", "surface_gravity", "metallicity"]

        # Radial velocity?
        if self.metadata["redshift"]:
            parameters.append("redshift")

        # Smoothing?
        if self.metadata["smoothing"]:
            parameters.append("smoothing")

        # Continuum?
        parameters += ["c{}".format(i) \
            for i in range(self.metadata["continuum_order"] + 1)]

        self._parameter_names = tuple(parameters)
        self._parameter_bounds = {}
        return True


    def _initial_guess(self, **kwargs):
        """ Return an initial guess about the model parameters. """

        p0 = list(np.mean(self.stellar_parameters, axis=1))

        if self.metadata["redshift"]:
            p0.append(0)

        if self.metadata["smoothing"]:
            p0.append(5)

        if self.metadata["continuum_order"] > -1:
            p0.append(1)
            p0.extend([0] * self.metadata["continuum_order"])

        return np.array(p0)


    def mask(self, spectrum):
        """
        Return a pixel mask based on the metadata and existing mask information
        available.

        :param spectrum:
            A spectrum to generate a mask for.
        """

        mask = (spectrum.dispersion >= self.wavelengths[0]) \
             * (spectrum.dispersion <= self.wavelengths[-1])

        # Any masked ranges specified in the metadata?
        for start, end in self.metadata.get("mask", []):
            mask *= ~((spectrum.dispersion >= start) \
                    * (spectrum.dispersion <= end))
        return mask


    def interpolate(self, stellar_parameters, fill_value=1.0):
        """
        Interpolate a Balmer line profile at the given stellar parameters.

        :param stellar_parameters:
            A list-like object containing the stellar parameters.

        :returns:
            An array of fluxes the same length as `self.wavelength`.
        """

        flux = self._interpolator(*stellar_parameters)
        if fill_value is not None:
            flux[~np.isfinite(flux)] = fill_value
        return flux


    def __call__(self, x, parameters):
        """
        Generate data at x given the model and specified parameters.

        :param x:
            An array of dispersion points.

        :param parameters:
            A list of model parameter values.
        """

        normalized_flux = self.interpolate(*parameters[:3])

        # Continuum.
        O = self.metadata["continuum_order"]
        continuum = np.ones_like(x) \
            if 0 > O else np.polyval(parameters[-(O + 1):], x)    

        y = continuum * normalized_flux

        # Smoothing.
        try:
            index = self.parameter_names.index("smoothing")

        except ValueError:
            None

        else:
            y = ndimage.gaussian_filter(y, parameters[index])

        # Redshift.
        try:
            index = self.parameter_names.index("redshift")

        except ValueError:
            None

        else:
            z = parameters[index]
            y = np.interp(x, self.wavelengths * (1 + z), y, left=1, right=1)

        return y
        

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



    def fit(self, spectrum):

        # Get initial point.

        # Fit with mask.

        # Return optimized, cov, + metadata


        raise NotImplementedError



if __name__ == "__main__":



    import scipy.interpolate
    from glob import glob

    base_model = BalmerLineModel(glob("data/gam_*.prf"))

    raise a
    # For each point, remove it and create a new thing.
    sampled_paths = []
    sampled_fluxes = []
    interpolated_fluxes = []

    shuffle(base_model.paths)
    diffs = []

    N = len(base_model.paths)
    for i, path in enumerate(base_model.paths):
        
        test_paths = [] + base_model.paths
        removed_point = test_paths.pop(i)

        # Is this a valid point?
        synthesised_spectrum = utils.parse_spectrum(removed_point)
        if synthesised_spectrum.size == 0: continue

        sampled_paths.append(path)

        # Create model
        test_model = BalmerLineModel(test_paths)

        f = scipy.interpolate.LinearNDInterpolator(
            test_model.stellar_parameters, test_model.fluxes, rescale=True)

        # Interpolate to our test point.
        stellar_parameters = utils.parse_stellar_parameters(removed_point)

        # Put our interpolated flux onto the same wavelengths as the point removed.
        resampled_flux = np.interp(synthesised_spectrum[:,0],
            test_model.wavelength, f(*stellar_parameters))

        diff = synthesised_spectrum[:, 1] - resampled_flux
        diffs.append(diff)

        print("{0}/{1}: {2} (mean {3:.1e}, std {4:.1e} max abs {5:.1e})".format(
            i, N, path, np.nanmean(diff), np.nanstd(diff), np.nanmax(np.abs(diff))))

        #interpolated_fluxes.append(
        #    np.array([synthesised_spectrum[:,0], resampled_flux]).T)
        #sampled_fluxes.append(synthesised_spectrum)

    # Calculate deltas.




