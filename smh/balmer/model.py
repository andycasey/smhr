#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Fit Balmer lines. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from scipy import interpolate

import utils

__all__ = ["BalmerLineModel"]


class BalmerLineModel(object):

    _minimum_required_spectra = 10

    def __init__(self, paths):
        """
        Initiate a Balmer line model to fit spectra.

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
        self.wavelength = spectra[index][:, 0]

        # Skip over empty spectra and resample the others.
        
        fluxes = []
        usable_paths = []
        for i, (path, spectrum) in enumerate(zip(paths, spectra)):
            if spectrum.size > 0:
                usable_paths.append(path)

                # Interpolate the spectrum onto the common dispersion map.
                fluxes.append(np.interp(self.wavelength, spectrum[:, 0], 
                    spectrum[:, 1], left=None, right=None))

        self.fluxes = np.array(fluxes)
        self.paths = usable_paths

        # Parse the stellar parameters from the usable paths.
        self.stellar_parameters = np.array(
            [utils.parse_stellar_parameters(path) for path in self.paths])

        # Create an interpolator.
        self._interpolator = interpolate.griddata(
            self.stellar_parameters, self.fluxes, rescale=True)

        return None



    def interpolate(self, stellar_parameters, fill_value=1.0):
        """
        Interpolate a Balmer line profile at the given stellar parameters.

        :param stellar_parameters:
            A list-like object containing the stellar parameters.

        :returns:
            An array of fluxes the same length as `self.wavelength`.
        """

        print(stellar_parameters)

        flux = self._interpolator(*stellar_parameters)
        if fill_value is not None:
            flux[~np.isfinite(flux)] = fill_value
        return flux




if __name__ == "__main__":



    import scipy.interpolate
    from glob import glob

    base_model = BalmerLineModel(glob("data/gam_*.prf"))

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




