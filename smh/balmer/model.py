#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Fit Balmer lines. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from scipy import (ndimage, optimize as op)

import utils

__all__ = ["BalmerLineModel"]



def mask_common_region(observed_spectrum, model_wavelengths, mask=None):
    """
    Mask the region in common between the observed and model spectra, and apply
    masks given.

    """

    common_region = [
        np.max([observed_spectrum.dispersion[0], model_wavelengths[0]]),
        np.min([observed_spectrum.dispersion[-1], model_wavelengths[-1]])
    ]



def _fit_single_model(observed_spectrum, model_wavelength, model_flux, mask=None,
    redshift=True, smoothing=True, continuum_order=-1):
    raise a

class BalmerLineModel(object):

    def __init__(self, paths, **kwargs):
        """
        A model for fitting Balmer lines.

        :param path:
            The disk location of a pre-computed Balmer line profile.
        """

        # Put all the spectra onto the most densely-sampled spectrum.
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

        self.normalized_fluxes = np.array(fluxes)
        self.paths = usable_paths

        # Parse the stellar parameters from the usable paths.
        self.stellar_parameters = np.array(
            [utils.parse_stellar_parameters(path) for path in self.paths])

        self.metadata = {
            "mask": [],
            "redshift": False,
            "smoothing": False,
            "continuum_order": -1
        }
        for key, value in kwargs.items():
            if key in self.metadata:
                self.metadata[key] = value

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

        parameters = []

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



    def __call__(self, x, *parameters):

        # Marginalize over all models.

        # Continuum.
        O = self.metadata["continuum_order"]
        continuum = 1.0 if 0 > O \
            else np.polyval(parameters[-(O + 1)::][::-1], self.wavelengths)

        y = np.atleast_2d(self.normalized_fluxes * continuum)

        # Smoothing?
        try:
            index = self.parameter_names.index("smoothing")
        except ValueError:
            None
        else:
            y = ndimage.gaussian_filter(y, parameters[index], axis=1)

        # Redshift?
        try:
            index = self.parameter_names.index("redshift")
        except ValueError:
            z = 0
        else:
            z = parameters[index]

        return np.array([np.interp(x, self.wavelengths * (1 + z), yi, 
            left=1, right=1) for yi in y])



    @property
    def initial_guess(self):
        """ Generate an uninformed guess for the model parameters. """

        defaults = {
            "c0": 1,
            "redshift": 0,
            "smoothing": 5
        }
        return tuple([defaults.get(p, 0) for p in self.parameter_names])


    def fit(self, observed_spectrum):
        """
        Fit this Balmer line model (and the nuisance parameters) to the data.

        :param observed_spectrum:
            An observed spectrum.
        """

        self._update_parameters()

        show = (observed_spectrum.dispersion >= self.wavelengths[0]) \
             * (observed_spectrum.dispersion <= self.wavelengths[-1])

        mask = show * _generate_mask(observed_spectrum.dispersion, self.metadata["mask"])

        observed_dispersion = observed_spectrum.dispersion
        observed_flux = observed_spectrum.flux
        observed_ivar = observed_spectrum.ivar
        observed_ivar[~mask] = 0


        if 1 > len(self.parameter_names):
            chi_sq = (self(observed_dispersion, None) - observed_flux)**2 \
                * observed_ivar

            return ((), None, chi_sq)


        observed_sigma = np.sqrt(1.0/observed_ivar)
        ok = np.isfinite(observed_flux * observed_sigma)


        def marginalize(x, *parameters):
            all_model_spectra = self(x, *parameters)

            # Get best.
            chi_sq = np.nansum(
                (all_model_spectra - observed_flux[ok])**2 * observed_ivar[ok],
                axis=1)
            index = np.argmin(chi_sq)
            print(index, parameters)
            return all_model_spectra[index]


        op_parameters, cov = op.curve_fit(marginalize, observed_dispersion[ok],
            observed_flux[ok], sigma=observed_sigma[ok],
            p0=self.initial_guess, absolute_sigma=True)

        fig, ax = plt.subplots()

        ax.plot(observed_dispersion[show], observed_flux[show], c='k')
        ax.plot(observed_dispersion[ok], marginalize(observed_dispersion[ok], *op_parameters), c='r')

        raise a
        



def _generate_mask(x, masked_regions=None):
    if masked_regions is None:
        masked_regions = []

    mask = np.ones(x.size, dtype=bool)
    for start, end in masked_regions:
        mask *= ~((x >= start) * (x <= end))
    return mask



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







if __name__ == "__main__":



    from glob import glob

    paths = glob("data/gam_*.prf")
    model = BalmerLineModel(paths, mask=[
        [1000, 4330],
        [4336.84, 4337.31],
        [4337.48, 4338.09],
        [4339.31, 4339.53],
        [4340.03, 4340.98],
        [4341.2, 4341.59],
        [4344.11, 4344.72],
        [4350, 10000]
    ],
    redshift=True,
    continuum_order=2
    )


    # Get a rest-frame normalized spectrum...
    from smh.specutils import Spectrum1D
    hd122563 = Spectrum1D.read("hd122563.fits")


