#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Fit Balmer lines. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from scipy import (integrate, ndimage, optimize as op)

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
            "continuum_order": -1,
            "marginalization_boundaries": {
                "smoothing": [0, 10],
                "redshift": [-30/299792.458, +30,299792.458],
            }
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



    def __call__(self, model_index, x, *parameters):

        # Marginalize over all models.

        # Continuum.
        O = self.metadata["continuum_order"]
        continuum = 1.0 if 0 > O \
            else np.polyval(parameters[-(O + 1)::][::-1], self.wavelengths)

        if model_index is not None:
            y = np.atleast_2d(self.normalized_fluxes[model_index, :] * continuum)
        else:
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




    def _slice_spectrum(self, observed_spectrum):

        idx = observed_spectrum.dispersion.searchsorted([
            self.wavelengths[0], self.wavelengths[-1]])

        # Account for 'right' hand behaviour.
        idx[-1] = np.clip(idx[1] + 1, 0, observed_spectrum.dispersion.size)

        return (
            observed_spectrum.dispersion.__getslice__(*idx),
            observed_spectrum.flux.__getslice__(*idx),
            observed_spectrum.ivar.__getslice__(*idx)
        )


    def _optimize_nuisance_parameters(self, dispersion, flux, sigma, grid_index,
        **kwargs):
        """
        Optimize nuisance parameters at a single grid point. Note that no
        slicing and masking is performed to the data by this function.

        :param dispersion:
            The observed dispersion points.

        :param flux:
            The observed flux values.

        :param sigma:
            The 1-sigma uncertainty of the observed flux values.

        :param grid_index:
            The grid index.
        """

        kwds = {
            "f": lambda x, *theta: self(grid_index, x, *theta).flatten(),
            "xdata": dispersion,
            "ydata": flux,
            "sigma": sigma,
            "absolute_sigma": True,
            "p0": self.initial_guess,
        }
        kwds.update(kwargs)
        return op.curve_fit(**kwds)


    def optimize_over_all_models(self, data):
        """
        Optimize the (nuisance) model parameters over all grid points.

        :param data:
            The observed spectrum.
        """

        self._update_parameters()

        # Slice spectrum.
        obs_dispersion, obs_flux, obs_ivar = self._slice_spectrum(data)
        obs_sigma = np.sqrt(1.0/obs_ivar)

        # Apply masks.
        mask = _generate_mask(obs_dispersion, self.metadata["mask"]) \
             * np.isfinite(obs_flux * obs_sigma)

        obs_dispersion = obs_dispersion[mask]
        obs_flux = obs_flux[mask]
        obs_ivar = obs_ivar[mask]
        obs_sigma = obs_sigma[mask]


        def selector(x, *parameters):

            model_fluxes = self(None, x, *parameters)

            # Get the best point
            chi_sq = np.sum((model_fluxes - obs_flux)**2 * obs_ivar, axis=1)
            grid_index = np.argmin(chi_sq)

            print(grid_index, parameters, chi_sq[grid_index])
            return chi_sq[grid_index]


        op_parameters, cov = op.curve_fit(selector, 
            obs_dispersion, obs_flux, sigma=obs_sigma, 
            absolute_sigma=True, p0=self.initial_guess, maxfev=10000)

        # Stitch with grid parameters?

        raise a



    @property
    def initial_guess(self):
        """ Generate an uninformed guess for the model parameters. """

        defaults = {
            "c0": 1,
            "redshift": 0,
            "smoothing": 5
        }
        return tuple([defaults.get(p, 0) for p in self.parameter_names])


    def _nll(self, grid_index, theta, obs_dispersion, obs_flux, obs_ivar):
        """
        Calculate the negative log-likelihood for the data, given the model and
        parameters.
        """
        model = self(grid_index, obs_dispersion, *theta).flatten()
        return -0.5 * np.sum((model - obs_flux)**2 * obs_ivar)



    def infer(self, data, **kwargs):
        """
        Marginalize over all grid and nuisance parameters and produce posterior
        distributions on model parameters.

        """

        N_samples = kwargs.pop("N_samples", 30) # This is *per* grid point.       
        scale_uncertainties = kwargs.pop("scale_uncertainties", 5)

        self._update_parameters()

        obs_dispersion, obs_flux, obs_ivar = self._slice_spectrum(data)
        obs_sigma = np.sqrt(1.0/obs_ivar)

        # Apply masks.
        mask = _generate_mask(obs_dispersion, self.metadata["mask"]) \
             * np.isfinite(obs_flux * obs_sigma)

        obs_dispersion = obs_dispersion[mask]
        obs_flux = obs_flux[mask]
        obs_ivar = obs_ivar[mask]
        obs_sigma = obs_sigma[mask]


        N_stars, N_pixels = self.normalized_fluxes.shape
        N_theta = self.stellar_parameters.shape[1] + len(self.parameter_names)

        nlls = np.inf * np.ones((N_pixels * N_samples))
        sampled_theta = np.nan * np.ones((N_pixels * N_samples, N_theta))
        for index, stellar_parameters in enumerate(self.stellar_parameters):

            # At each grid point, optimize the parameters.
            op_theta, cov = self._optimize_nuisance_parameters(
                obs_dispersion, obs_flux, obs_sigma, index)

            scaled_cov = np.copy(cov)
            scaled_cov[np.diag_indices(op_theta.size)] *= scale_uncertainties**2
            samples = np.random.multivariate_normal(op_theta, cov, size=N_samples)

            # Calculate the likelihood at all nuisance parameter samples.
            si, ei = (index * N_samples, (index + 1) * N_samples)
            sampled_theta[si:ei, :] = np.hstack([
                np.repeat(stellar_parameters, N_samples).reshape(N_samples, -1),
                samples
            ])
            nlls[si:ei] \
                = [self._nll(index, theta, obs_dispersion, obs_flux, obs_ivar) \
                    for theta in samples]


            foo = integrate.quad(lambda v: self._nll(index, [v] + list(theta[1:]),
                obs_dispersion, obs_flux, obs_ivar), op_theta[0] - 3 * np.sqrt(cov[0,0]),
            op_theta[0] + 3 * np.sqrt(cov[0,0]))

            raise a

            # 5 sigma = -0.027855
            # 3 sigma = -0.016702

            print(index, N_pixels)

            if index > 10:
                raise a
        raise a

        

        # Now integrate +/- some reasonable limits.

        # Sum the likelihoods.

        # Look at the grou



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

        observed_sigma = np.sqrt(1.0/observed_ivar)
        finite = np.isfinite(observed_flux * observed_sigma)


        # Marginalize over grid points.
        N_pixels = self.normalized_fluxes.shape[0]
        chi_sqs = np.inf * np.ones(N_pixels)

        for grid_index in range(N_pixels):
            print(grid_index, N_pixels)

            f = lambda x, *p: self(grid_index, x, *p).flatten()

            try:
                op_parameters, cov = op.curve_fit(f, observed_dispersion[finite],
                    observed_flux[finite], sigma=observed_sigma[finite],
                    p0=self.initial_guess, absolute_sigma=True)

            except RuntimeError:
                continue


            model_flux = f(observed_dispersion[finite], *op_parameters)
            chi_sqs[grid_index] = np.sum((observed_flux[finite] - model_flux)**2 \
                * observed_ivar[finite])

        raise a


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


