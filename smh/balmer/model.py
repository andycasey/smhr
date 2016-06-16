#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Fit Balmer lines. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import numpy as np
from scipy import (integrate, interpolate, ndimage, optimize as op)

import utils

__all__ = ["BalmerLineModel"]


logger = logging.getLogger(__name__)


class BalmerLineModel(object):

    _default_metadata = {
        "mask": [],
        "redshift": False,
        "smoothing": False,
        "continuum_order": -1
    }

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

        self.metadata = {}
        self.metadata.update(self._default_metadata)
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

        assert len(parameters) == len(self.parameter_names)

        # Continuum.
        O = self.metadata["continuum_order"]
        continuum = 1.0 if 0 > O \
            else np.polyval(parameters[-(O + 1):][::-1], self.wavelengths)

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
            kernel = abs(parameters[index])
            y = ndimage.gaussian_filter1d(y, kernel, axis=-1)


        # Redshift?
        try:
            index = self.parameter_names.index("redshift")
        except ValueError:
            z = 0
        else:
            v = parameters[index]
            z = v/299792.458 # km/s

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


    def _optimize_nuisance_parameters(self, dispersion, flux, ivar, grid_index,
        **kwargs):
        """
        Optimize nuisance parameters at a single grid point. Note that no
        slicing and masking is performed to the data by this function.

        :param dispersion:
            The observed dispersion points.

        :param flux:
            The observed flux values.


        :param grid_index:
            The grid index.
        """

        kwds = {
            "f": lambda x, *theta: self(grid_index, x, *theta).flatten(),
            "xdata": dispersion,
            "ydata": flux,
            "sigma": ivar,
            "absolute_sigma": False,
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
            return model_fluxes[grid_index]

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
            "smoothing": 20
        }
        return tuple([defaults.get(p, 0) for p in self.parameter_names])


    def _chi_sq(self, grid_index, theta, obs_dispersion, obs_flux, obs_ivar):
        """
        Calculate the chi_sq for the data, given the model and
        parameters.
        """
        model = self(grid_index, obs_dispersion, *theta).flatten()
        return np.sum((model - obs_flux)**2 * obs_ivar)



    def infer(self, data, **kwargs):
        """
        Marginalize over all grid and nuisance parameters and produce posterior
        distributions on model parameters.

        """

        self._update_parameters()

        obs_dispersion, obs_flux, obs_ivar = self._slice_spectrum(data)
        
        # Apply masks.
        mask = _generate_mask(obs_dispersion, self.metadata["mask"]) \
             * np.isfinite(obs_flux * obs_ivar)

        obs_dispersion = obs_dispersion[mask]
        obs_flux, obs_ivar = obs_flux[mask], obs_ivar[mask]
        
        K = len(self.parameter_names)
        S = self.normalized_fluxes.shape[0]
        
        previous_cov = None

        likelihoods = np.nan * np.ones(S)
        integrals = np.nan * np.ones(S)
        optimized_theta = np.nan * np.ones((S, K))
        covariances = np.nan * np.ones((S, K, K))

        for index, stellar_parameters in enumerate(self.stellar_parameters):

            # At each grid point, optimize the parameters.
            try:
                op_theta, cov = self._optimize_nuisance_parameters(
                    obs_dispersion, obs_flux, obs_ivar, index, maxfev=50000)

            except RuntimeError:
                logger.exception(
                    "Could not optimize nuisance parameters at grid point {}"\
                    .format(stellar_parameters))
                continue

            # TODO revisit this -- maybe we should just force bounds on params
            if np.all(np.isfinite(cov)):
                previous_cov = cov
            else:
                cov = previous_cov

            optimized_theta[index] = op_theta
            covariances[index, :, :] = cov

            model_flux = self(index, obs_dispersion, *op_theta).flatten()
            likelihoods[index] \
                = np.exp(-0.5 * np.sum((model_flux - obs_flux)**2 * obs_ivar))

            # Use a multi-Gaussian approximation of the nuisance parameters.
            integrals[index] = (2*np.pi)**(K/2.) * np.abs(np.linalg.det(cov))**(-0.5)

            # DEBUG show progress.
            print(index, S, dict(zip(self.parameter_names, op_theta)))
        

        """
        fig, ax = plt.subplots(3)
        ax[0].scatter(self.stellar_parameters[:, 0], self.stellar_parameters[:, 1],
            c=likelihoods, cmap=matplotlib.cm.afmhot)
        ax[0].set_title("likelihoods")
        ax[1].scatter(self.stellar_parameters[:,0], self.stellar_parameters[:, 1],
            c=integrals, cmap=matplotlib.cm.afmhot)
        ax[1].set_title("integrals")
        ax[2].scatter(self.stellar_parameters[:, 0], self.stellar_parameters[:,1],
            c=likelihoods * integrals, cmap=matplotlib.cm.afmhot)
        # TODO No priors here... yet
        """

        # Marginalize over the grid parameters.
        posteriors = likelihoods * integrals
        marginalized_posteriors = posteriors / np.nansum(posteriors)

        return (marginalized_posteriors, likelihoods, integrals, optimized_theta,
            covariances)

        """


        fig, axes = plt.subplots(3, 2)
        for i in range(3):
            x = np.unique(self.stellar_parameters[:,i])
            yvals = np.array([np.sum(marginalized_posteriors[self.stellar_parameters[:, i] == xi]) for xi in x])
            axes[i, 0].scatter(x, yvals)


            
            # Fit spline
            percentiles = np.array([np.sum(yvals[:j+1]) for j in range(yvals.size)])
            tck = interpolate.splrep(x, percentiles)

            axes[i, 1].scatter(x, percentiles)

            x_resampled = np.linspace(min(x), max(x), 1000)
            axes[i, 1].plot(x_resampled, intepolate.splev(x_resampled, tck), c='b')



        # Show the best thing.
        fig, ax = plt.subplots()
        best_index = np.nanargmax(normalized_metrics)

        ax.plot(obs_dispersion, obs_flux, c='k')
        ax.plot(obs_dispersion, self(best_index, obs_dispersion, *optimized_theta[best_index]).flatten(),
            c='r')


        return marginalized_posteriors

        """


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

    model = BalmerLineModel(glob("data/gam_*.prf"), mask=[
        [1000, 4328],
        [4336.84, 4337.31],
        [4337.48, 4338.09],
        [4339.31, 4339.53],
        [4340.03, 4340.98],
        [4341.2, 4341.59],
        [4344.11, 4344.72],
        [4350, 10000]
    ],
    redshift=True,
    smoothing=False,
    continuum_order=2
    )

    """
    model = BalmerLineModel(paths,
        mask=[
            [1000, 4841],
            [4859.51, 4859.84],
            [4860, 4862], # core.
            [4871, 4873],
            [4878, 4879],
            [4881, 10000]
        ],
        redshift=True, smoothing=False, continuum_order=2)
    """

    # Get a rest-frame normalized spectrum...
    from smh.specutils import Spectrum1D
    hd122563 = Spectrum1D.read("hd122563.fits")
    hd140283 = Spectrum1D.read("../../hd140283.fits")
    hd140283._ivar = 1e-5 * np.ones(hd140283.dispersion.size)
    hd122563._ivar = 1e-5 * np.ones(hd122563.dispersion.size)
 
