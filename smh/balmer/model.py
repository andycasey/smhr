#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Fit Balmer lines. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import (integrate, interpolate, ndimage, optimize as op)
from matplotlib.ticker import MaxNLocator

# HACK REMOVE TODO
try:
    from . import utils
except ValueError:
    import utils

__all__ = ["BalmerLineModel"]


logger = logging.getLogger(__name__)

def unique_indices(a):
    b = np.ascontiguousarray(a).view(
        np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)

    return idx

class BalmerLineModel(object):

    _minimum_valid_model_file_size = 32
    _minimum_grid_points_required = 1
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

        # Check file sizes to verify they are probably usable.
        usable_paths = []
        for path in paths:
            if not os.path.exists(path):
                logger.warn("Skipping path '{}' because it does not exist."\
                    .format(path))
                continue

            elif self._minimum_valid_model_file_size >= os.stat(path).st_size:
                logger.warn(
                    "Skipping path '{}' because it is a nearly-empty file."\
                    .format(path))
                continue

            usable_paths.append(path)

        if len(usable_paths) < self._minimum_grid_points_required:
            raise ValueError("insufficient number of grid points ({} < {})"\
                .format(len(usable_paths), self._minimum_grid_points_required))

        self.paths = usable_paths

        # Get representative wavelength limits.
        wavelengths = utils.parse_spectrum(self.paths[0])[:, 0]
        self.wavelength_limits = (min(wavelengths), max(wavelengths))

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



    def __call__(self, x, model_dispersion, model_normalized_flux, *parameters):

        # Marginalize over all models.

        assert len(parameters) == len(self.parameter_names)

        # Continuum.
        O = self.metadata["continuum_order"]
        continuum = 1.0 if 0 > O \
            else np.polyval(parameters[-(O + 1):][::-1], model_dispersion)

        y = model_normalized_flux * continuum

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

        return np.interp(x, model_dispersion * (1 + z), y)




    def _slice_spectrum(self, observed_spectrum):

        idx = observed_spectrum.dispersion.searchsorted(self.wavelength_limits)

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

        :param ivar:
            The inverse variance values for the observed fluxes.

        :param grid_index:
            The grid index.
        """

        return_model_spectrum = kwargs.pop("__return_model_spectrum", False)

        # Load the spectrum from the grid.
        model_disp, model_flux = utils.parse_spectrum(self.paths[grid_index]).T

        kwds = {
            "f": lambda x, *theta: self(x, model_disp, model_flux, *theta),
            "xdata": dispersion,
            "ydata": flux,
            "sigma": ivar,
            "absolute_sigma": False,
            "p0": self.initial_guess,
            "full_output": True
        }
        kwds.update(kwargs)
        result = op.curve_fit(**kwds)

        return tuple(list(result) + [model_disp, model_flux]) \
            if return_model_spectrum else result


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

        # Any multiprocessing queue to send progress through?
        mp_queue = kwargs.pop("mp_queue", None)

        self._update_parameters()

        obs_dispersion, obs_flux, obs_ivar = self._slice_spectrum(data)
        
        # Apply masks.
        mask = _generate_mask(obs_dispersion, self.metadata["mask"]) \
             * np.isfinite(obs_flux * obs_ivar)

        if 0 in (obs_dispersion.size, mask.sum()):
            raise ValueError("no overlapping spectra with finite flux/ivar")

        obs_dispersion = obs_dispersion[mask]
        obs_flux, obs_ivar = obs_flux[mask], obs_ivar[mask]
    
        K, S = (len(self.parameter_names), len(self.paths))
        
        ln_likelihood = np.nan * np.ones(S)
        optimized_theta = np.nan * np.ones((S, K))
        covariances = np.nan * np.ones((S, K, K))

        kwds = {}
        kwds.update(kwargs)
        kwds["__return_model_spectrum"] = True
        for index, stellar_parameters in enumerate(self.stellar_parameters):

            # At each grid point, optimize the parameters.
            try:
                op_theta, cov, meta, mesg, ier, model_disp, model_flux \
                    = self._optimize_nuisance_parameters(
                        obs_dispersion, obs_flux, obs_ivar, index, **kwds)

            except RuntimeError:
                logger.exception(
                    "Could not optimize nuisance parameters at grid point {}"\
                    .format(stellar_parameters))
                if mp_queue is not None:
                    mp_queue.put([1 + index, S])

                continue

            meta.update(mesg=mesg, ier=ier)

            optimized_theta[index] = op_theta
            if cov is not None:
                covariances[index, :, :] = cov
            
            model = self(obs_dispersion, model_disp, model_flux, *op_theta)
            ln_likelihood[index] \
                = -0.5 * np.sum((model - obs_flux)**2 * obs_ivar)

            assert np.isfinite(ln_likelihood[index])

            # DEBUG show progress.
            print(index, S, dict(zip(self.parameter_names, op_theta)), ln_likelihood[index])

            if mp_queue is not None:
                mp_queue.put([1 + index, S])

        # Use a multi-Gaussian approximation of the nuisance parameters.    
        integrals = (2*np.pi)**(K/2.) * np.abs(np.linalg.det(covariances))**(-0.5)
        approximate_integrals = ~np.isfinite(integrals)
        integrals[approximate_integrals] = np.nanmin(integrals)

        # Marginalize over the grid parameters.
        posteriors = ln_likelihood * integrals
        marginalized_posteriors = posteriors / np.nansum(posteriors)

        # Save information to the class.
        self._inference_result = (marginalized_posteriors, ln_likelihood,
            integrals, optimized_theta, covariances)

        if mp_queue is not None:
            mp_queue.put(self._inference_result)

        return (marginalized_posteriors, ln_likelihood, integrals, optimized_theta,
            covariances)


    def marginalize(self, parameters):
        """
        Return posterior distributions marginalized over the given parameters.

        :param parameters:
            A tuple of parameters to marginalize over.
        """

        N = self.stellar_parameters.shape[1]
        grid_parameters = ("TEFF", "LOGG", "MH", "ALPHA_MH")[:N]

        remaining_parameters = [p for p in grid_parameters if p not in parameters]
        if not remaining_parameters:
            raise ValueError("too many parameters!")

        M = len(remaining_parameters)

        # Get the unique set of rows.
        column_indices \
            = np.array([grid_parameters.index(p) for p in remaining_parameters])

        row_indices = unique_indices(self.stellar_parameters[:, column_indices])

        points = self.stellar_parameters[row_indices][:, column_indices]
        pdf = np.nan * np.ones(len(points))

        posteriors = self._inference_result[1] * self._inference_result[2]
        for i, point in enumerate(points):
            match = np.all(
                self.stellar_parameters[:, column_indices] == point, axis=1)
            pdf[i] = np.nansum(posteriors[match]) / match.sum()

        pdf /= np.nansum(pdf)

        return (points, pdf)



    def draw_samples(self, N=1000):
        """
        Draw samples from the posterior in a very crude way.
        """

        x, pdf = model.marginalized_posteriors["TEFF"][4]



        samples = list(np.sort(np.random.uniform(x[0], x[-1], size=10)))
        while len(samples) < N:

            a, b = np.random.uniform(x[0], x[-1], size=2)

            a_index = np.argmin(np.abs(x - a))
            b_index = np.argmin(np.abs(x - b))

            if pdf[a_index] > pdf[b_index]:
                samples.append(a)

            else:
                samples.append(b)

            print("took {}".format(samples[-1]))




        raise a




    def plot_corner(self, fig=None):

        # Plot a corner distribution for all parameters.
        max_n_ticks = 5
        labels = None
        label_kwargs = {}


        K = self.stellar_parameters.shape[1]
        # Credit to DFM! TODO
        factor = 2.0
        lbdim = 0.5 * factor
        trdim = 0.2 * factor
        whspace = 0.05
        plotdim = factor * K + factor * (K - 1) * whspace
        dim = lbdim + plotdim + trdim

        if fig is None:
            fig, axes = plt.subplots(K, K, figsize=(dim, dim))
        else:
            try:
                axes = np.array(fig.axes).reshape((K, K))
            except:
                raise ValueError("Incorrect number of axes")

        lb = lbdim / dim
        tr = (lbdim + plotdim)/dim
        fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                            wspace=whspace, hspace=whspace)


        parameters = ("TEFF", "LOGG", "MH", "ALPHA_MH")[:K]
        for i, px in enumerate(parameters):

            ax = axes if K == 1 else axes[i, i]

            x, pdf = self.marginalize(set(parameters).difference([px]))
            ax.scatter(x, pdf, facecolor="k")

            ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks, prune="lower"))
            if i < K - 1:
                ax.set_xticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_xticklabels()]
            
            ax.set_yticklabels([])
            ax.set_ylim(0, ax.get_ylim()[1])

            for j, py in enumerate(parameters):
                
                ax = axes[i, j]

                if j > i:
                    ax.set_frame_on(False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    continue

                elif j == i:
                    continue
                

                xy, pdf = self.marginalize(set(parameters).difference([px, py]))
                ax.scatter(xy[:, 0], xy[:,1], c=pdf, s=80)

                ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks, prune="lower"))
                ax.yaxis.set_major_locator(MaxNLocator(max_n_ticks, prune="lower"))
                
                if i < K - 1:
                    ax.set_xticklabels([])

                else:
                    [l.set_rotation(45) for l in ax.get_xticklabels()]
                    if labels is not None:
                        ax.set_xlabel(labels[i], **label_kwargs)
                        ax.xaxis.set_label_coords(0.5, -0.3)

                if j > 0:
                    ax.set_yticklabels([])
                else:
                    [l.set_rotation(45) for l in ax.get_yticklabels()]
                    if labels is not None:
                        ax.set_ylabel(labels[i], **label_kwargs)
                        ax.yaxis.set_label_coords(-0.3, 0.5)


        raise a
        return fig




    def plot_projection(self, data):

        obs_dispersion, obs_flux, obs_ivar = self._slice_spectrum(data)
        
        # Apply masks.
        mask = _generate_mask(obs_dispersion, self.metadata["mask"]) \
             * np.isfinite(obs_flux * obs_ivar)

        if 0 in (obs_dispersion.size, mask.sum()):
            raise ValueError("no overlapping spectra with finite flux/ivar")

        #obs_dispersion = obs_dispersion[mask]
        #obs_flux, obs_ivar = obs_flux[mask], obs_ivar[mask]
        obs_flux[~mask] = np.nan

        fig, ax = plt.subplots()

        index = np.nanargmax(self._inference_result[0])
        op_theta = self._inference_result[3][index]

        model_disp, model_flux = utils.parse_spectrum(self.paths[index]).T
        y = self(obs_dispersion, model_disp, model_flux, *op_theta)
        y[~mask] = np.nan

        ax.plot(obs_dispersion, obs_flux, c='k')
        ax.plot(obs_dispersion, y, c='r')

        raise a




    def fit(self, observed_spectrum):
        """
        Fit this Balmer line model (and the nuisance parameters) to the data.

        :param observed_spectrum:
            An observed spectrum.
        """

        self._update_parameters()

        raise YoureADinasour

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






def unique_indices(a):
    b = np.ascontiguousarray(a).view(
        np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)

    return idx


if __name__ == "__main__":



    # Show results from some previous inference.



    from glob import glob

    """
    model = BalmerLineModel(
        glob("smh/balmer/models/metpoor*_bet/*.prf"), 
        mask=[
            [1000, 4841],
            [4859.51, 4859.84],
            [4860, 4862], # core.
            [4871, 4873],
            [4878, 4879],
            [4846.26, 4852],
            [4881, 10000]
        ],
        redshift=True,
        smoothing=False,
        continuum_order=3
    )
    """
    model = BalmerLineModel(glob("smh/balmer/models/metpoor*bet*/*.prf"),
        mask=[
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
        continuum_order=2)

    #import cPickle as pickle
    #with open("temp_inference.pkl", "rb") as fp:
    #    model._inference_result = pickle.load(fp)

    model = BalmerLineModel(glob("smh/balmer/models/metpoor*bet*/*.prf"),
        mask=[
            [1000, 4841],
            [4859.51, 4859.84],
            [4860, 4862], # core.
            [4871, 4873],
            [4878, 4879],
            [4881, 10000]
        ],
        redshift=True, smoothing=False, continuum_order=2)

    # Get a rest-frame normalized spectrum...
    from smh.specutils import Spectrum1D
    hd122563 = Spectrum1D.read("smh/balmer/hd122563.fits")
    hd140283 = Spectrum1D.read("hd140283.fits")
    hd140283._ivar = 1e-5 * np.ones(hd140283.dispersion.size)
    hd122563._ivar = 1e-5 * np.ones(hd122563.dispersion.size)
 
