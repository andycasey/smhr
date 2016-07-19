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


def fill_between_steps(ax, x, y1, y2=0, h_align='mid', **kwargs):
    """
    Fill between for step plots in matplotlib.

    **kwargs will be passed to the matplotlib fill_between() function.
    """

    # If no Axes opject given, grab the current one:

    # First, duplicate the x values
    xx = x.repeat(2)[1:]
    # Now: the average x binwidth
    xstep = np.repeat((x[1:] - x[:-1]), 2)
    xstep = np.concatenate(([xstep[0]], xstep, [xstep[-1]]))
    # Now: add one step at end of row.
    xx = np.append(xx, xx.max() + xstep[-1])

    # Make it possible to chenge step alignment.
    if h_align == 'mid':
        xx -= xstep / 2.
    elif h_align == 'right':
        xx -= xstep

    # Also, duplicate each y coordinate in both arrays
    y1 = y1.repeat(2)#[:-1]
    if type(y2) == np.ndarray:
        y2 = y2.repeat(2)#[:-1]

    # now to the plotting part:
    return ax.fill_between(xx, y1, y2=y2, **kwargs)



class BalmerLineModel(object):

    _minimum_valid_model_file_size = 32
    _minimum_grid_points_required = 1
    _default_metadata = {
        "mask": [],
        "redshift": False,
        "smoothing": False,
        "continuum_order": -1,
        "bounds": {}
    }

    def __init__(self, paths, **kwargs):
        """
        A model for fitting Balmer lines.

        :param path:
            The disk location of a pre-computed Balmer line profile.
        """

        # Check file sizes to verify they are probably usable.
        self.paths = self._get_usable_paths(paths)

        # Get representative wavelength limits.
        self.wavelength_limits = self._get_wavelength_limits(self.paths[0])

        # Parse the stellar parameters from the usable paths.
        self.stellar_parameters = self._get_stellar_parameters(self.paths)

        self.metadata = {}
        self.metadata.update(self._default_metadata)
        for key, value in kwargs.items():
            if key in self.metadata:
                self.metadata[key] = value

        self._update_parameters()

        return None


    def _get_usable_paths(self, paths):
        """
        Get usable paths.

        :param paths:
            A list of Balmer-line model paths.
        """

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

        return usable_paths


    def _get_wavelength_limits(self, path):

        wavelengths = utils.parse_spectrum(self.paths[0])[:, 0]
        return (min(wavelengths), max(wavelengths))


    def _get_stellar_parameters(self, paths):
        return np.array(
            [utils.parse_stellar_parameters(path) for path in self.paths])


    def __getstate__(self):
        """ Return the model state for serialization. """

        return (self.paths, self.metadata, self._inference_result)


    def __setstate__(self, state):

        paths, metadata, inference_result = state

        self.paths = self._get_usable_paths(paths)
        self.wavelength_limits = self._get_wavelength_limits(self.paths[0])
        self.stellar_parameters = self._get_stellar_parameters(self.paths)

        self.metadata = {}
        self.metadata.update(metadata)
        self._update_parameters()
        self._inference_result = inference_result

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
        for parameter, (lower, upper) in self.metadata.get("bounds", {}).items():
            if parameter not in parameters: continue
            lower, upper = (lower, upper) if upper > lower else (upper, lower)
            self._parameter_bounds[parameter] = (lower, upper)
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


    def fitting_function(self, x, model_disp, model_flux, *parameters):

        # Return nans if outside bounds.
        for parameter_name, (lower, upper) in self.parameter_bounds.items():
            value = parameters[self.parameter_names.index(parameter_name)]
            if not (upper >= value and value >= lower):
                return np.nan * np.ones_like(x)

        return self.__call__(x, model_disp, model_flux, *parameters)


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
            "f": lambda x, *theta: \
                self.fitting_function(x, model_disp, model_flux, *theta),
            "xdata": dispersion,
            "ydata": flux,
            "sigma": ivar,
            "absolute_sigma": False,
            "p0": self.initial_guess(dispersion, flux),
            "full_output": True,
            "maxfev": 10000
        }
        kwds.update(kwargs)
        result = op.curve_fit(**kwds)

        return tuple(list(result) + [model_disp, model_flux]) \
            if return_model_spectrum else result




    def initial_guess(self, dispersion, flux):
        """ Generate an uninformed guess for the model parameters. """

        defaults = {
            "redshift": 0,
            "smoothing": 20
        }
        N = self.metadata["continuum_order"]
        if N > -1:
            defaults.update(dict(zip(
                ["c{}".format(i) for i in range(N + 1)],
                np.polyfit(dispersion, flux, N)[::-1])))

        return tuple([defaults[p] for p in self.parameter_names])


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
        
        likelihood = np.nan * np.ones(S)
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
            likelihood[index] \
                = np.sum(np.exp(-0.5 * ((model - obs_flux)**2 * obs_ivar)))

            assert np.isfinite(likelihood[index])

            # DEBUG show progress.
            print(index, S, dict(zip(self.parameter_names, op_theta)), likelihood[index])

            if mp_queue is not None:
                mp_queue.put([1 + index, S])

        # Use a multi-Gaussian approximation of the nuisance parameters.    
        integrals = (2*np.pi)**(K/2.) * np.abs(np.linalg.det(covariances))**(-0.5)
        
        # Give approximate (bad) covariances to non-finite integrals.
        approximate_integrals = ~np.isfinite(integrals)
        idx = np.nanargmin(integrals)
        covariances[approximate_integrals] = covariances[idx]
        integrals[approximate_integrals] = integrals[idx]

        posterior = likelihood * integrals
        
        # Save information to the class.
        self._inference_result \
            = (posterior, likelihood, integrals, optimized_theta, covariances)

        if mp_queue is not None:
            mp_queue.put(self._inference_result)

        return self._inference_result


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

        posteriors = self._inference_result[0]
        for i, point in enumerate(points):
            match = np.all(
                self.stellar_parameters[:, column_indices] == point, axis=1)
            pdf[i] = np.sum(posteriors[match]) / np.nansum(match)
            #print(i, point, pdf[i])

        if points.shape[1] == 1:
            xi = np.argsort(points.flatten())
            points = points[xi]
            pdf = pdf[xi]

        pdf /= np.nansum(pdf) # Scaling.

        return (points, pdf)


    def plot_likelihoods(self, axes=None):

        if axes is None:
            fig, axes = plt.subplots(1, 3)
        else:
            fig = axes[0].figure
            assert len(axes) == 3


        for i, ax in enumerate(axes):
            ax.scatter(
                self.stellar_parameters[:, 0], self.stellar_parameters[:, 1],
                c=self._inference_result[i])


        return fig


    def plot_projections_best_wrt_teff(self, spectrum, ax=None):

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        ax.plot(spectrum.dispersion, spectrum.flux)

        # For each unique value of teff, find the most likely result.
        unique_teffs = np.unique(self.stellar_parameters[:, 0])

        indices = []
        colors = (
            "#1ABC9C",
            "#2ECC71",
            "#3498DB",
            "#9B59B6",
            "#34495E",
            "#F1C40F",
            "#E67E22", 
            "#E74C3C",
            )
        posterior = self._inference_result[0]
        for i, unique_teff in enumerate(unique_teffs):

            match_indices = np.where(self.stellar_parameters[:, 0] == unique_teff)[0]

            index = match_indices[np.argmax(posterior[match_indices])]

            self.plot_projection(spectrum, ax=ax, index=index,
                c=colors[i % len(colors)], 
                label="{:.0f} ({:.0f})".format(
                    self.stellar_parameters[index, 0], index))

        ax.legend()
        return fig





    @property
    def marginalized_posteriors(self):

        interpolation_factor = 100
        N = self.stellar_parameters.shape[1]
        grid_parameters = ("TEFF", "LOGG", "MH", "ALPHA_MH")[:N]

        posteriors = {}
        for parameter in grid_parameters:

            x, pdf = self.marginalize(set(grid_parameters).difference([parameter]))

            xi = np.argsort(x.flatten())
            x, pdf = x[xi], pdf[xi]

            if x.size > 2:
                # Estimate a MAP value

                tck = interpolate.splrep(x, pdf, k=min(3, x.size - 1))
                x_interp = np.linspace(x[0], x[-1], x.size * interpolation_factor)
                pdf_interp = interpolate.splev(x_interp, tck)

                # Estimate a MAP value and +/- uncertainties.
                map_value = x_interp[pdf_interp.argmax()]

                cdf = [pdf_interp[:i].sum() for i in range(pdf_interp.size)]
                cdf = np.array(cdf) / np.max(cdf)

                # Get quartiles.
                li, ci, ui = cdf.searchsorted([0.16, 0.5, 0.84])
                neg = x_interp[li] - x_interp[ci]
                pos = x_interp[ui] - x_interp[ci]

            else:
                map_value = x.flatten()[np.argmax(pdf)]
                pos, neg = (np.nan, np.nan)

            posteriors[parameter] = (map_value, pos, neg), (x, pdf)
            continue

        return posteriors




    def draw_samples(self, N=1000):
        """
        Draw samples from the posterior in a very crude way.
        """

        raise a
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




    def plot_projection(self, data, ax=None, index=None, sample_cov=0,  **kwargs):

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        obs_dispersion, obs_flux, obs_ivar = self._slice_spectrum(data)
        
        # Apply masks.
        mask = _generate_mask(obs_dispersion, self.metadata["mask"]) \
             * np.isfinite(obs_flux * obs_ivar)

        if 0 in (obs_dispersion.size, mask.sum()):
            raise ValueError("no overlapping spectra with finite flux/ivar")

        #obs_dispersion = obs_dispersion[mask]
        #obs_flux, obs_ivar = obs_flux[mask], obs_ivar[mask]

        _ = np.where(mask)[0]
        si, ei = _[0], _[-1]
        
        # Show uncertainties.
        obs_sigma = np.sqrt(1.0/obs_ivar)
        fill_between_steps(ax,
            obs_dispersion[si:ei],
            obs_flux[si:ei] - obs_sigma[si:ei],
            obs_flux[si:ei] + obs_sigma[si:ei],
            facecolor="#AAAAAA", edgecolor="none", alpha=1)

        # Limit to the edge of what is OK.
        ax.plot(obs_dispersion[si:ei], obs_flux[si:ei], c="#444444", drawstyle="steps-mid")

        obs_flux[~mask] = np.nan
        ax.plot(obs_dispersion, obs_flux, c='k', drawstyle="steps-mid")


        # Get the MAP value.
        if index is None:
            index = np.nanargmax(self._inference_result[1])
        op_theta = self._inference_result[3][index]


        model_disp, model_flux = utils.parse_spectrum(self.paths[index]).T

        y = self(obs_dispersion, model_disp, model_flux, *op_theta)
        y[~mask] = np.nan

        c = kwargs.pop("c", "r")
        ax.plot(obs_dispersion, y, c=c, **kwargs)

        # Get the covariance matrix?
        if sample_cov > 0:
            cov = self._inference_result[4][index]
            print(np.sqrt(np.diag(cov)))

            # Sample values from the cov matrix and project them.
            draws = np.random.multivariate_normal(
                self._inference_result[3][index],
                self._inference_result[4][index],
                size=sample_cov)

            for draw in draws:
                y_draw = self(obs_dispersion, model_disp, model_flux, *draw)
                y_draw[~mask] = np.nan

                ax.plot(obs_dispersion, y_draw, c=c, alpha=10.0/sample_cov)

        # Draw fill_between in y?
        ax.set_title("Index {}: {}".format(index,
            self.stellar_parameters[index]))

        return fig




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


