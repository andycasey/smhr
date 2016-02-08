#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Interpolate model photospheres. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import os
import logging
import cPickle as pickle
from pkg_resources import resource_stream

# Third-party.
import astropy.table
import numpy as np
import scipy.interpolate as interpolate
from scipy import __version__ as scipy_version

from .photosphere import Photosphere

major, minor = map(int, str(scipy_version).split(".")[:2])
has_scipy_requirements = (major > 0 or minor >= 14)

# Create logger.
logger = logging.getLogger(__name__)

# Ignore divide by warnings.
np.seterr(divide="ignore", invalid="ignore")

class BaseInterpolator(object):
    
    opacity_scale = None
    logarithmic_photosphere_quantities = []

    def __init__(self, pickled_photospheres, neighbours=30, method="linear",
        rescale=True, live_dangerously=True):
        """
        Create a class to interpolate photospheric quantities.

        :param pickled_photospheres:
            The kind of photospheres to interpolate. 

        :type pickled_photospheres:
            str
        """

        if os.path.exists(pickled_photospheres):
            with open(pickled_photospheres, "rb") as fp:
                _ = pickle.load(fp)

        else:
            try:
                with resource_stream(__name__, pickled_photospheres) as fp:
                    _ = pickle.load(fp)
            except:
                raise ValueError("photosphere filename '{}' does not exist"\
                    .format(pickled_photospheres))

        stellar_parameters, photospheres, photospheric_quantities, meta = _

        # Look for duplicate stellar parameter rows.
        if stellar_parameters.dtype.names is None:
            raise ValueError("no stellar parameter names given -- the pickled "
                "stellar parameters are expected to be a record array")
        array_view = stellar_parameters.view(float).reshape(
            stellar_parameters.size, -1)
        _ = np.ascontiguousarray(array_view).view(np.dtype((np.void,
            array_view.dtype.itemsize * array_view.shape[1])))
        _, idx = np.unique(_, return_index=True)

        if idx.size != stellar_parameters.size:
            raise ValueError("{} duplicate stellar parameters found".format(
                stellar_parameters.size - idx.size))

        self.live_dangerously = live_dangerously
        self.stellar_parameters = stellar_parameters
        self.photospheres = photospheres
        self.photospheric_quantities = photospheric_quantities
        self.meta = meta
        self.method = method
        self.rescale = rescale
        if self.rescale and not has_scipy_requirements:
            logger.warn("scipy >= 0.14.0 is required for auto-rescaling points "
                "before interpolation")
            self.rescale = False
        self.neighbours = neighbours

        # Create unique copies of stellar parameters for faster access
        names = stellar_parameters.dtype.names
        self._stellar_parameters = dict(zip(names,
            [np.unique(stellar_parameters[name]) for name in names]))
        self._boundaries = \
            [(stellar_parameters[name].min(), stellar_parameters[name].max()) \
                for name in names]

    def __call__(self, *args, **kwargs):
        """ Alias to Interpolator.interpolate """
        return self.interpolate(*args, **kwargs)


    def _return_photosphere(self, stellar_parameters, quantities):
        """ 
        Prepare the interpolated photospheric quantities (with correct columns,
        units, metadata, etc).
        """

        meta = self.meta.copy()
        meta["common_optical_scale"] = self.opacity_scale
        meta["stellar_parameters"] = \
            dict(zip(self.stellar_parameters.dtype.names, stellar_parameters))
        units = meta.pop("photospheric_units", None)

        photosphere = Photosphere(data=quantities, meta=meta,
            names=self.photospheric_quantities)
        if units is not None:
            for name, unit in zip(self.photospheric_quantities, units):
                photosphere[name].unit = unit
        return photosphere


    def nearest_neighbours(self, point, n):
        """
        Return the indices of the n nearest neighbours to the point.
        """

        stellar_parameters = _recarray_to_array(self.stellar_parameters)
        distances = np.sum(((point - stellar_parameters) \
            / np.ptp(stellar_parameters, axis=0))**2, axis=1)
        return distances.argsort()[:n]


    def nearest(self, *point):
        logger.warn("Living dangerously!")
        return self._return_photosphere(point, 
            self.photospheres[self.nearest_neighbours(point, 1)[0]])


    def _test_interpolator(self, *point):
        """
        Interpolate a photosphere to the given point, and ignore the nearest
        stellar model.
        """
        return self.interpolate(*point, __ignore_nearest=True)


    def interpolate(self, *point, **kwargs):
        """
        Interpolate the photospheric structure at the given stellar parameters.
        """

        # Is the point actually within the grid?
        point = np.array(point)

        # point will contain: effective temperature, surface gravity, metallicity
        if 0 >= point[0]:
            raise ValueError("effective temperature must be positive")

        __ignore_nearest = kwargs.pop("__ignore_nearest", False)

        grid = self.stellar_parameters.view(float).reshape(
            len(self.stellar_parameters), -1)
        grid_index = np.all(grid == point, axis=1)
        if np.any(grid_index) and not __ignore_nearest:
            grid_index = np.where(grid_index)[0][0]
            return self._return_photosphere(point, self.photospheres[grid_index])

        # Work out what the optical depth points will be in our (to-be)-
        # interpolated photosphere.
        if __ignore_nearest:
            neighbours = self.nearest_neighbours(point, self.neighbours + 1)[1:]
        else:
            neighbours = self.nearest_neighbours(point, self.neighbours)
        stellar_parameters = _recarray_to_array(self.stellar_parameters)

        # Shapes required for griddata:
        # points: (N, ndim)
        # values: shape (N, )
        # xi: shape (M, ndim)

        # Protect Qhull from columns with a single value.
        cols = _protect_qhull(stellar_parameters[neighbours])
        shape = [self.neighbours] + list(self.photospheres.shape[1:])
        if self.opacity_scale is not None:
            opacity_index = self.photospheric_quantities.index(self.opacity_scale)
            kwds = {
                "xi": point[cols].reshape(1, len(cols)),
                "points": stellar_parameters[neighbours][:, cols],
                "values": self.photospheres[neighbours, :, opacity_index],
                "method": self.method,
                "rescale": self.rescale
            }
            common_opacity_scale = interpolate.griddata(**kwds)

            if np.all(~np.isfinite(common_opacity_scale)):
                if self.live_dangerously: return self.nearest(*point)
                raise ValueError("cannot interpolate {0} photosphere at {1}"\
                    .format(self.meta["kind"], point))

            # At the neighbouring N points, create splines of all the values
            # with respect to their own opacity scales, then calcualte the 
            # photospheric quantities on the common opacity scale.

            #photospheres.shape = (N_model, N_depth, N_quantities)
            neighbour_quantities = np.zeros(shape)
            for i, neighbour in enumerate(neighbours):
                neighbour_quantities[i, :, :] = \
                    resample_photosphere(common_opacity_scale,
                        self.photospheres[neighbour, :, :], opacity_index)

        else:
            opacity_index = None
            common_opacity_scale = np.arange(self.photospheres.shape[1])
            neighbour_quantities = self.photospheres[neighbours, :, :]

        # Logify/unlogify any quantities.
        for quantity in self.logarithmic_photosphere_quantities:
            try:
                index = self.photospheric_quantities.index(quantity)
            except ValueError:
                logger.warn("Could not find logarithmic photospheric quantity "\
                    "'{}'".format(quantity))
                continue
            else:
                neighbour_quantities[:, :,  index] = \
                    np.log10(neighbour_quantities[:, :, index])

        # Now interpolate the photospheric quantities.
        kwds = {
            "xi": point[cols].reshape(1, len(cols)),
            "points": stellar_parameters[neighbours][:, cols],
            "values": neighbour_quantities,
            "method": self.method,
            "rescale": self.rescale
        }
        interpolated_quantities = interpolate.griddata(**kwds).reshape(shape[1:])

        if np.all(~np.isfinite(interpolated_quantities)):
            if self.live_dangerously: return self.nearest(*point)
            raise ValueError("cannot interpolate {0} photosphere at {1}".format(
                self.meta["kind"], point))

        # Logify/unlogify any quantities.
        for quantity in self.logarithmic_photosphere_quantities:
            try:
                index = self.photospheric_quantities.index(quantity)
            except ValueError:
                continue
            else:
                interpolated_quantities[:, index] = \
                    10**interpolated_quantities[:, index]

        return self._return_photosphere(point, interpolated_quantities)


def resample_photosphere(opacities, photosphere, opacity_index):
    """ Resample photospheric quantities onto a new opacity scale. """

    if opacity_index is None:
        return photosphere

    resampled_photosphere = np.zeros(photosphere.shape)
    n_quantities = photosphere.shape[1]
    for i in range(n_quantities):
        if i == opacity_index: continue
        # Create spline function.
        tk = \
            interpolate.splrep(photosphere[:, opacity_index], photosphere[:, i])

        # Evaluate photospheric quantities at the new opacities
        resampled_photosphere[:, i] = interpolate.splev(opacities.flatten(), tk)

    # Update photosphere with new opacities
    resampled_photosphere[:, opacity_index] = opacities
    return resampled_photosphere


def _recarray_to_array(a, dtype=float):
    return a.view(dtype).reshape(len(a), -1)

def _protect_qhull(a):
    return np.where([np.unique(a[:, i]).size > 1 for i in range(a.shape[1])])[0]

"""
if __name__ == "__main__":


    import matplotlib.pyplot as plt

    point = [5777., 4.445, 0.0]

    marcs = Interpolator()
    sun_interp = marcs.interpolate(*point)
    neighbours = marcs.neighbours(*point)

    fig, axes = plt.subplots(4)


    quantities = ("Depth", "T", "Pe", "Pg")
    for i, (ax, quantity) in enumerate(zip(axes, quantities)):
    
        quantity_ranges = marcs.photospheres[neighbours, :, i+1]
        min_range = np.min(quantity_ranges, axis=0)
        max_range = np.max(quantity_ranges, axis=0)

        ax.fill_between(sun_interp["logTau5000"], min_range, max_range,
            facecolor="#cccccc")
        ax.plot(sun_interp["logTau5000"], min_range, c="#666666")
        ax.plot(sun_interp["logTau5000"], max_range, c="#666666")
        
        ax.plot(sun_interp["logTau5000"], sun_interp[quantity], 'k', lw=2)

        ax.set_xlabel("tau(5000)")
        ax.set_ylabel(quantity)

        if i > 1:
            ax.set_yscale('log')

    print(sun_interp)

    plt.show()
    
"""


