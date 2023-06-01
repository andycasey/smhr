#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functions for dealing with MARCS model photospheres. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import gzip
import logging
import warnings

# Third party.
import numpy as np

# Module-specific.
from .interpolator import BaseInterpolator

# Create logger.
logger = logging.getLogger(__name__)

# Warn about missing variance arrays, but only once.
class StandardCompositionAssumed(Warning):
    pass
warnings.simplefilter("once", StandardCompositionAssumed)

class Interpolator(BaseInterpolator):

    opacity_scale = "lgTau5"
    logarithmic_photosphere_quantities = ["Pe", "Pg"]

    def __init__(self, **kwargs):
        """
        A class to interpolate spherical and plane-parallel MARCS model
        photospheres.

        We use the standard-composition 1 Solar mass models with microturbulence
        of 1 km/s in plane-parallel models and 2 km/s in spherical models.

        """
        return super(self.__class__, self).__init__("marcs-2011-standard.pkl",
            **kwargs)
        

    @property
    def stellar_parameters_grid(self):
        """
        Overload this function so that the radii are considered 0 or 1 for the
        purposes of interpolation.
        """

        grid = self.stellar_parameters.view(float).reshape(
            len(self.stellar_parameters), -1).copy()
        grid[:, -1] = np.clip(grid[:, -1], 0, 1)
        return grid


    def _spherical_or_plane_parallel(self, *point):

        point = list(point) + [0.5] # equi-spaced from plane-parallel/spherical
        neighbours = self.nearest_neighbours(point, 8) # 8 = 2**3
        return np.round(np.median(self.stellar_parameters_grid[neighbours, -1]))


    def interpolate(self, *point, **kwargs):
        """ 
        Return the interpolated photospheric quantities on a common opacity
        scale.
        """

        if len(point) > 3:
            point = list(point)
            warnings.warn(
                "Dropping alpha from interpolation")
            point = point[0:3]
        
        # Try either spherical / plane parallel, and if that fails, switch.
        geometry = int(self._spherical_or_plane_parallel(*point))

        p = list(point) + [geometry]
        print("p",p)
        try:
            return super(self.__class__, self).interpolate(*p, **kwargs)

        except ValueError:
            # Ok, switch geometry and try again.
            new_geometry = (1, 0)[geometry > 0]
            human_name = ["plane-parallel", "spherical"] # 1 = spherical
            logger.exception("Failed to produce photosphere with {0} geometry "\
                "at parameters {1}. Trying a {2} photosphere..".format(
                    human_name[geometry], point, human_name[new_geometry]))

            p = list(point) + [new_geometry]
            return super(self.__class__, self).interpolate(*p, **kwargs)



    def _prepare_photosphere(self, stellar_parameters, quantities, 
        neighbour_indices):
        """ 
        Prepare the interpolated photospheric quantities (with correct columns,
        units, metadata, etc).
        """

        radius = np.mean(self.stellar_parameters["radius"][neighbour_indices])
        return self._return_photosphere(
            stellar_parameters, quantities, meta=dict(radius=radius))




def parse_filename(filename, full_output=False):
    """
    Return the basic stellar parameters from the filename.
    """

    basename = filename.split("/")[-1]
    teff = basename[1:5]
    logg = basename.split("_")[1][1:]
    feh = basename.split("_")[5][1:]

    is_spherical = int(basename[0].lower() == "s")
    if is_spherical:
        # Need to open the file and get the radius from line 8
        f = gzip.open if filename.lower().endswith(".gz") else open
        with f(filename, "r") as fp:
            content = fp.readlines()
        radius = content[7].split()[0]

    else:
        radius = 0

    parameters = np.array([teff, logg, feh, radius]).astype(float)

    if full_output:
        names = ("effective_temperature", "surface_gravity", "metallicity",
            "radius")
        return (parameters, names)
    return parameters


def parse_photospheric_structure(filename, ndepth=56, line=25,
    columns=("lgTau5", "Depth", "T", "Pe", "Pg", "Prad"), full_output=False):
    """
    Parse the photospheric structure (optical depths, temperatures, electron
    and gas pressures) from the filename provided.
    """

    opener = gzip.open if filename[-3:].lower() == ".gz" else open
    with opener(filename, "r") as fp:
        contents = fp.readlines()

    all_columns = ["k", "lgTauR", "lgTau5", "Depth", "T", "Pe", "Pg", "Prad",
        "Pturb"]
    data = np.array(map(float, 
        "".join(contents[line:line+ndepth]).split())).reshape(ndepth, -1)

    # Splice by the columns we want.
    indices = np.array([all_columns.index(c) for c in columns])
    data = data[:, indices]

    if full_output:
        return (data, columns)
    return data
