#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functions for dealing with MARCS model photospheres. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import gzip
import logging
from collections import Counter

# Third party.
import numpy as np

# Module-specific.
from oracle.photospheres.interpolator import BaseInterpolator

# Create logger.
logger = logging.getLogger(__name__)


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
        

    def _spherical_or_plane_parallel(self, *point):

        point = list(point) + [0.5] # equi-spaced from plane-parallel/spherical
        neighbours = self.nearest_neighbours(point, 8) # 8 = 2**3

        sph_or_pp = self.stellar_parameters.view(float).reshape(
            len(self.stellar_parameters), -1)[:, -1]
        return np.round(np.median(sph_or_pp[neighbours]))


    def interpolate(self, *point, **kwargs):
        """ 
        Return the interpolated photospheric quantities on a common opacity
        scale.
        """

        # Try either spherical / plane parallel, and if that fails, switch.
        geometry = int(self._spherical_or_plane_parallel(*point))

        p = list(point) + [geometry]
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




def parse_filename(filename, full_output=False):
    """
    Return the basic stellar parameters from the filename.
    """

    basename = filename.split("/")[-1]
    teff = basename[1:5]
    logg = basename.split("_")[1][1:]
    feh = basename.split("_")[5][1:]
    parameters = map(float, [teff, logg, feh, int(basename[0].lower() == "s")])

    if full_output:
        names = ("effective_temperature", "surface_gravity", "metallicity",
            "is_spherical?")
        return (parameters, names)
    return parameters


def parse_photospheric_structure(filename, ndepth=56, line=25,
    columns=("lgTau5", "Depth", "T", "Pe", "Pg"), full_output=False):
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
