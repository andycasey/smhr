#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functions for dealing with Castelli-Kurucz model photospheres. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import gzip
import logging
import warnings

# Third party.
import numpy as np

# Module-specific.
from oracle.photospheres.interpolator import BaseInterpolator

# Create logger.
logger = logging.getLogger(__name__)

# Warn about missing variance arrays, but only once.
class StandardCompositionAssumed(Warning):
    pass
warnings.simplefilter("once", StandardCompositionAssumed)


class Interpolator(BaseInterpolator):

    logarithmic_photosphere_quantities = ["P", "XNE"]

    def __init__(self, **kwargs):
        return super(self.__class__, self).__init__("castelli-kurucz-2004.pkl",
            **kwargs)


    def interpolate(self, *point, **kwargs):
        """ 
        Return the interpolated photospheric quantities on a common opacity
        scale.
        """

        # Assume zero alpha enhancement if not given.
        if len(point) == 3:
            point = [] + list(point) + [0]
            warnings.warn("Assuming standard [alpha/Fe] = 0 composition unless "
                "otherwise specified.", StandardCompositionAssumed)

        return super(self.__class__, self).interpolate(*point, **kwargs)


def parse_filename(filename, full_output=False):
    """
    Return the basic stellar parameters from the filename.
    """

    basename = filename.split("/")[-1]
    teff = float(basename.split("t")[1].split("g")[0])
    logg = float(basename.split("g")[1].split("k")[0])/10.
    feh = float(basename[1:4].replace("p", "").replace("m", "-"))/10.
    alpha = [0, 0.4][basename[4] == "a"]
    parameters = [teff, logg, feh, alpha]

    if full_output:
        names = ("effective_temperature", "surface_gravity", "metallicity",
            "alpha_enhancement")
        return (parameters, names)
    return parameters


def parse_photospheric_structure(filename, ndepth=None, line=23,
    columns=("RHOX", "T", "P", "XNE", "ABROSS"), full_output=False):
    """
    Parse the photospheric structure (optical depths, temperatures, electron
    and gas pressures) from the filename provided.
    """

    opener = gzip.open if filename[-3:].lower() == ".gz" else open
    with opener(filename, "r") as fp:
        contents = fp.readlines()
    if ndepth is None: # Read it from file if it was not given
        ndepth = int(contents[22].split()[2])

    all_columns = ("RHOX", "T", "P", "XNE", "ABROSS", "ACCRAD", "VTURB",
        "FLXCNV", "VCONV", "VELSND")
    data = np.array(map(float, 
        "".join(contents[line:line+ndepth]).split())).reshape(ndepth, -1)
    indices = np.array([all_columns.index(c) for c in columns])

    data = data[:, indices]

    if full_output:
        return (data, columns)
    return data



