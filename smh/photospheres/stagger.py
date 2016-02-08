#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functions for dealing with the Stagger model photospheres. """

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

    opacity_scale = "logtau"

    def __init__(self, filename, **kwargs):
        return super(self.__class__, self).__init__(filename, **kwargs)


def pickle_from_tsv_file(filename, depth_scale="optical", skiprows=72,
    delimiter=";"):
    """
    Pickle the Stagger-grid models from TSV-formatted filename.

    :param filename:
        The path of the TSV-formatted file.

    :type filename:
        str

    :param depth_scale: [optional, optical assumed]
        Which horizontal averaging method to use. Available options are:
        optical, mass density, Rosseland, or geometric height

    :type depth_scale:
        str

    :param skiprows: [optional]
        The number of rows at the top of the file before the header information.

    :type skiprows:
        int

    :param delimiter: [optional]
        The delimiting character between columns.

    :type delimiter:
        str
    """

    depth_scale_hint = depth_scale.lower()[0] # work it out from the first letter
    if depth_scale_hint not in ("o", "m", "r", "z", "g", "h"): # zgh are the same 
        raise ValueError("depth scale expected to be 'optical', 'mass density',"
            " Rosseland, or geometric height")
    if depth_scale_hint in ("g", "h"):
        depth_scale_hint = "z"
    elif depth_scale_hint == "r":
        depth_scale_hint = "R"

    depth_scale = {
        "o": "optical",
        "m": "mass density",
        "R": "Rossland opacity",
        "z": "geometric height",
    }[depth_scale_hint]

    with open(filename, "r") as fp:
        contents = fp.readlines()[skiprows + 1:]
    if contents[-1] == "\n": contents.pop(-1)

    # Number of extra columns in each row.
    n = 4

    # First three lines are for headers
    names = contents[0].strip().split(delimiter)
    units = contents[1].strip().split(delimiter)
    contents = contents[3:]

    num_models = len(set([row.split(delimiter)[n - 1] for row in contents]))
    parameters = np.nan * np.ones((num_models, n - 1))

    # Assume they all have the same number of depth points.
    assert (len(contents) % num_models) == 0
    num_depth_points = int(len(contents) / num_models)
    num_photospheric_quantitites = len(names) - n
    photospheres = np.nan * np.ones(
        (num_models, num_depth_points, num_photospheric_quantitites))

    for i in range(num_models):
        # The '4:' arises from the first four columns being the model parameters
        parameters[i, :] = \
            map(float, contents[i*num_depth_points].split(delimiter)[:n-1])
        photospheres[i, :, :] = np.array(
            [map(float, map(str.strip, _.split(delimiter)[n:])) \
                for _ in contents[i*num_depth_points:(i + 1)*num_depth_points]])

    names, units = names[n:], units[n:]
    # Replace dimensionless columns with "" for astropy.
    
    # Which depth scale do we want?
    indices = np.array([0] + [i for i, name in enumerate(names) \
        if name.endswith("({})".format(depth_scale_hint))])
    names = [names[0]] + [names[i][:-3] for i in indices[1:]]

    units = [units[i].replace("[-]", "") for i in indices]
    photospheres = photospheres[:, :, indices]

    meta = {
        "kind": "Stagger",
        "source_path": filename,
        "horizontal_averaging": depth_scale,
        "photospheric_units": units
    }
    assert np.all(np.isfinite(parameters))
    assert np.all(np.isfinite(photospheres))

    parameters = np.core.records.fromarrays(parameters.T,
        names=("effective_temperature", "surface_gravity", "metallicity"))

    return (parameters, photospheres, names, meta)


