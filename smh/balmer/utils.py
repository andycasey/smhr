#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Utility functions. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
import numpy as np

__all__ = ["parse_spectrum", "parse_stellar_parameters"]


def parse_stellar_parameters(path):
    """
    Parse the stellar parameters from the path name.

    :param path:
        The path containing computed spectra:
        (e.g. "alf_synth_t5500g2.90m-2.00x1.5.prf")
    """

    basename = os.path.basename(path).split("_")[-1]
    parent_folder = path.split("/")[-2]
    
    teff, logg, mh = (float(each) for each in \
        (basename[1:5],  basename[6:10], basename[11:16].rstrip("x")))
    alpha_mh = 0.4 if "alpha04" in parent_folder.lower() else 0.0

    return (teff, logg, mh, alpha_mh)



def parse_spectrum(path):
    """
    Parse the spectrum from the provided path.

    :param path:
        The path containing computed spectra.
    """

    with open(path, "r") as fp:
        contents = fp.read().split("\n")[3:-3]

    return np.array([point.split() for point in contents], dtype=float)
