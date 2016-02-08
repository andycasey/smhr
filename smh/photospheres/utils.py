#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Atmosphere utilities. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"


def estimate_microturbulence(effective_temperature, surface_gravity):
    """
    Estimate microtubulence from relations between effective temperature and
    surface gravity. For giants (logg < 3.5) the relationship employed is from
    Kirby et al. (2008, ) and for dwarfs (logg >= 3.5) the Reddy et al. (2003)
    relation is used.

    :param effective_temperature:
        The effective temperature of the star in Kelvin.

    :type effective_temperature:
        float

    :param surface_gravity:
        The surface gravity of the star.

    :type surface_gravity:
        float

    :returns:
        The estimated microturbulence (km/s) from the given stellar parameters.

    :rtype:
        float
    """

    if surface_gravity >= 3.5:
        return 1.28 + 3.3e-4 * (effective_temperature - 6000) \
            - 0.64 * (surface_gravity - 4.5)
    else:
        return 2.70 - 0.509 * surface_gravity