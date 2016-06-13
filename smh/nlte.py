#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Approximate corrections to non-local thermodynamic equilibrium conditions. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import os

from astropy.table import Table
from glob import glob
from scipy import interpolate


def parse_stellar_parameters(path):
    """
    Parse the element and stellar parameters from a pre-computed set of non-LTE
    corrections.

    :param path:
        The filename path (e.g., "EW_out.fe.nlte.marcs_4000_+1.0_-2.00_00.MULTI")

    :returns:
        A four-length tuple containing the temperature, surface gravity,
        metallicity, and element.
    """

    basename = os.path.basename(path)
    element = basename.split(".")[1].title()

    _, __, teff, logg, feh, ___ = basename.split("_")
    return (float(teff), float(logg), float(feh), element)


def parse_transitions(path, fast=False, **kwargs):
    """
    Parse transitions from a pre-computed file of non-LTE corrections.

    :param fast: [optional]
        If `True`, `numpy.loadtxt` is used to load the data rather than 
        returning back an `astropy.Table`.
    """

    # Load in the transitions.
    if fast:
        return np.loadtxt(path, **kwargs)

    return Table.read(path, names=("wavelength", "expot", "g_lower_level", 
        "EW_nlte", "EW_lte", "EW_delta", "departure_coefficient"), format="ascii")


class ApproximateNLTECorrectionsBase(object):
    pass


class ApproximateNLTECorrections(ApproximateNLTECorrectionsBase):

    _minimum_paths = 1

    def __init__(self, precomputed_corrections_wildmask):
        """
        Initialize an object to approximate (interpolate) between grids of pre-
        computed non-LTE corrections.

        :param precomputed_corrections_wildmask:
            A wildmask for paths matching files that contain non_LTE corrections
            to include in this class.
        """

        self._wildmask = precomputed_corrections_wildmask
        self._paths = glob(self._wildmask)

        if len(self._paths) < self._minimum_paths:
            raise ValueError(
                "could not find enough computed corrections: ({} < {})".format(
                    len(self._paths), self._minimum_paths))

        # Load the gridpoints
        self._grid = []
        self._element = []

        for path in self._paths:
            teff, logg, mh, element = parse_stellar_parameters(path)
            self._grid.append([teff, logg, mh])
            self._element.append(element)

        if len(set(self._element)) > 1:
            raise ValueError(
                "separate ApproximateNLTECorrections classes must be set up for"
                " each element")

        self._element, self._grid = self._element[0], np.array(self._grid)

        # Load the transitions from one file.
        self._transitions = parse_transitions(self._paths[0])

        # Get the EW_delta for each transition.
        self._ew_delta = np.array([
            parse_transitions(path, fast=True, usecols=(5, )) \
            for path in self._paths])

        return None



    def match_transition(self, wavelength, expot, wavelength_tolerance=0.01, 
        expot_tolerance=0.01, **kwargs):

        return \
            (np.abs(self._transitions["wavelength"] - wavelength) < wavelength_tolerance) \
          * (np.abs(self._transitions["expot"] - expot) < expot_tolerance)



    def neighbours(self, effective_temperature, surface_gravity, metallicity, N,
        scales=None):
        """
        Return indices of the `N`th-nearest neighbours in the grid. The three
        parameters are scaled by the peak-to-peak range in the grid, unless
        `scales` are indicates.

        :param effective_temperature:
            The effective temperature of the star.

        :param surface_gravity:
            The surface gravity of the star.

        :param metallicity:
            The metallicity of the star.

        :param N:
            The number of neighbouring indices to return.

        :returns:
            An array of length `N` that contains the indices of the closest
            neighbours in the grid.
        """

        point = np.array([effective_temperature, surface_gravity, metallicity])
        if scales is None:
            scales = np.ptp(self._grid, axis=0)

        distance = np.sum(((self._grid - point)/scales)**2, axis=1)

        return np.argsort(distance)[:N]


    def __call__(self, effective_temperature, surface_gravity, metallicity,
        wavelength, excitation_potential, equivalent_width=None, **kwargs):
        """
        Calculate the interpolated non-LTE correction for a single transition
        given the stellar parameters.

        :param effective_temperature:
            The effective temperature of the star.

        :param surface_gravity:
            The surface gravity of the star.

        :param metallicity:
            The overall metallicity of the star.

        :param wavelength:
            The wavelength of the transition (in Angstroms).

        :param excitation_potential:
            The excitation potential of the transition (in eV).

        :param equivalent_width: [optional]
            The measured equivalent width (in milliAngstroms).

        """

        # Is this transition in our grid?
        matches = self.match_transition(wavelength, excitation_potential,
            **kwargs)

        if not any(matches):
            # No correction available.
            return np.nan

        assert 2 > matches.sum()

        ews = self._ew_delta[:, matches]

        # Find N closest in the grid to interpolate.
        N = kwargs.pop("N", 30)
        neighbouring_indices = self.neighbours(
            effective_temperature, surface_gravity, metallicity, N)

        raise NotImplementedError


    def interpolate(self, stellar_photosphere, transitions):
        """
        Interpolate non-LTE corrections to the equivalent widths of many
        transitions for a single stellar photosphere.

        :param stellar_photosphere:
            A stellar atmosphere model.

        :param transitions:
            A table of atomic transitions.
        """


        # A convenience function.
        raise NotImplementedError



if __name__ == "__main__":

    foo = ApproximateNLTECorrections("EW*")
