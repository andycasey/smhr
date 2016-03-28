#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" General utilities for dealing with spectra. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["calculate_fractional_overlap", "find_overlaps", "stitch"]

import numpy as np
from scipy import interpolate

from . import spectrum

def calculate_fractional_overlap(interest_range, comparison_range):
    """
    Calculate how much of the range of interest overlaps with the comparison
    range.
    """

    if not (interest_range[-1] >= comparison_range[0] \
        and comparison_range[-1] >= interest_range[0]):
        return 0.0 # No overlap

    elif   (interest_range[0] >= comparison_range[0] \
        and interest_range[-1] <= comparison_range[-1]):
        return 1.0 # Total overlap 

    else:
        # Some overlap. Which side?
        if interest_range[0] < comparison_range[0]:
            # Left hand side
            width = interest_range[-1] - comparison_range[0]

        else:
            # Right hand side
            width = comparison_range[-1] - interest_range[0]

        return width/np.ptp(interest_range)




def find_overlaps(spectra, dispersion_range, return_indices=False):
    """
    Find spectra that overlap with the dispersion range given. Spectra are
    returned in order of how much they overlap with the dispersion range.

    :param spectra:
        A list of spectra.

    :param dispersion_range:
        A two-length tuple containing the start and end wavelength.

    :param return_indices: [optional]
        In addition to the overlapping spectra, return their corresponding
        indices.

    :returns:
        The spectra that overlap with the dispersion range provided, ranked by
        how much they overlap. Optionally, also return the indices of those
        spectra.
    """

    fractions = np.array([
        calculate_fractional_overlap(s.dispersion, dispersion_range) \
            for s in spectra])

    N = (fractions > 0).sum()
    indices = np.argsort(fractions)[::-1]
    overlaps = [spectra[index] for index in indices[:N]]

    """
    A faster, alternative method if sorting is not important:
    # http://stackoverflow.com/questions/325933/determine-whether-two-date-ranges-overlap/325964#325964    
    overlaps, indices = zip(*[(spectrum, i) \
        for i, spectrum in enumerate(spectra) \
            if  spectrum.dispersion[-1] >= dispersion_range[0] \
            and dispersion_range[-1] >= spectrum.dispersion[0]])
    """

    return overlaps if not return_indices else (overlaps, indices[:N])



def stitch(spectra, wl_ranges=None, mode='average'):
    """Stitches together two overlapping spectra.

    Inputs
    ----
    spectra : list of `Spectrum1D` objects.
        A list of spectra to stitch.

    wl_ranges : list of tuples with length two, optional
        The wavelength regions to use for stitching each spectrum. If no wavelength
        region is specified then the entire spectrum is used for stitching. If only
        wavelength region for some of the spectra are given, then None should be
        supplied when to use the entire spectrum.

    mode : str
        How to stack the spectra on a pixel-by-pixel basis. Available: average, sum, quadrature

    Returns
    ----
    stitched_spectrum : `Spectrum1D`
        A linear stitched spectrum.
    """

    raise NotImplementedError("must be refactored given new Spectrum1D obj")

    available = ('average', 'sum', 'quadrature', )
    if mode not in available:
        raise ValueError("available modes: {}".format(", ".join(available)))

    # Ensure that our spectra go from blue to red
    wl_starts = []
    wl_ends = []
    min_disp_step = 999
    for spectrum in spectra:
        wl_starts.append(spectrum.disp[0])
        wl_ends.append(spectrum.disp[-1])

        if np.min(np.diff(spectrum.disp)) < min_disp_step:
            min_disp_step = np.min(np.diff(spectrum.disp))

    sorted_indices = np.argsort(wl_starts)

    # Let's generate our new dispersion map
    new_disp = np.arange(spectra[sorted_indices[0]].disp[0], spectra[sorted_indices[-1]].disp[-1], min_disp_step)
    new_flux = np.zeros(len(new_disp))
    
    # And an array with the number of contributing flux points for each pixel point
    num_flux = np.zeros(len(new_disp))

    # Maintain headers, give preference to blue
    headers = {}
    for index in sorted_indices[::-1]:
        headers.update(spectra[index].headers)
        

    for index in sorted_indices:
        spectrum = spectra[index]


        spectrum.flux[np.where(~np.isfinite(spectrum.flux))] = 0.0

        if wl_ranges is None or wl_ranges[index] is None:
            wl_range = [spectrum.disp[0], spectrum.disp[-1]]

        else:
            wl_range = wl_ranges[index]
            

        # Get the subset dispersion map
        dispersion_indices = np.searchsorted(new_disp, wl_range)
        subset_disp = new_disp[dispersion_indices[0]:dispersion_indices[1]]

        # Interpolate the spectra
        f = interpolate.interp1d(spectrum.disp, spectrum.flux, kind='linear', copy=False,
                                 bounds_error=False, fill_value=np.nan)

        subset_flux = f(subset_disp)

        # Put the subset flux into our master arrays
        additional_flux = pow(subset_flux, 2) if mode == 'quadrature' else subset_flux

        good_indices = np.where(additional_flux > 0.)[0]

        new_flux[good_indices + dispersion_indices[0]] += additional_flux[good_indices]
        num_flux[good_indices + dispersion_indices[0]] += 1

    if mode == 'average':
        new_flux /= num_flux

    elif mode == 'quadrature':
        new_flux = pow(new_flux, 0.5)

    # Remove out non-finite values on the edge
    idx_l = np.searchsorted(np.isfinite(new_flux), True)
    idx_r = len(new_flux) - np.searchsorted(np.isfinite(new_flux)[::-1], True)

    return spectrum.Spectrum1D(new_disp[idx_l:idx_r], new_flux[idx_l:idx_r],
        headers=headers)