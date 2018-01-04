#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" General utilities for dealing with spectra. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["calculate_fractional_overlap", "find_overlaps"]

import numpy as np
import time


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
        return width/np.ptp(interest_range) # Fractional overlap


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

def calculate_noise(y, window, method='madstd', verbose=True, full_output=False):
    """
    Calculate SNR a posteriori.
    Very inefficient but who cares!
    :param y:
        array of flux values to calculate SNR.
        Assumes equal or slowly changing x sampling
    :param window:
        pixel window to use
    :param method:
        madstd (default): use astropy.stats.mean_absolute_deviation
        betasigma: use PyAstronomy BSEqSamp
    :param verbose:
        if True (default), prints how long it takes to calculate the SNR.
    """
    assert method in ['madstd', 'betasigma']

    from scipy import signal
    ## Signal
    S = signal.medfilt(y, window)
    
    ## Noise Estimator
    if method == 'betasigma':
        from PyAstronomy import pyasl
        beq = pyasl.BSEqSamp()
        N = 1; j = 1
        def noise_estimator(x):
            return beq.betaSigma(x, N, j, returnMAD=True)[0]
    elif method == 'madstd':
        try:
            from astropy.stats import mean_absolute_deviation
        except ImportError:
            def noise_estimator(x):
                return np.nanmedian(np.abs(x-np.nanmedian(x)))
        else:
            noise_estimator = lambda x: mean_absolute_deviation(x, ignore_nan=True)
    
    ## Calculate noise
    noise = np.zeros(len(y))
    if verbose: start = time.time()
    # middle: calculate rolling value
    for i in range(window, len(noise)-window):
        noise[i] = noise_estimator(y[i-window:i+window])
    # edges: just use the same noise value
    noise[:window] = noise_estimator(y[:window])
    noise[-window:] = noise_estimator(y[-window:])
    
    if verbose:
        print("Noise estimate of {} points with window {} took {:.1f}s".format(len(noise), window, time.time()-start))

    # https://en.wikipedia.org/wiki/Median_absolute_deviation
    # sigma = 1.4826 * MAD
    noise *= 1.4826
    
    if full_output:
        SNR = S/noise
        return noise, S, SNR
    return noise

def calculate_snr(*args, **kwargs):
    kwargs = kwargs.copy()
    kwargs["full_output"] = True
    noise, S, SNR = calculate_noise(*args, **kwargs)
    return SNR
