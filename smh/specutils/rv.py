#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Functions related to radial velocity measurement and correction. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["cross_correlate"]

import logging
import numpy as np
from scipy import interpolate
from scipy.optimize import leastsq

from . import spectrum

logger = logging.getLogger(__name__)


def cross_correlate(observed_spectrum, template_spectrum, dispersion_range=None,
    apodize=0, resample="template"):
    """
    Cross-correlate the observed spectrum against a rest-frame template spectrum
    and measure the radial velocity of the source.

    :param observed_spectrum:
        The observed spectrum.

    :type observed_spectrum:
        :class:`specutils.Spectrum1D`

    :param template_spectrum:
        A rest-frame template spectrum.

    :type template_spectrum:
        :class:`specutils.Spectrum1D`

    :param dispersion_range: [optional]
        A two-length tuple containing the start and end wavelengths to consider.
        If not provided, then the overlap of the `observed_spectrum` and the
        `template_spectrum` will be used.

    :param apodize: [optional]
        The fraction of pixels to apodize on either side of the spectrum.

    :param resample: [optional]
        Whether to resample the 'template' onto the observed spectrum, or
        resample the 'observed' spectrum onto the template.

    :returns:
        The radial velocity, uncertainty in radial velocity, and the CCF.
    """

    if not isinstance(observed_spectrum, spectrum.Spectrum1D):
        raise TypeError(
            "observed_spectrum must be a `specutils.Spectrum1D` object")

    if not isinstance(template_spectrum, spectrum.Spectrum1D):
        raise TypeError(
            "template_spectrum must be a `spectuils.Spectrum1D` object")

    if dispersion_range is None:
        # Use the common ranges.
        dispersion_range = (
            np.max([
                observed_spectrum.dispersion[0],
                template_spectrum.dispersion[0]
            ]),
            np.min([
                observed_spectrum.dispersion[-1],
                template_spectrum.dispersion[-1]
            ])
        )

    if not isinstance(dispersion_range, (tuple, list, np.ndarray)) \
    or len(dispersion_range) != 2:
        raise TypeError("wavelength region must be a two length list-type")

    if apodize != 0:
        raise NotImplementedError("apodization not implemented yet")
        
    resample = resample.lower()
    if resample == "template":
        idx = np.searchsorted(observed_spectrum.dispersion, dispersion_range)
        finite = np.isfinite(observed_spectrum.flux[idx[0]:idx[1]])

        dispersion = observed_spectrum.dispersion[idx[0]:idx[1]][finite]
        observed_flux = observed_spectrum.flux[idx[0]:idx[1]][finite]

        func = interpolate.interp1d(
            template_spectrum.dispersion, template_spectrum.flux,
            bounds_error=False, fill_value=0.0)
        template_flux = func(dispersion)

    elif resample == "observed":
        raise NotImplementedError("why would you do this?")

    else:
        raise ValueError("resample must be 'template' or 'observed'")


    # Perform the cross-correlation
    padding = observed_flux.size + template_flux.size
    # Is this necessary?: # TODO
    x_norm = observed_flux - np.mean(observed_flux[np.isfinite(observed_flux)])
    y_norm = template_flux - np.mean(template_flux[np.isfinite(template_flux)])

    Fx = np.fft.fft(x_norm, padding, )
    Fy = np.fft.fft(y_norm, padding, )
    iFxy = np.fft.ifft(Fx.conj() * Fy).real
    varxy = np.sqrt(np.inner(x_norm, x_norm) * np.inner(y_norm, y_norm))

    fft_result = iFxy/varxy

    # Put around symmetry axis.
    num = len(fft_result) - 1 if len(fft_result) % 2 else len(fft_result)

    fft_y = np.zeros(num)
    fft_y[:num/2] = fft_result[num/2:num]
    fft_y[num/2:] = fft_result[:num/2]

    fft_x = np.arange(num) - num/2

    # Get initial guess of peak.
    p0 = np.array([fft_x[np.argmax(fft_y)], np.max(fft_y), 10])

    gaussian = lambda p, x: p[1] * np.exp(-(x - p[0])**2 / (2.0 * p[2]**2))
    errfunc = lambda p, x, y: y - gaussian(p, x)

    try:
        p1, ier = leastsq(errfunc, p0.copy(), args=(fft_x, fft_y))

    except:
        logger.exception("Exception in measuring peak of CCF:")
        raise


    # Create functions for interpolating back onto the dispersion map
    fft_points = (0, p1[0])
    interp_x = np.arange(num/2) - num/4

    wl_points = []
    for point in fft_points:
        idx = np.searchsorted(interp_x, point)
        f = interpolate.interp1d(interp_x[idx-3:idx+3], dispersion[idx-3:idx+3],
            bounds_error=True, kind='cubic')
        wl_points.append(f(point))

    # Calculate velocity 
    c = 299792458e-3 # km/s
    f, g = wl_points
    rv = c * (1 - g/f)

    # Uncertainty
    rv_uncertainty = np.nan # TODO

    # Create a CCF spectrum.
    ccf = np.array([fft_x * (rv/p1[0]), fft_y])

    return (rv, rv_uncertainty, ccf)



