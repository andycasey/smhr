#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["wavelength_to_hex"]

from numpy import clip


def wavelength_to_hex(wavelength):
    """
    Convert a wavelength into its approximate hexadecimal color code.

    :param wavelength:
        Wavelength value in Angstroms.
    """

    gamma, intensity = 0.8, 255

    if wavelength >= 3800 and wavelength < 4400:
        red = -(wavelength - 4400) / (4400 - 3800)
        green, blue = 0, 1

    elif wavelength >= 4400 and wavelength < 4900:
        green = (wavelength - 4400) / (4900 - 4400)
        red, blue = 0, 1

    elif wavelength >= 4900 and wavelength < 5100:
        blue = -(wavelength - 5100) / (5100 - 4900)
        red, green = 0, 1

    elif wavelength >= 5100 and wavelength < 5800:
        red = (wavelength - 5100) / (5800 - 5100)
        green, blue = 1, 0

    elif wavelength >= 5800 and wavelength < 6450:
        green = -(wavelength - 6450) / (6450 - 5800)
        red, blue = 1, 0

    elif wavelength >= 6450 and wavelength < 7810:
        red, green, blue = 1, 0, 0

    else:
        red, green, blue = 0, 0, 0

    if wavelength >= 3800 and wavelength < 4200:
        factor = 0.3 + 0.7 * (wavelength - 3800)/(4200 - 3800)

    elif wavelength >= 4200 and wavelength < 7010:
        factor = 1.0

    elif wavelength >= 7010 and wavelength < 7810:
        factor = 0.3 + 0.7 * (7800 - wavelength) / (7800 - 7000)

    else:
        factor = 0

    red = int(clip(round(intensity * (red * factor)**gamma), 0, intensity))
    green = int(clip(round(intensity * (green * factor)**gamma), 0, intensity))
    blue = int(clip(round(intensity * (blue * factor)**gamma), 0, intensity))

    return "#{0:02x}{1:02x}{2:02x}".format(red, green, blue)