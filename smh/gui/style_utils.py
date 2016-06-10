#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["wavelength_to_hex", "relim_axes", "fill_between_steps"]

import numpy as np
from matplotlib.collections import PathCollection

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

    red = int(np.clip(round(intensity * (red * factor)**gamma), 0, intensity))
    green = int(np.clip(round(intensity * (green * factor)**gamma), 0, intensity))
    blue = int(np.clip(round(intensity * (blue * factor)**gamma), 0, intensity))

    return "#{0:02x}{1:02x}{2:02x}".format(red, green, blue)


def fill_between_steps(ax, x, y1, y2=0, h_align='mid', **kwargs):
    """
    Fill between for step plots in matplotlib.

    **kwargs will be passed to the matplotlib fill_between() function.
    """

    # If no Axes opject given, grab the current one:

    # First, duplicate the x values
    xx = x.repeat(2)[1:]
    # Now: the average x binwidth
    xstep = np.repeat((x[1:] - x[:-1]), 2)
    xstep = np.concatenate(([xstep[0]], xstep, [xstep[-1]]))
    # Now: add one step at end of row.
    xx = np.append(xx, xx.max() + xstep[-1])

    # Make it possible to chenge step alignment.
    if h_align == 'mid':
        xx -= xstep / 2.
    elif h_align == 'right':
        xx -= xstep

    # Also, duplicate each y coordinate in both arrays
    y1 = y1.repeat(2)#[:-1]
    if type(y2) == np.ndarray:
        y2 = y2.repeat(2)#[:-1]

    # now to the plotting part:
    return ax.fill_between(xx, y1, y2=y2, **kwargs)




def relim_axes(axes, percent=20):
    """
    Generate new axes for a matplotlib axes based on the collections present.

    :param axes:
        The matplotlib axes.

    :param percent: [optional]
        The percent of the data to extend past the minimum and maximum data
        points.

    :returns:
        A two-length tuple containing the lower and upper limits in the x- and
        y-axis, respectively.
    """

    data = np.vstack([item.get_offsets() for item in axes.collections \
        if isinstance(item, PathCollection)])

    if data.size == 0:
        return (None, None)
    
    data = data.reshape(-1, 2)
    x, y = data[:,0], data[:, 1]

    # Only use finite values.
    finite = np.isfinite(x*y)
    x, y = x[finite], y[finite]

    if x.size > 1:
        xlim = [
            np.min(x) - np.ptp(x) * percent/100.,
            np.max(x) + np.ptp(x) * percent/100.,
        ]

    elif x.size == 0:
        xlim = None

    else:
        xlim = (x[0] - 1, x[0] + 1)


    if y.size > 1:
        ylim = [
            np.min(y) - np.ptp(y) * percent/100.,
            np.max(y) + np.ptp(y) * percent/100.
        ]

    elif y.size == 0: 
        ylim = None

    else:
        ylim = (y[0] - 1, y[0] + 1)

    axes.set_xlim(xlim)
    axes.set_ylim(ylim)

    return (xlim, ylim)

