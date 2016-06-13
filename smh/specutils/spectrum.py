#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" An object for dealing with one-dimensional spectra. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Spectrum1D", "stitch"]

import logging
import numpy as np
import os
import re
from collections import OrderedDict
from hashlib import md5

from astropy.io import fits
from scipy import interpolate, ndimage, polyfit, poly1d, optimize as op

logger = logging.getLogger(__name__)

class Spectrum1D(object):
    """ A one-dimensional spectrum. """

    def __init__(self, dispersion, flux, ivar, metadata=None):
        """
        Initialie a `Spectrum1D` object with the given dispersion, flux and
        inverse variance arrays.

        :param dispersion:
            An array containing the dispersion values for every pixel.

        :param flux:
            An array containing flux values for all pixels.

        :param ivar:
            An array containing the inverse variances for the flux values in all
            pixels.

        :param metadata: [optional]
            A dictionary containing metadata for this spectrum.
        """

        dispersion = np.array(dispersion)
        flux = np.array(flux)
        ivar = np.array(ivar)

        if max([len(_.shape) for _ in (dispersion, flux, ivar)]) > 1:
            raise ValueError(
                "dispersion, flux and ivar must be one dimensional arrays")
        
        if flux.size != dispersion.size:
            raise ValueError(
                "dispersion and flux arrays must have the same "
                "size ({0} != {1})".format(dispersion.size, flux.size, ))

        if ivar.size != dispersion.size:
            raise ValueError(
                "dispersion and ivar arrays must have the same "
                "size ({0} != {1})".format(dispersion.size, ivar.size, ))

        self.metadata = metadata or {}
        self._dispersion = dispersion
        self._flux = flux
        self._ivar = ivar

        return None
    

    @property
    def dispersion(self):
        """ Return the dispersion points for all pixels. """
        return self._dispersion


    @property
    def flux(self):
        """ Return the flux for all pixels. """
        return self._flux


    @property
    def ivar(self):
        """ Return the inverse variance of the flux for all pixels. """
        return self._ivar


    @classmethod
    def read(cls, path, **kwargs):
        """
        Create a Spectrum1D class from a path on disk.

        :param path:
            The path of the filename to load.
        """

        if not os.path.exists(path):
            raise IOError("filename '{}' does not exist".format(path))

        # Try multi-spec first since this is currently the most common use case.
        methods = (
            cls.read_fits_multispec,
            cls.read_fits_spectrum1d,
            cls.read_ascii_spectrum1d
        )

        for method in methods:
            try:
                dispersion, flux, ivar, metadata = method(path, **kwargs)

            except:
                continue

            else:
                orders = [cls(dispersion=d, flux=f, ivar=i, metadata=metadata) \
                    for d, f, i in zip(dispersion, flux, ivar)]
                break
        else:
            raise ValueError("cannot read spectrum from path {}".format(path))

        # If it's a single order, just return that instead of a 1-length list.
        orders = orders if len(orders) > 1 else orders[0]
        return orders


    @classmethod
    def read_fits_multispec(cls, path, flux_ext=None, ivar_ext=None, **kwargs):
        """
        Create multiple Spectrum1D classes from a multi-spec file on disk.

        :param path:
            The path of the multi-spec filename to load.

        :param flux_ext: [optional]
            The zero-indexed extension number containing flux values.

        :param ivar_ext: [optional]
            The zero-indexed extension number containing the inverse variance of
            the flux values.
        """

        image = fits.open(path)

        # Merge headers into a metadata dictionary.
        metadata = OrderedDict()
        for key, value in image[0].header.items():

            # NOTE: In the old SMH we did a try-except block to string-ify and
            #       JSON-dump the header values, and if they could not be
            #       forced to a string we didn't keep that header.

            #       I can't remember what types caused that problem, but it was
            #       to prevent SMH being unable to save a session.
            
            #       Since we are pickling now, that shouldn't be a problem
            #       anymore, but this note is here to speed up debugging in case
            #       that issue returns.

            if key in metadata:
                metadata[key] += value
            else:
                metadata[key] = value

        metadata["smh_read_path"] = path

        flux = image[0].data

        assert metadata["CTYPE1"].upper().startswith("MULTISPE") \
            or metadata["WAT0_001"].lower() == "system=multispec"

        # Join the WAT keywords for dispersion mapping.
        i, concatenated_wat, key_fmt = (1, str(""), "WAT2_{0:03d}")
        while key_fmt.format(i) in metadata:
            # .ljust(68, " ") had str/unicode issues across Python 2/3
            value = metadata[key_fmt.format(i)]
            concatenated_wat += value + (" "  * (68 - len(value)))
            i += 1

        # Split the concatenated header into individual orders.
        order_mapping = np.array([map(float, each.rstrip('" ').split()) \
                for each in re.split('spec[0-9]+ ?= ?"', concatenated_wat)[1:]])

        # Parse the order mapping into dispersion values.
        dispersion = np.array(
            [compute_dispersion(*mapping) for mapping in order_mapping])

        # Get the flux and inverse variance arrays.
        # NOTE: Most multi-spec data previously used with SMH have been from
        #       Magellan/MIKE, and reduced with CarPy.
        md5_hash = md5(";".join([v for k, v in metadata.items() \
            if k.startswith("BANDID")])).hexdigest()
        is_carpy_product = (md5_hash == "0da149208a3c8ba608226544605ed600")

        if is_carpy_product:
            # CarPy gives a 'noise' spectrum, which we must convert to an
            # inverse variance array
            flux_ext = flux_ext or 1
            noise_ext = ivar_ext or 2

            flux = image[0].data[flux_ext]
            ivar = image[0].data[noise_ext]**(-2)

        else:
            ivar = np.nan * np.ones_like(flux)
            #raise ValueError("could not identify flux and ivar extensions")

        dispersion = np.atleast_2d(dispersion)
        flux = np.atleast_2d(flux)
        ivar = np.atleast_2d(ivar)

        # Ensure dispersion maps from blue to red direction.
        if np.min(dispersion[0]) > np.min(dispersion[-1]):

            dispersion = dispersion[::-1]
            if len(flux.shape) > 2:
                flux = flux[:, ::-1]
                ivar = ivar[:, ::-1]
            else:
                flux = flux[::-1]
                ivar = ivar[::-1]

        # Do something sensible regarding zero or negative fluxes.
        ivar[0 >= flux] = 0
        flux[0 >= flux] = np.nan

        return (dispersion, flux, ivar, metadata)


    @classmethod
    def read_fits_spectrum1d(cls, path, **kwargs):
        """
        Read Spectrum1D data from a binary FITS file.

        :param path:
            The path of the FITS filename to read.
        """

        image = fits.open(path)

        # Merge headers into a metadata dictionary.
        metadata = OrderedDict()
        for key, value in image[0].header.items():
            if key in metadata:
                metadata[key] += value
            else:
                metadata[key] = value
        metadata["smh_read_path"] = path

        # Find the first HDU with data in it.
        for hdu_index, hdu in enumerate(image):
            if hdu.data is not None: break

        if len(image) == 2 and hdu_index == 1:            
            dispersion = image[hdu_index].data["dispersion"]
            flux = image[hdu_index].data["flux"]
            ivar = image[hdu_index].data["ivar"]

        else:
            # Build a simple linear dispersion map from the headers.
            # See http://iraf.net/irafdocs/specwcs.php
            crval = image[0].header["CRVAL1"]
            naxis = image[0].header["NAXIS1"]
            crpix = image[0].header.get("CRPIX1", 0)
            cdelt = image[0].header["CDELT1"]
            ltv = image[0].header.get("LTV1", 0)

            dispersion = \
                crval + (np.arange(naxis) - crpix) * cdelt - ltv * cdelt

            flux = image[0].data
            if len(image) == 1:
                ivar = np.ones_like(flux)*1e+5 # HACK S/N ~300 just for training/verification purposes
            else:
                ivar = image[1].data

        dispersion = np.atleast_2d(dispersion)
        flux = np.atleast_2d(flux)
        ivar = np.atleast_2d(ivar)

        return (dispersion, flux, ivar, metadata)


    @classmethod
    def read_ascii_spectrum1d(cls, path, **kwargs):
        """
        Read Spectrum1D data from an ASCII-formatted file on disk.

        :param path:
            The path of the ASCII filename to read.
        """

        kwds = kwargs.copy()
        kwds.update({
            "unpack": True
        })
        kwds.setdefault("usecols", (0, 1, 2))

        try:
            dispersion, flux, ivar = np.loadtxt(path, **kwds)
        except:
            # Try by ignoring the first row.
            kwds.setdefault("skiprows", 1)
            dispersion, flux, ivar = np.loadtxt(path, **kwds)

        dispersion = np.atleast_2d(dispersion)
        flux = np.atleast_2d(flux)
        ivar = np.atleast_2d(ivar)
        metadata = { "smh_read_path": path }
        
        return (dispersion, flux, ivar, metadata)


    def write(self, filename):
        """ Write spectrum to disk. """

        # TODO: only ascii atm.
        a = np.array([self.dispersion, self.flux, self.ivar]).T
        np.savetxt(filename, a)


    def redshift(self, v=None, z=None):
        """
        Redshift the spectrum.

        :param v:
            The velocity in km/s.

        :param z:
            A redshift.
        """

        if (v is None and z is None) or (v is not None and z is not None):
            raise ValueError("either v or z must be given, but not both")

        z = z or v/299792458e-3
        self._dispersion *= 1 + z
        return True


    # State functionality for serialization.
    def __getstate__(self):
        """ Return the spectrum state. """
        return (self.dispersion, self.flux, self.ivar, self.metadata)


    def __setstate__(self, state):
        """
        Set the state of the spectrum.

        :param state:
            A four-length tuple containing the dispersion array, flux array, the
            inverse variance of the fluxes, and a metadata dictionary.
        """
        
        dispersion, flux, ivar, metadata = state
        self._dispersion = dispersion
        self._flux = flux
        self._ivar = ivar
        self.metadata = metadata
        return None


    def copy(self):
        """
        Create a copy of the spectrum.
        """
        return self.__class__(
            dispersion=self.dispersion.copy(),
            flux=self.flux.copy(), ivar=self.ivar.copy(),
            metadata=self.metadata.copy())


    def gaussian_smooth(self, fwhm, **kwargs):
        
        profile_sigma = fwhm / (2 * (2*np.log(2))**0.5)
        
        # The requested FWHM is in Angstroms, but the dispersion between each
        # pixel is likely less than an Angstrom, so we must calculate the true
        # smoothing value
        
        true_profile_sigma = profile_sigma / np.median(np.diff(self.disp))
        smoothed_flux = ndimage.gaussian_filter1d(self.flux, true_profile_sigma, **kwargs)
        
        return self.__class__(self.disp, smoothed_flux)
        
    
    def smooth(self, window_len=11, window='hanning'):
        """Smooth the data using a window with requested size.
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        input:
            window_len: the dimension of the smoothing window; should be an odd integer
            window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.

        output:
            the smoothed spectra
            
        example:

        (This function from Heather Jacobson)
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        """

        if self.flux.size < window_len:
            raise ValueError("input vector must be bigger than the window size")

        if window_len<3:
            return self

        available = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
        if window not in available:
            raise ValueError("window is one of {}".format(", ".join(available)))


        s = np.r_[self.flux[window_len-1:0:-1], self.flux, self.flux[-1:-window_len:-1]]
        
        if window == 'flat': #moving average
            w = np.ones(window_len, 'd')

        else:
            w = eval('np.' + window + '(window_len)', {'__builtins__': None})

        smoothed_flux = np.convolve(w/w.sum(), s,mode='valid')
        smoothed_flux = smoothed_flux[(window_len/2-1):-(window_len/2)]

        return self.__class__(self.disp, smoothed_flux, headers=self.headers)
    


    def fit_continuum(self, knot_spacing=200, low_sigma_clip=1.0, \
        high_sigma_clip=0.2, max_iterations=3, order=3, exclude=None, \
        include=None, additional_points=None, function='spline', scale=1.0,
        **kwargs):
        """
        Fits the continuum for a given `Spectrum1D` spectrum.
        
        Parameters
        ----
        knot_spacing : float or None, optional
            The knot spacing for the continuum spline function in Angstroms. Optional.
            If not provided then the knot spacing will be determined automatically.
        
        sigma_clip : a tuple of two floats, optional
            This is the lower and upper sigma clipping level respectively. Optional.
            
        max_iterations : int, optional
            Maximum number of spline-fitting operations.
            
        order : int, optional
            The order of the spline function to fit.
            
        exclude : list of tuple-types containing floats, optional
            A list of wavelength regions to always exclude when determining the
            continuum. Example:
            
            >> exclude = [
            >>    (3890.0, 4110.0),
            >>    (4310.0, 4340.0)
            >>  ]
            
            In the example above the regions between 3890 A and 4110 A, as well as
            4310 A to 4340 A will always be excluded when determining the continuum
            regions.

        function: only 'spline' or 'poly'

        scale : float
            A scaling factor to apply to the normalised flux levels.
            
        include : list of tuple-types containing floats, optional
            A list of wavelength regions to always include when determining the
            continuum.

        :param rv:
            A radial velocity correction (in km/s) to apply to the spectrum.
        """

        dispersion = self.dispersion.copy()

        exclusions = []
        continuum_indices = range(len(self.flux))

        # Snip left and right
        finite_positive_flux = np.isfinite(self.flux) * self.flux > 0

        function = str(function).lower()
        left = np.where(finite_positive_flux)[0][0]
        right = np.where(finite_positive_flux)[0][-1]

        # See if there are any regions we need to exclude
        if exclude is not None and len(exclude) > 0:
            exclude_indices = []
            
            if isinstance(exclude[0], float) and len(exclude) == 2:
                # Only two floats given, so we only have one region to exclude
                exclude_indices.extend(range(*np.searchsorted(dispersion, exclude)))
                
            else:
                # Multiple regions provided
                for exclude_region in exclude:
                    exclude_indices.extend(
                        range(*np.searchsorted(dispersion, exclude_region)))
        
            continuum_indices = np.sort(list(set(continuum_indices).difference(
                np.sort(exclude_indices))))
            
        # See if there are any regions we should always include
        if include is not None and len(include) > 0:
            include_indices = []
            
            if isinstance(include[0], float) and len(include) == 2:
                # Only two floats given, so we can only have one region to include
                include_indices.extend(range(*np.searchsorted(dispersion, include)))
                
            else:
                # Multiple regions provided
                for include_region in include:
                    include_indices.extend(range(*np.searchsorted(
                        dispersion, include_region)))
        
        # We should exclude non-finite values from the fit.
        non_finite_indices = np.where(~np.isfinite(self.flux * self.ivar))[0]
        continuum_indices = np.sort(list(set(continuum_indices).difference(
            non_finite_indices)))

        # We should also exclude zero or negative flux points from the fit
        zero_flux_indices = np.where(0 >= self.flux)[0]
        continuum_indices = np.sort(list(set(continuum_indices).difference(
            zero_flux_indices)))

        if 1 > continuum_indices.size:
            no_continuum = np.nan * np.ones_like(dispersion)
            failed_spectrum = self.__class__(dispersion=dispersion,
                flux=no_continuum, ivar=no_continuum, metadata=self.metadata)

            if kwargs.get("full_output", False):
                return (failed_spectrum, no_continuum, 0, dispersion.size - 1)

            return failed_spectrum        

        original_continuum_indices = continuum_indices.copy()

        if knot_spacing is None or knot_spacing == 0:
            knots = []

        else:
            knot_spacing = abs(knot_spacing)
            
            end_spacing = ((dispersion[-1] - dispersion[0]) % knot_spacing) /2.
            if knot_spacing/2. > end_spacing: end_spacing += knot_spacing/2.
                
            knots = np.arange(dispersion[0] + end_spacing, 
                dispersion[-1] - end_spacing + knot_spacing, 
                knot_spacing)

            if len(knots) > 0 and knots[-1] > dispersion[continuum_indices][-1]:
                knots = knots[:knots.searchsorted(dispersion[continuum_indices][-1])]
                
            if len(knots) > 0 and knots[0] < dispersion[continuum_indices][0]:
                knots = knots[knots.searchsorted(dispersion[continuum_indices][0]):]

        # TODO: Use inverse variance array when fitting polynomial/spline.
        for iteration in range(max_iterations):
            
            if 1 > continuum_indices.size:

                no_continuum = np.nan * np.ones_like(dispersion)
                failed_spectrum = self.__class__(dispersion=dispersion,
                    flux=no_continuum, ivar=no_continuum, metadata=self.metadata)

                if kwargs.get("full_output", False):
                    return (failed_spectrum, no_continuum, 0, dispersion.size - 1)

                return failed_spectrum

            splrep_disp = dispersion[continuum_indices]
            splrep_flux = self.flux[continuum_indices]
            splrep_weights = self.ivar[continuum_indices]**0.5

            median_weight = np.nanmedian(splrep_weights)

            # We need to add in additional points at the last minute here
            if additional_points is not None and len(additional_points) > 0:

                for point, flux, weight in additional_points:

                    # Get the index of the fit
                    insert_index = int(np.searchsorted(splrep_disp, point))
                    
                    # Insert the values
                    splrep_disp = np.insert(splrep_disp, insert_index, point)
                    splrep_flux = np.insert(splrep_flux, insert_index, flux)
                    splrep_weights = np.insert(splrep_weights, insert_index, 
                        median_weight * weight)

            if function == 'spline':
                if order > 5:
                    logger.warn("Splines can only have a maximum order of 5. "
                        "Limiting order value to 5.")
                    order = 5

                try:
                    tck = interpolate.splrep(splrep_disp, splrep_flux,
                        k=order, task=-1, t=knots, w=splrep_weights)

                except:
                    logger.exception("Exception in fitting continuum:")
                    continuum = np.nan * np.ones_like(dispersion)

                else:
                    continuum = interpolate.splev(dispersion, tck)

            elif function in ("poly", "polynomial"):

                coeffs = polyfit(splrep_disp, splrep_flux, order)

                popt, pcov = op.curve_fit(lambda x, *c: np.polyval(c, x), 
                    splrep_disp, splrep_flux, coeffs, 
                    sigma=self.ivar[continuum_indices], absolute_sigma=False)
                continuum = np.polyval(popt, dispersion)

            else:
                raise ValueError("Unknown function type: only spline or poly "\
                    "available ({} given)".format(function))
            
            difference = continuum - self.flux
            sigma_difference = difference \
                / np.std(difference[np.isfinite(self.flux)])

            # Clipping
            upper_exclude = np.where(sigma_difference > high_sigma_clip)[0]
            lower_exclude = np.where(sigma_difference < -low_sigma_clip)[0]
            
            exclude_indices = list(upper_exclude)
            exclude_indices.extend(lower_exclude)
            exclude_indices = np.array(exclude_indices)
            
            if len(exclude_indices) is 0: break
            
            exclusions.extend(exclude_indices)
            
            # Before excluding anything, we must check to see if there are regions
            # which we should never exclude
            if include is not None:
                exclude_indices \
                    = set(exclude_indices).difference(include_indices)
            
            # Remove regions that have been excluded
            continuum_indices = np.sort(list(set(continuum_indices).difference(
                exclude_indices)))
        
        # Snip the edges based on exclude regions
        if exclude is not None and len(exclude) > 0:

            # If there are exclusion regions that extend past the left/right,
            # then we will need to adjust left/right accordingly

            left = np.max([left, np.min(original_continuum_indices)])
            right = np.min([right, np.max(original_continuum_indices)])
            
        # Apply flux scaling
        continuum *= scale

        normalized_spectrum = self.__class__(
            dispersion=dispersion[left:right],
            flux=(self.flux/continuum)[left:right], 
            ivar=(continuum * self.ivar * continuum)[left:right],
            metadata=self.metadata)

        # Return a normalized spectrum.
        continuum[:left] = np.nan
        continuum[right:] = np.nan
        if kwargs.get("full_output", False):
            return (normalized_spectrum, continuum, left, right)
        return normalized_spectrum


def compute_dispersion(aperture, beam, dispersion_type, dispersion_start,
    mean_dispersion_delta, num_pixels, redshift, aperture_low, aperture_high,
    weight=1, offset=0, function_type=None, order=None, Pmin=None, Pmax=None,
    *coefficients):
    """
    Compute a dispersion mapping from a IRAF multi-spec description.

    :param aperture:
        The aperture number.

    :param beam:
        The beam number.

    :param dispersion_type:
        An integer representing the dispersion type:

        0: linear dispersion
        1: log-linear dispersion
        2: non-linear dispersion

    :param dispersion_start:
        The value of the dispersion at the first physical pixel.

    :param mean_dispersion_delta:
        The mean difference between dispersion pixels.

    :param num_pixels:
        The number of pixels.

    :param redshift:
        The redshift of the object. This is accounted for by adjusting the
        dispersion scale without rebinning:

        >> dispersion_adjusted = dispersion / (1 + redshift)

    :param aperture_low:
        The lower limit of the spatial axis used to compute the dispersion.

    :param aperture_high:
        The upper limit of the spatial axis used to compute the dispersion.

    :param weight: [optional]
        A multiplier to apply to all dispersion values.

    :param offset: [optional]
        A zero-point offset to be applied to all the dispersion values.

    :param function_type: [optional]
        An integer representing the function type to use when a non-linear 
        dispersion mapping (i.e. `dispersion_type = 2`) has been specified:

        1: Chebyshev polynomial
        2: Legendre polynomial
        3: Cubic spline
        4: Linear spline
        5: Pixel coordinate array
        6: Sampled coordinate array

    :param order: [optional]
        The order of the Legendre or Chebyshev function supplied.

    :param Pmin: [optional]
        The minimum pixel value, or lower limit of the range of physical pixel
        coordinates.

    :param Pmax: [optional]
        The maximum pixel value, or upper limit of the range of physical pixel
        coordinates.

    :param coefficients: [optional]
        The `order` number of coefficients that define the Legendre or Chebyshev
        polynomial functions.

    :returns:
        An array containing the computed dispersion values.
    """

    if dispersion_type in (0, 1):
        # Simple linear or logarithmic spacing
        dispersion = \
            dispersion_start + np.arange(num_pixels) * mean_dispersion_delta

        if dispersion_start == 1:
            dispersion = 10.**dispersion

    elif dispersion_type == 2:
        # Non-linear mapping.
        if function_type is None:
            raise ValueError("function type required for non-linear mapping")
        elif function_type not in range(1, 7):
            raise ValueError(
                "function type {0} not recognised".format(function_type))

        if function_type == 1:
            order = int(order)
            n = np.linspace(-1, 1, Pmax - Pmin + 1)
            temp = np.zeros((Pmax - Pmin + 1, order), dtype=float)
            temp[:, 0] = 1
            temp[:, 1] = n
            for i in range(2, order):
                temp[:, i] = 2 * n * temp[:, i-1] - temp[:, i-2]
            
            for i in range(0, order):
                temp[:, i] *= coefficients[i]

            dispersion = temp.sum(axis=1)


        elif function_type == 2:
            # Legendre polynomial.
            if None in (order, Pmin, Pmax, coefficients):
                raise TypeError("order, Pmin, Pmax and coefficients required "
                                "for a Chebyshev or Legendre polynomial")

            Pmean = (Pmax + Pmin)/2
            Pptp = Pmax - Pmin
            x = (np.arange(num_pixels) + 1 - Pmean)/(Pptp/2)
            p0 = np.ones(num_pixels)
            p1 = mean_dispersion_delta

            dispersion = coefficients[0] * p0 + coefficients[1] * p1
            for i in range(2, int(order)):
                if function_type == 1:
                    # Chebyshev
                    p2 = 2 * x * p1 - p0
                else:
                    # Legendre
                    p2 = ((2*i - 1)*x*p1 - (i - 1)*p0) / i

                dispersion += p2 * coefficients[i]
                p0, p1 = (p1, p2)

        elif function_type == 3:
            # Cubic spline.
            if None in (order, Pmin, Pmax, coefficients):
                raise TypeError("order, Pmin, Pmax and coefficients required "
                                "for a cubic spline mapping")
            s = (np.arange(num_pixels, dtype=float) + 1 - Pmin)/(Pmax - Pmin) \
              * order
            j = s.astype(int).clip(0, order - 1)
            a, b = (j + 1 - s, s - j)
            x = np.array([
                a**3,
                1 + 3*a*(1 + a*b),
                1 + 3*b*(1 + a*b),
                b**3])
            dispersion = np.dot(np.array(coefficients), x.T)

        else:
            raise NotImplementedError("function type not implemented yet")

    else:
        raise ValueError(
            "dispersion type {0} not recognised".format(dispersion_type))

    # Apply redshift correction.
    dispersion = weight * (dispersion + offset) / (1 + redshift)
    return dispersion



def common_dispersion_map(spectra, full_output=True):
    """
    Produce a common dispersion mapping for (potentially overlapping) spectra
    and minimize the number of resamples required. Pixel bin edges are
    preferred from bluer to redder wavelengths.

    :param spectra:
        A list of spectra to produce the mapping for.

    :param full_output: [optional]
        Optinally return the indices, and the sorted spectra.

    :returns:
        An array of common dispersion values. If `full_output` is set to `True`,
        then a three-length tuple will be returned, containing the common
        dispersion pixels, the common dispersion map indices for each spectrum,
        and a list of the sorted spectra.
    """

    # Make spectra blue to right.
    spectra = sorted(spectra, key=lambda s: s.dispersion[0])

    common = []
    discard_bluest_pixels = None
    for i, blue_spectrum in enumerate(spectra[:-1]):
        red_spectrum = spectra[i + 1]

        # Do they overlap?
        if blue_spectrum.dispersion[-1] >= red_spectrum.dispersion[0] \
        and red_spectrum.dispersion[-1] >= blue_spectrum.dispersion[0]:
            
            # Take the "entire" blue spectrum then discard some blue pixels from
            # the red spectrum.
            common.extend(blue_spectrum.dispersion[discard_bluest_pixels:])
            discard_bluest_pixels = red_spectrum.dispersion.searchsorted(
                blue_spectrum.dispersion[-1])

        else:
            # Can just extend the existing map, modulo the first N-ish pixels.
            common.extend(blue_spectrum.dispersion[discard_bluest_pixels:])
            discard_bluest_pixels = None

    # For the last spectrum.
    if len(spectra) > 1:
        common.extend(red_spectrum.dispersion[discard_bluest_pixels:])

    else:
        common = spectra[0].dispersion.copy()

    # Ensure that we are contiguous from blue to red.
    common = np.array(common)
    common = common[np.argsort(common)]

    assert np.all(np.diff(common) > 0), \
        "Spectra must be contiguous from blue to red"

    if full_output:
        indices = [common.searchsorted(s.dispersion) for s in spectra]
        return (common, indices, spectra)

    return common


def stitch(spectra):
    """
    Stitch spectra together, some of which may have overlapping dispersion
    ranges. This is a crude (knowingly incorrect) approximation: we interpolate
    fluxes without ensuring total conservation.

    :param spectra:
        A list of potentially overlapping spectra.
    """

    # Create common mapping.
    N = len(spectra)
    dispersion, indices, spectra = common_dispersion_map(spectra, True)
    common_flux = np.zeros((N, dispersion.size))
    common_ivar = np.zeros_like(common_flux)

    for i, (j, spectrum) in enumerate(zip(indices, spectra)):
        common_flux[i, j] = np.interp(
            dispersion[j], spectrum.dispersion, spectrum.flux,
            left=0, right=0)
        common_ivar[i, j] = spectrum.ivar

    finite = np.isfinite(common_flux * common_ivar)
    common_flux[~finite] = 0
    common_ivar[~finite] = 0

    numerator = np.sum(common_flux * common_ivar, axis=0)
    denominator = np.sum(common_ivar, axis=0)

    flux, ivar = (numerator/denominator, denominator)

    # Create a spectrum with no header provenance.
    return Spectrum1D(dispersion, flux, ivar)
    

