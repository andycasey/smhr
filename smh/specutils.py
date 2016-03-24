#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Spectroscopy-related utilities. """


from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Spectrum", "Spectrum1D"]

import os
import json
import logging
import re
from collections import OrderedDict
from hashlib import md5
from shutil import copyfile

import numpy as np
from astropy.io import fits

from scipy import interpolate, ndimage, polyfit, poly1d
from scipy.optimize import leastsq

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
            cls._read_multispec_fits,
            cls._read_spectrum1d_fits,
            cls._read_spectrum1d_ascii
        )

        for method in methods:
            try:
                dispersion, flux, ivar, metadata = method(path, **kwargs)

            except:
                logger.exception("Exception in trying to load {0} with {1}:"\
                    .format(path, method))
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
    def _read_multispec_fits(cls, path, flux_ext=None, ivar_ext=None, **kwargs):
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
        N_pixels = flux.shape[-1]
        N_orders = 1 if len(flux.shape) == 1 else flux.shape[-2]

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
            raise ValueError("could not identify flux and ivar extensions")

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

        return (dispersion, flux, ivar, metadata)


    @classmethod
    def _read_spectrum1d_fits(cls, path, **kwargs):
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
                # HACK TODO: Issue 4
                ivar = np.nan * np.ones_like(flux)
            else:
                ivar = image[1].data

        dispersion = np.atleast_2d(dispersion)
        flux = np.atleast_2d(flux)
        ivar = np.atleast_2d(ivar)

        return (dispersion, flux, ivar, metadata)


    @classmethod
    def _read_spectrum1d_ascii(cls, path, **kwargs):
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


    def gaussian_smooth(self, fwhm, **kwargs):
        
        profile_sigma = fwhm / (2 * (2*np.log(2))**0.5)
        
        # The requested FWHM is in Angstroms, but the dispersion between each
        # pixel is likely less than an Angstrom, so we must calculate the true
        # smoothing value
        
        true_profile_sigma = profile_sigma / np.median(np.diff(self.disp))
        print("true profile sigma", profile_sigma, np.diff(self.disp))
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
            raise ValueError, "Input vector needs to be bigger than window size."

        if window_len<3:
            return self

        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


        s = numpy.r_[self.flux[window_len-1:0:-1], self.flux, self.flux[-1:-window_len:-1]]
        
        if window == 'flat': #moving average
            w = numpy.ones(window_len, 'd')

        else:
            w = eval('numpy.' + window + '(window_len)', {'__builtins__': None})

        smoothed_flux = numpy.convolve(w/w.sum(), s,mode='valid')
        smoothed_flux = smoothed_flux[(window_len/2-1):-(window_len/2)]

        return self.__class__(self.disp, smoothed_flux, headers=self.headers)
    


    def doppler_correct(self, v=None, z=None, interpolate=True):
        """Performs a Doppler correction on the given `Spectrum1D` object by the
        amount required. Either velocities or redshifts can be provided.
        
        Parameters
        ----
        v : float
            The velocity (in km/s) to correct the `Spectrum1D` object by.
            
        z : float
            The redshift to correct the `Spectrum1D` object by.
            
        Raises
        ----
        ValueError
            When both a velocity `v` and a redshift `z` are provided for Doppler
            corrections.
        """
        
        if isinstance(v, float) and isinstance(z, float):
            raise ValueError('Both a velocity and redshift was supplied')
            
        c = 299792458e-3 # km/s

        if v is not None:
            new_disp = self.disp * np.sqrt((1+v/c)/(1-v/c))
            
        elif z is not None:
            new_disp = self.disp * (1 + z)
            
        else:
            raise TypeError('No velocity or redshift supplied.')
        
        rest_spectrum = self.__class__(new_disp, self.flux)
          
        if interpolate:
            return rest_spectrum.interpolate(self.disp)
            
        else:
            return rest_spectrum
        
    
    def fit_continuum(self, knot_spacing=200, sigma_clip=(1.0, 0.2), \
                      max_iterations=3, order=3, exclude=None, include=None, \
                      additional_points=None, function='spline', scale=1.0, **kwargs):
        """Fits the continuum for a given `Spectrum1D` spectrum.
        
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
        """
        
        exclusions = []
        continuum_indices = range(len(self.flux))

        # Snip left and right
        finite_positive_flux = np.isfinite(self.flux) * self.flux > 0

        #print "finite flux", np.any(finite_positive_flux), finite_positive_flux
        #print "where flux", np.where(finite_positive_flux)
        #print "flux is...", self.flux
        left_index = np.where(finite_positive_flux)[0][0]
        right_index = np.where(finite_positive_flux)[0][-1]

        # See if there are any regions we need to exclude
        if exclude is not None and len(exclude) > 0:
            exclude_indices = []
            
            if isinstance(exclude[0], float) and len(exclude) == 2:
                # Only two floats given, so we only have one region to exclude
                exclude_indices.extend(range(*np.searchsorted(self.dispersion, exclude)))
                
            else:
                # Multiple regions provided
                for exclude_region in exclude:
                    exclude_indices.extend(range(*np.searchsorted(self.dispersion, exclude_region)))
        
            continuum_indices = np.sort(list(set(continuum_indices).difference(np.sort(exclude_indices))))
            
        # See if there are any regions we should always include
        if include is not None and len(include) > 0:
            include_indices = []
            
            if isinstance(include[0], float) and len(include) == 2:
                # Only two floats given, so we can only have one region to include
                include_indices.extend(range(*np.searchsorted(self.dispersion, include)))
                
            else:
                # Multiple regions provided
                for include_region in include:
                    include_indices.extend(range(*np.searchsorted(self.dispersion, include_region)))
        

        # We should exclude non-finite numbers from the fit
        non_finite_indices = np.where(~np.isfinite(self.flux))[0]
        continuum_indices = np.sort(list(set(continuum_indices).difference(non_finite_indices)))

        # We should also exclude zero or negative flux points from the fit
        zero_flux_indices = np.where(0 >= self.flux)[0]
        continuum_indices = np.sort(list(set(continuum_indices).difference(zero_flux_indices)))

        original_continuum_indices = continuum_indices.copy()

        if knot_spacing is None or knot_spacing == 0:
            knots = []

        else:
            knot_spacing = abs(knot_spacing)
            
            end_spacing = ((self.dispersion[-1] - self.dispersion[0]) % knot_spacing) /2.
        
            if knot_spacing/2. > end_spacing: end_spacing += knot_spacing/2.
                
            knots = np.arange(self.dispersion[0] + end_spacing, self.dispersion[-1] - end_spacing + knot_spacing, knot_spacing)
            if len(knots) > 0 and knots[-1] > self.dispersion[continuum_indices][-1]:
                knots = knots[:knots.searchsorted(self.dispersion[continuum_indices][-1])]
                
            if len(knots) > 0 and knots[0] < self.dispersion[continuum_indices][0]:
                knots = knots[knots.searchsorted(self.dispersion[continuum_indices][0]):]

        for iteration in xrange(max_iterations):
            
            splrep_disp = self.dispersion[continuum_indices]
            splrep_flux = self.flux[continuum_indices]

            splrep_weights = np.ones(len(splrep_disp))

            # We need to add in additional points at the last minute here
            if additional_points is not None and len(additional_points) > 0:

                for point, flux, weight in additional_points:

                    # Get the index of the fit
                    insert_index = int(np.searchsorted(splrep_disp, point))
                    
                    # Insert the values
                    splrep_disp = np.insert(splrep_disp, insert_index, point)
                    splrep_flux = np.insert(splrep_flux, insert_index, flux)
                    splrep_weights = np.insert(splrep_weights, insert_index, weight)

            if function == 'spline':
                order = 5 if order > 5 else order
                tck = interpolate.splrep(splrep_disp, splrep_flux,
                    k=order, task=-1, t=knots, w=splrep_weights)

                continuum = interpolate.splev(self.dispersion, tck)

            elif function in ("poly", "polynomial"):
            
                p = poly1d(polyfit(splrep_disp, splrep_flux, order))
                continuum = p(self.dispersion)

            else:
                raise ValueError("Unknown function type: only spline or poly available (%s given)" % (function, ))
            
            difference = continuum - self.flux
            sigma_difference = difference / np.std(difference[np.isfinite(self.flux)])

            # Clipping
            upper_exclude = np.where(sigma_difference > sigma_clip[1])[0]
            lower_exclude = np.where(sigma_difference < -sigma_clip[0])[0]
            
            exclude_indices = list(upper_exclude)
            exclude_indices.extend(lower_exclude)
            exclude_indices = np.array(exclude_indices)
            
            if len(exclude_indices) is 0: break
            
            exclusions.extend(exclude_indices)
            
            # Before excluding anything, we must check to see if there are regions
            # which we should never exclude
            if include is not None:
                exclude_indices = set(exclude_indices).difference(include_indices)
            
            # Remove regions that have been excluded
            continuum_indices = np.sort(list(set(continuum_indices).difference(exclude_indices)))
        
        # Snip the edges based on exclude regions
        if exclude is not None and len(exclude) > 0:

            # If there are exclusion regions that extend past the left_index/right_index,
            # then we will need to adjust left_index/right_index accordingly

            left_index = np.max([left_index, np.min(original_continuum_indices)])
            right_index = np.min([right_index, np.max(original_continuum_indices)])
            
        # Apply flux scaling
        continuum *= scale
        return continuum

        #return self.__class__(disp=self.dispersion[left_index:right_index], flux=continuum[left_index:right_index], headers=self.headers)
    

    
    def interpolate(self, new_disp, mode='linear', bounds_error=False,
                    fill_value=np.nan):
        """Interpolate the `Spectrum1D` onto a new dispersion map.
        
        Parameters
        ----
        new_disp : np.array
            An array of floating-point types containing the new dispersion points.
            
        mode : str
            Interpolation mode. See `scipy.interpolate.interp1d` for available
            options.
        
        bounds_error : bool
            See `scipy.interpolate.interp1d` for details.
        
        fill_value : float-type
            See `scipy.interpolate.interp1d`
        """
        
        f = interpolate.interp1d(self.disp, self.flux, kind=mode, copy=False,
                                 bounds_error=bounds_error, fill_value=fill_value)
        
        return self.__class__(new_disp, f(new_disp))
        


def ccf(x, y, axis=None):    
    """Computes the cross-correlation function of two series `x` and `y`.
    Note that the computations are performed on anomalies (deviations from 
    average).
    Returns the values of the cross-correlation at different lags.
    Lags are given as [0,1,2,...,n,n-1,n-2,...,-2,-1].
     
    :Parameters:
    `x` : 1D MaskedArray
        Time series.
    `y` : 1D MaskedArray
        Time series.
    `axis` : integer *[None]*
        Axis along which to compute (0 for rows, 1 for cols).
        If `None`, the array is flattened first.
    """

    if axis is None:    
        if x.ndim > 1:
            x = x.ravel()
            y = y.ravel()
        npad = x.size + y.size
        xanom = (x - x.mean(axis=None))
        yanom = (y - y.mean(axis=None))
        Fx = np.fft.fft(xanom, npad, )
        Fy = np.fft.fft(yanom, npad, )
        iFxy = np.fft.ifft(Fx.conj()*Fy).real
        varxy = np.sqrt(np.inner(xanom, xanom) * np.inner(yanom,yanom))

    else:
        npad = x.shape[axis] + y.shape[axis]
        if axis == 1:
            if x.shape[0] != y.shape[0]:
                raise ValueError("Arrays should have the same length!")
            xanom = (x - x.mean(axis=1)[:, None])
            yanom = (y - y.mean(axis=1)[:, None])
            varxy = np.sqrt((xanom*xanom).sum(1) * (yanom*yanom).sum(1))[:,None]
        else:
            if x.shape[1] != y.shape[1]:
                raise ValueError("Arrays should have the same width!")
            xanom = (x - x.mean(axis=0))
            yanom = (y - y.mean(axis=0))
            varxy = np.sqrt((xanom*xanom).sum(0) * (yanom*yanom).sum(0))
        Fx = np.fft.fft(xanom, npad, axis=axis)
        Fy = np.fft.fft(yanom, npad, axis=axis)
        iFxy = np.fft.ifft(Fx.conj()*Fy, n=npad, axis=axis).real
    
    return iFxy/varxy


def cross_correlate(observed, template, wl_region, full_output=False):
    """Performs a cross-correlation between the observed and template spectrum and
    provides a radial velocity and associated uncertainty.

    Parameters
    ----------
    observed : `Spectrum1D`
        The normalised observed spectrum.

    template : `Spectrum1D`
        The normalised template spectrum.

    wl_region : two length list containing floats [start, end]
        The starting and end wavelength to perform the cross-correlation on.

    full_output : `bool`, default False
        Whether or not to return the full output of the cross-correlation. If set to True
        then the output is as follows:

        v_rad, v_err, fft, profile

        where fft is a `np.ndarray` of shape (2, *) containing the Fourier transform
        and profile is a length 3 list containing the central peak point, peak height, and
        standard deviation.
    """

    if not isinstance(observed, Spectrum1D):
        raise TypeError("input observed spectrum must be a `specutils.Spectrum1D` object")

    if not isinstance(template, Spectrum1D):
        raise TypeError("template spectrum must be a `specutils.Spectrum1D` object")

    if not isinstance(wl_region, (tuple, list, np.ndarray)) or len(wl_region) != 2:
        raise TypeError("wavelength region must be a two length list-type")

    try:
        wl_region = map(float, wl_region)

    except:
        raise TypeError("wavelength regions must be float-like")

    # The following line of code will be supported until the end of the universe.
    c = 299792458e-3 # km/s

    # Splice the observed spectrum
    idx = np.searchsorted(observed.disp, wl_region)
    finite_values = np.isfinite(observed.flux[idx[0]:idx[1]])

    observed_slice = Spectrum1D(disp=observed.disp[idx[0]:idx[1]][finite_values], flux=observed.flux[idx[0]:idx[1]][finite_values])


    # Ensure the template and observed spectra are on the same scale
    template_func = interpolate.interp1d(template.disp, template.flux, bounds_error=False, fill_value=0.0)
    template_slice = Spectrum1D(disp=observed_slice.disp, flux=template_func(observed_slice.disp))

    # Perform the cross-correlation
    padding = observed_slice.flux.size + template_slice.flux.size
    x_norm = (observed_slice.flux - observed_slice.flux[np.isfinite(observed_slice.flux)].mean(axis=None))
    y_norm = (template_slice.flux - template_slice.flux[np.isfinite(template_slice.flux)].mean(axis=None))

    Fx = np.fft.fft(x_norm, padding, )
    Fy = np.fft.fft(y_norm, padding, )
    iFxy = np.fft.ifft(Fx.conj() * Fy).real
    varxy = np.sqrt(np.inner(x_norm, x_norm) * np.inner(y_norm, y_norm))

    fft_result = iFxy/varxy

    # Put around symmetry
    num = len(fft_result) - 1 if len(fft_result) % 2 else len(fft_result)

    fft_y = np.zeros(num)

    fft_y[:num/2] = fft_result[num/2:num]
    fft_y[num/2:] = fft_result[:num/2]

    fft_x = np.arange(num) - num/2

    # Get initial guess of peak
    p0 = np.array([fft_x[np.argmax(fft_y)], np.max(fft_y), 10])

    gaussian_profile = lambda p, x: p[1] * np.exp(-(x - p[0])**2 / (2.0 * p[2]**2))
    errfunc = lambda p, x, y: y - gaussian_profile(p, x)

    try:
        p1, ier = leastsq(errfunc, p0.copy(), args=(fft_x, fft_y))

    except:
        raise

    # Uncertainty
    sigma = np.mean(2.0*(fft_y.real)**2)**0.5

    # Create functions for interpolating back onto the dispersion map
    points = (0, p1[0], sigma)
    interp_x = np.arange(num/2) - num/4

    functions = []
    for point in points:
        idx = np.searchsorted(interp_x, point)
        f = interpolate.interp1d(interp_x[idx-3:idx+3], observed_slice.disp[idx-3:idx+3], bounds_error=True, kind='cubic')

        functions.append(f)

    # 0, p1, sigma
    f, g, h = [func(point) for func, point in zip(functions, points)]


    # Calculate velocity 
    measured_vrad = c * (1 - g/f)

    # Uncertainty
    measured_verr = np.abs(c * (1 - h/f))

    
    if full_output:
        results = [measured_vrad, measured_verr, np.vstack([fft_x, fft_y])]
        results.extend(p1)

        return results

    return [measured_vrad, measured_verr]


def compute_dispersion(aperture, beam, dispersion_type, dispersion_start,
    mean_dispersion_delta, num_pixels, redshift, aperture_low, aperture_high,
    weight=1, offset=0, function_type=None, order=None, Pmin=None, Pmax=None,
    coefficients=None):
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

        if function_type in (1, 2):
            # Chebyshev or Legendre polynomial.
            if None in (order, Pmin, Pmax, coefficients):
                raise TypeError("order, Pmin, Pmax and coefficients required "
                                "for a Chebyshev or Legendre polynomial")

            Pmean = (Pmax + Pmin)/2
            Pptp = Pmax - Pmin
            x = (np.arange(num_pixels) + 1 - Pmean)/(Pptp/2)
            p0 = np.ones(N_pixels)

            dispersion = coefficients[0] * p0 + coefficients[1] * p1
            for i in range(2, order):
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

    if mode not in ('average', 'sum', 'quadrature', ):
        raise ValueError, "Mode must be either 'average', 'quadrature' or 'sum'."


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

    return Spectrum1D(new_disp[idx_l:idx_r], new_flux[idx_l:idx_r], headers=headers)

