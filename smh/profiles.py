# coding: utf-8

""" Code to measure equivalent widths of atomic absorption lines """

__author__ = "Andy Casey <andy@astrowizici.st>"

# Standard library
import logging
import numpy as np
import operator

from itertools import groupby
from time import time

# Third party
from scipy import polyval, polyfit, optimize, integrate, c_
from scipy.special import wofz

__all__ = ['measure_transition_with_continuum_points', 'measure_transition']

logger = logging.getLogger(__name__)

def gaussian_profile(x, p, continuum):
    """Evaluates a Gaussian absorption profile across the
    `x` values with a given local continuum.

    Inputs
    ------
    x : `np.array`
        The values to evaluate the Gaussian profile on.

    p : list (length 3)
        The input parameters for the Gaussian profile:

        [position, width, amplitude]

    continuum : `np.array` or call-able function
        The continuum over the given `x` region.
    """

    pos, fwhm, amp = p

    # Is the continuum call-able?
    if hasattr(continuum, "__call__"): continuum = continuum(x)
    profile = amp * np.exp(-(x - pos)**2 / (2.0 * fwhm**2))

    return continuum - profile


def lorentzian_profile(x, p, continuum):
    """Evaluates a Lorentzian absorption profile across the
    `x` values with a given local continuum.

    L(x) = (1/pi) * (p[1]/((x - p[0])**2 + p[1]**2))
    
    Inputs
    ------
    x : `np.array`
        The values to evaluate the Lorentzian profile for.

    p : list (length 3)
        The input parameters for the profile:

        [pos, width, amp]

    continuum : `np.array` or callable function
        The continuum over the given `x` region.
    """

    pos, width, amp = p

    # Is the continuum call-able?
    if hasattr(continuum, "__call__"): continuum = continuum(x)

    profile = (amp/np.pi) * (width/((x - pos)**2 + width**2))

    return continuum - profile



def voigt_profile(x, p, continuum):
    """Evaluates a Voigt absorption profile across the `x`
    values with a given local continuum.

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))

    Inputs
    ------
    x : `np.array`
        The values to evaluate the Voigt profile for.

    p : list (length 4)
        The input parameters for the profile:

            [pos, fwhm, amp, shape]

    continuum : `np.array` or callable function
        The continuum over the given region.
    """

    pos, fwhm, amp, shape = p
    
    # Is the continuum call-able?
    if hasattr(continuum, '__call__'): continuum = continuum(x)
    n = len(x) if not isinstance(x, float) else 1

    profile = 1 / wofz(np.zeros((n)) + 1j * np.sqrt(np.log(2.0)) * shape).real
    profile = profile * amp * wofz(2*np.sqrt(np.log(2.0)) * (x - pos)/fwhm + 1j * np.sqrt(np.log(2.0))*shape).real

    return continuum - profile


def measure_transition_with_continuum_points(spectrum, line, fwhm, lhs_cont, rhs_cont):
    """Measure a line profile with two given continuum points; one on either side
    of the profile. This will only let you measure a line interactively with a
    Gaussian profile.
    
    Parameters
    ----------
    spectrum : `specutils.Spectrum1D` object
        The spectrum to measure the profile from.
        
    line : float-type
        The rest wavelength of the line to measure.
        
    fwhm : float-type
        Initial FWHM estimate of the line profile (in Angstroms).
        
    lhs_cont : tuple of float-type with length two
        Contains the (x, y) point for the left hand side continuum point.
        
    rhs_cont : tuple of float-types with length two
        Contains the (x, y) point for the right hand side continuum point.
    """
    
    if line + 1 > spectrum.disp[-1] or spectrum.disp[0] > line - 1:
        return [False, 0, 0, 'outside of spectral range']
        
    tol_profile = 1.0 # Angstroms
    tol_peak_difference = 0.3 # Angstroms
    integrate_sigmas = 5.
    
    # Get the line centroid region
    idx = np.searchsorted(spectrum.disp, [line - tol_peak_difference, tol_peak_difference + line])
    
    # Get the trough point
    trough = spectrum.flux[np.searchsorted(spectrum.disp, line)]
    
    # Get some wl_start and wl_end
    
    # Sort the continuum points (TODO - do we even need to do this?)
    if lhs_cont[0] > rhs_cont[0]:
        # The left hand side is actually on the right hand side
        rhs_x, rhs_y = lhs_cont
        lhs_x, lhs_y = rhs_cont
    else:
        rhs_x, rhs_y = rhs_cont
        lhs_x, lhs_y = lhs_cont
    
    wl_start, wl_end = lhs_x, rhs_x
    idx = []
    idx.append(np.searchsorted(spectrum.disp, wl_start, side='left'))
    idx.append(np.searchsorted(spectrum.disp, wl_end, side='right'))
    
    # Determine the continuum line
    m = (rhs_y - lhs_y)/(rhs_x - lhs_x)
    c = rhs_y - m * rhs_x
    
    continuum_x = spectrum.disp[idx[0]:idx[1]]
    continuum_y = m * continuum_x + c
    
    continuum_func = lambda x: m * x + c
            
    errfunc = lambda p, x, y, c: y - gaussian_profile(x, p, c)
    
    # Prepare initial guess and arguments for the least-sq optimization algorithm
    p0 = c_[line, fwhm, trough]
    args = (spectrum.disp[idx[0]:idx[1]], spectrum.flux[idx[0]:idx[1]], continuum_y)
    
    
    try:
        # Fit a Gaussian absorption profile to the data
        p1, cov_p, infodict, mesg, ier = optimize.leastsq(errfunc, p0.copy()[0], args=args, full_output=True, epsfcn=0.)
    
    except:
        
        # If we couldn't fit a line
        return [False, wl_start, wl_end, 'exception in fitting profile']
        
    else:

        # Was the fitting successful?
        if ier not in range(5):
            return [False, wl_start, wl_end, mesg]
        
        # Is the FWHM reasonable?
        if p1[1] > 15.:
            return [False, wl_start, wl_end, 'measured FWHM is greater than 15 Angstroms']
        
        # Force positive FWHM
        p1[1] = np.abs(p1[1])
        
        # Calculate chi squared value
        chi_sq = np.sum(pow(errfunc(p1, *args), 2)/(len(args[0]) - len(p1) - 1))
        
        # Integrate to +/- 5 sigmas
        idx = []
        idx.append(np.searchsorted(spectrum.disp, p1[0] - 5 * p1[1], side='left'))
        idx.append(np.searchsorted(spectrum.disp, p1[0] + 5 * p1[1], side='right'))

        continuum_x = spectrum.disp[idx[0]:idx[1]]
        continuum_y = m * continuum_x + c

        # Draw the profile
        fitted_y = gaussian_profile(continuum_x, p1, continuum_y)

        # Set the integrate function because it will be useful later
        def integrate_function(x, p, continuum):
            if hasattr(continuum, "__call__"): continuum = continuum(x)
            return continuum - gaussian_profile(x, p, continuum)

        # Calculate equivalent width
        # Integrate to get the equivalent width
        ew = integrate.quad(integrate_function, p1[0] - integrate_sigmas * p1[1], p1[0] + integrate_sigmas * p1[1], args=(p1, lambda x: m * x + c))
        
        # Convert profile sigma to FWHM: FWHM = profile_sigma * 2(2log2)^0.5
        p1[1] = p1[1] * 2 * (2 * np.log(2))**0.5

        p1 = list(p1)
        p1.insert(0, True)
        p1.extend([wl_start, wl_end, ew[0], chi_sq, ier, continuum_x, fitted_y, None])
        
        return p1


def measure_transition(spectrum, rest_wavelength, fwhm_guess, function='gaussian',
    local_continuum=True, use_central_weighting=False, window=10, order=2,
    detection_sigma=1.0, max_iter=4, shape_guess=1.0, **kwargs):
    """Fits a Gaussian absorption profile to the pixels in a given spectrum at the
    specified rest wavelength.

    Inputs
    ------
    spectrum : `Spectrum1D` object
        The spectrum to measure the line profile in.

    rest_wavelength : float-type
        The rest wavelength of the line to measure.

    fwhm_guess : float-type
        The initial FWHM guess for the absorption line. The fitting routine is fairly
        insensitive to this guess, but having an accurate guess will result in a faster
        fit to the data.

    function : string, either 'gaussian', 'voigt', or 'lorentzian'
        The function to fit the absorption profile with.

    local_continuum : bool, default True
        Whether to determine the local continuum surrounding the line when measuring
        the absorption profile. If this is not turned on, the continuum is assumed to
        be at unity.

    use_central_weighting : bool, default False
        Whether to weight the central line pixels more than the surrounding regions.

    window : float, default 10 (Angstroms)
        Window on either side of the line to use during the local continuum fitting stage.
        This is mute if `local_continuum` is not True.

    order : int, greater than 0
        The polynomial  order to use during the local continuum fitting process. Mute
        if `local_continuum` is False.

    detection_sigma : float, positive
        The standard deviation necessary to identify groups of neighbouring pixels that
        require us to fit a Gaussian absorption line. Mute if `local_continuum` is False.

    max_iter : int, positive
        The maximum number of iterations to perform during local continuum fitting. Mute
        if `local_continuum` is False.

    shape_guess : float
        The initial guess for the profile shape. This value is only relevant for Voight
        profile fitting.


    Keyword arguments
    -----------------
    exclusion_ranges : `list` of 2-length tuples
        Specifies ranges of spectrum to exclude (in Angstroms). These segments of spectra
        will be excluded from the continuum determination.

    search_within_exclusion_ranges : bool
        If `exclusion_ranges` has been specified, then setting this to True will ensure
        that we don't re-look within the range of exclusion ranges for new regions to
        fit profiles for or exclude.

    """

    function = function.lower()
    if function not in ("gaussian", "voigt", "lorentzian"):
        raise ValueError("invalid fitting function provided. Must be gaussian, voigt or lorentzian.")
   
    if rest_wavelength + 1 > spectrum.disp[-1] or spectrum.disp[0] > rest_wavelength - 1:
        return [False, 0, 0, "outside of observed spectral range"]
    
    # Here's where the magic happens
    tol_profile = 1.0 # Angstroms
    tol_peak_difference = 0.3 # Angstroms
    integrate_sigmas = 5.

    profile_sigma_scale = 2 * (2 * np.log(2))**0.5
    
    if local_continuum:
        # Should we determine the local continuum? 

        # Determine the region surrounding the line profile
        idx = np.searchsorted(spectrum.disp,
            [
                rest_wavelength - 3 * fwhm_guess * profile_sigma_scale - window,
                rest_wavelength - 3 * fwhm_guess * profile_sigma_scale,
                rest_wavelength + 3 * fwhm_guess * profile_sigma_scale,
                rest_wavelength + 3 * fwhm_guess * profile_sigma_scale + window
            ])
        
        cont_flux = []
        cont_flux.extend(spectrum.flux[idx[0]:idx[1]])
        cont_flux.extend(spectrum.flux[idx[2]:idx[3]])
        
        # First pass at the continuum level
        continuum = np.median(cont_flux)

        exclusion_ranges = []

    else:
        # Assume perfectly normalised
        continuum = 1.0
        exclusion_ranges = None
        
        
    # We are looking for the rest_wavelength within +/- tol_peak_difference (Angstroms)
    idx = np.searchsorted(spectrum.disp, [rest_wavelength - tol_peak_difference, tol_peak_difference + rest_wavelength])

    # Get minima point at the provided rest_wavelength position as a first-guess
    trough = spectrum.flux[np.searchsorted(spectrum.disp, rest_wavelength)]
    
    # Get start and ending wavelength

    # This following line must be from an earlier attempt. Leaving here just in case you wonder "what changed one day"
    #idx = np.searchsorted(spectrum.disp, [rest_wavelength - tol_profile, tol_profile + rest_wavelength])

    idx = np.searchsorted(spectrum.disp, [rest_wavelength - integrate_sigmas * fwhm_guess, integrate_sigmas * fwhm_guess + rest_wavelength])

    wl_start, wl_end = spectrum.disp[idx[0]], spectrum.disp[idx[1] - 1]
    
    p0 = c_[rest_wavelength, fwhm_guess, trough]
    args = (spectrum.disp[idx[0]:idx[1]], spectrum.flux[idx[0]:idx[1]], continuum)

    # Set the profile function and update the initial guess if necessary
    if function == 'gaussian':
        profile_function = gaussian_profile

    elif function == "voigt":
        profile_function = voigt_profile

        # We require an additional 'shape' input for the voigt function
        p0 = c_[rest_wavelength, fwhm_guess, trough, shape_guess]

    elif function == "lorentzian":
        profile_function = lorentzian_profile

    else:
        return [False, wl_start, wl_end, "this is not the profile you're looking for"]

    # Set the integrate function because it will be useful later
    def integrate_function(x, p, continuum):
        if hasattr(continuum, "__call__"): continuum = continuum(x)
        return continuum - profile_function(x, p, continuum)

    # Define the error function based on whether we should use central weighting or not
    if use_central_weighting:
        errfunc = lambda p, x, y, c: (y - profile_function(x, p, c)) * (1 + np.exp(-(x - p[0])**2 / (4.0 * p[1]**2)))

    else:
        errfunc = lambda p, x, y, c: y - profile_function(x, p, c)


    # Now we can fit to our heart's content!    
    try:
        p1, cov_p, infodict, mesg, ier = optimize.leastsq(errfunc, p0.copy()[0], args=args, epsfcn=0.0, full_output=True)
        
    except:
        # If we couldn't fit a rest_wavelength
        return [False, wl_start, wl_end, 'exception in first pass']
        
    else:
        
        # Was the fitting successful?
        if ier not in range(5):
            return [False, wl_start, wl_end, mesg]
        
        # Is the fwhm_guess reasonable?
        if p1[1] > 2.:
            return [False, wl_start, wl_end, 'measured FWHM was greater than 2 Angstroms']
        
        # Is the answer reasonable?
        if np.abs(p1[0] - rest_wavelength) > tol_peak_difference:
            return [False, wl_start, wl_end, 'measured line centroid too far from initial guess: |%1.2f| > %1.2f' % (p1[0] - rest_wavelength, tol_peak_difference, )]
        
        # Force positive FWHM
        p1[1] = np.abs(p1[1])
        

        # We will integrate from +/- `integrate_sigmas` sigma of the rest_wavelength
        x = np.arange(p1[0] - integrate_sigmas * p1[1], p1[0] + integrate_sigmas * p1[1], np.min(np.diff(spectrum.disp)) / 4.)
        
        if local_continuum:

            # Get 3-sigma + the continuum window region on either side.
            cont_lhs, cont_rhs = np.searchsorted(spectrum.disp, [
                p1[0] - integrate_sigmas * p1[1] - window,
                p1[0] + integrate_sigmas * p1[1] + window
                ])

            # Indices will store which points to keep when determining the continuum
            # indices is referenced from 0 on spectrum.wave
            indices = list(range(cont_lhs, cont_rhs))
            
            # Exclude possible continuum points that are +/- 3-sigma from our currently measured rest_wavelength 
            indices = set(indices).difference(range(*np.searchsorted(spectrum.disp,
                [   p1[0] - integrate_sigmas * p1[1],
                    integrate_sigmas * p1[1] + p1[0]
                ]
                )))

            # Exclude continuum points specified in kwargs
            max_exclusion_ranges = None
            if 'exclusion_ranges' in kwargs.keys():

                for region in kwargs['exclusion_ranges']:
                    exclude = np.arange(*np.searchsorted(spectrum.disp, region))

                    indices = indices.difference(exclude)

                # Should we be searching within the exclusion regions? Default is yes.
                if 'search_within_exclusion_ranges' in kwargs.keys() \
                and not kwargs['search_within_exclusion_ranges'] \
                and len(kwargs['exclusion_ranges']) > 0:
                    all_exclusion_ranges = np.array(kwargs['exclusion_ranges']).flatten()
                    max_exclusion_ranges = np.min(all_exclusion_ranges), np.max(all_exclusion_ranges)

            indices = np.sort(list(indices))

            # Remove non-finite fluxes
            indices = indices[np.isfinite(spectrum.flux[indices])]

            # Iterate on the continuum fitting
            for iteration in xrange(max_iter):

                if len(indices) == 0:
                    return [False, 0, 0, "outside of observed spectral range"]

                # Fit a polynomial to the continuum points
                coeffs = polyfit(spectrum.disp[indices], spectrum.flux[indices], order)
                continuum = polyval(coeffs, spectrum.disp[indices])
                
                # Detect outliers based on continuum_rest_wavelength_detection (min. sigma)
                sigmas = (spectrum.flux[indices] - continuum)/np.std(spectrum.flux[indices])
                outliers = np.where(sigmas < -detection_sigma)[0]
                
                # Since indices is not a consecutive list, we need to reference
                # outliers back to spectrum.wave
                outliers = indices[outliers]
                
                # Let's exclude obvious outlier points before trying to fit profiles
                indices = np.sort(list(set(indices).difference(outliers)))
                
                # We need to group consecutive outlier points together and only fit
                # profiles once to each consecutive group
                outlier_groups = []
                for k, g in groupby(enumerate(outliers), lambda (i,x):i-x):
                    region = map(operator.itemgetter(1), g)

                    if max_exclusion_ranges is not None:
                        # Is this region within the pre-defined exclusion region?
                        if max_exclusion_ranges[1] >= region[0] >= max_exclusion_ranges[0] \
                        or (len(region) > 1 and max_exclusion_ranges[1] >= region[-1] >= max_exclusion_ranges[0]):
                            continue

                    outlier_groups.append(region)
                

                # See if we can fit a profile to each rest_wavelength
                for i, outlier_group in enumerate(outlier_groups):
                    
                    # We can guess that the rest_wavelength trough is in the middle of this outlier group points
                    # (Unless we have two blended rest_wavelengths next to each other)
                    
                    centroid = spectrum.disp[outlier_group[int(len(outlier_group) / 2)]]
                    
                    # Look for the rest_wavelength profile +/- 0.3 Angstroms from here
                    start, end = centroid - tol_peak_difference, centroid + tol_peak_difference
                    
                    trough = spectrum.flux[outlier_group[int(len(outlier_group) / 2)]]
                    
                    rest_wavelength_p0 = c_[centroid, fwhm_guess, trough]
                    
                    idx_ = np.searchsorted(spectrum.disp, [start, end])
                    args = (spectrum.disp[idx_[0]:idx_[1]], spectrum.flux[idx_[0]:idx_[1]], continuum)
                    
                    try: rest_wavelength_p1 = optimize.leastsq(errfunc, rest_wavelength_p0.copy()[0], args=args, full_output=True, epsfcn=0.)
                    except: continue
                    else:
                        
                        # Some criterion for being a good rest_wavelength
                        if rest_wavelength_p1[-1] < 5 and np.abs(rest_wavelength_p1[0][2]) < 0.5:
                            
                            # Grow out the exclusion region by +/- 3 sigmas of the fitted rest_wavelength
                            exclude = np.arange(*np.searchsorted(spectrum.disp, [rest_wavelength_p1[0][0] - 3 * np.abs(rest_wavelength_p1[0][1]), rest_wavelength_p1[0][0] + 3 * np.abs(rest_wavelength_p1[0][1])]))
                            indices = np.sort(list(set(indices).difference(exclude)))
                
                            # Add this exclusion range.
                            exclusion_ranges.append([
                                rest_wavelength_p1[0][0] - 3 * np.abs(rest_wavelength_p1[0][1]),
                                rest_wavelength_p1[0][0] + 3 * np.abs(rest_wavelength_p1[0][1])
                                ])
            
            # Fit our final continuum
            coeffs = polyfit(spectrum.disp[indices], spectrum.flux[indices], order)

            cont_lhs = int(np.max([cont_lhs, indices[0]]))
            cont_rhs = int(np.min([cont_rhs, indices[-1]]))
            
            continuum_x = spectrum.disp[range(cont_lhs, cont_rhs)]
            continuum_y = polyval(coeffs, continuum_x)
            
            continuum = lambda x: polyval(coeffs, x)

        
        else:
            # Assuming perfect normalisation
            continuum_x = x
            continuum_y = [1] * len(x)

        idx = np.searchsorted(spectrum.disp, [rest_wavelength - integrate_sigmas * p1[1], integrate_sigmas * p1[1] + rest_wavelength])
        
        # Re-calculate the initial rest_wavelength
        p0 = c_[p1[0], fwhm_guess, p1[2]]

        # If the function is a Voigt, we need a shape parameter too.
        if function == "voigt":
            p0 = c_[rest_wavelength, fwhm_guess, trough, shape_guess]

        # Update the arguments with our best guess of the continuum
        args = (spectrum.disp[idx[0]:idx[1]], spectrum.flux[idx[0]:idx[1]], continuum)

        try:
            # Update our fit of the profile
            p1, cov_p, infodict, mesg, ier = optimize.leastsq(errfunc, p0.copy()[0], args=args, full_output=True, epsfcn=0.)

        except:
            # If we couldn't fit a rest_wavelength
            return [False, wl_start, wl_end, 'exception in second pass']

        else:
            
            # Was the fitting successful?
            if ier not in range(5):
                return [False, wl_start, wl_end, mesg]
            
            # Is the answer reasonable?
            if np.abs(p1[0] - rest_wavelength) > tol_peak_difference:
                return [False, wl_start, wl_end, 'rest_wavelength pos tol |%1.2f| > %1.2f' % (p1[0] - rest_wavelength, tol_peak_difference, )]
            
            # Integrate to get the equivalent width
            ew = integrate.quad(integrate_function, p1[0] - integrate_sigmas * p1[1], p1[0] + integrate_sigmas * p1[1], args=(p1, continuum))
            
            # For the moment, leave chi_sq as nan
            chi_sq = np.sum(pow(errfunc(p1, *args), 2)/(len(args[0]) - len(p1) - 1))
            
            if ew[0] < 0:
                return [False, wl_start, wl_end, 'negative equivalent width (%1.2f)' % (ew[0], )]
            
            # We want to return the full profile across continuum_x
            fitted_y = profile_function(continuum_x, p1, continuum_y)

            p1 = list(p1)
            
            # If this is a Gaussian...
            # Convert profile sigma to FWHM: FWHM = profile_sigma * 2(2log2)^0.5
            if function == "gaussian":
                p1[1] = p1[1] * profile_sigma_scale

            p1.insert(0, True)
            
            # Check for the proper format in exclusion_ranges
            if exclusion_ranges is not None and len(exclusion_ranges) == 2 and isinstance(exclusion_ranges[0], float):
                exclusion_ranges = [exclusion_ranges]
            
            p1.extend([wl_start, wl_end, ew[0], chi_sq, ier, continuum_x, fitted_y, exclusion_ranges])
            return p1
