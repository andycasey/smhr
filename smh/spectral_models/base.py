#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Models for fitting spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["BaseSpectralModel"]

import numpy as np

from .quality_constraints import constraints
from ..linelists import LineList
from ..robust_polyfit import polyfit as rpolyfit
from astropy.table import Row
from smh.photospheres.abundances import asplund_2009 as solar_composition
from scipy import optimize as op
from scipy.linalg import svd, cholesky, solve_triangular, LinAlgError, inv

class BaseSpectralModel(object):

    def __init__(self, session, transitions, **kwargs):
        """
        Initialize a base class for modelling spectra.

        :param session:
            The session that this spectral model will be associated with.

        :param transitions:
            A linelist containing atomic data for this model.
        """

        if isinstance(transitions, Row):
            transitions = LineList(transitions)
        assert isinstance(transitions, LineList), type(transitions)
        

        self._session = session
        self._transitions = transitions

        self.metadata = {
            "is_upper_limit": False,
            "use_for_stellar_composition_inference": True,
            "use_for_stellar_parameter_inference": (
                "Fe I" in self.transitions["element"] or
                "Fe II" in self.transitions["element"]),
            "antimask_flag": False,
            "user_flag": 0
        }

        # Create a _repr_wavelength property.
        if len(self.transitions) == 1:
            self._repr_wavelength \
                = "{0:.1f}".format(self.transitions["wavelength"][0])
        else:
            self._repr_wavelength \
                = "~{0:.0f}".format(np.mean(self.transitions["wavelength"]))
        
        return None


    @property
    def wavelength(self):
        """
        Return a (sometimes approximate) wavelength for where this spectral line
        occurs.
        """

        if len(self.transitions) == 1: return float(self.transitions["wavelength"])
        if hasattr(self,"_wavelength"):
            return self._wavelength
        wavelength = np.mean(self.transitions["wavelength"])
        return int(wavelength)


    @property
    def is_acceptable(self):
        """ Return whether this spectral model is acceptable. """
        return self.metadata.get("is_acceptable", False)


    @is_acceptable.setter
    def is_acceptable(self, decision):
        """
        Mark the spectral model as acceptable or unacceptable.

        :param decision:
            A boolean flag.
        """
        decision = bool(decision)
        if not decision or (decision and "fitted_result" in self.metadata):
            self.metadata["is_acceptable"] = bool(decision)
        return None


    @property
    def is_upper_limit(self):
        """ Return whether this spectral model is acceptable. """
        return self.metadata.get("is_upper_limit", False)


    @is_upper_limit.setter
    def is_upper_limit(self, decision):
        """
        Mark the spectral model as acceptable or unacceptable.

        :param decision:
            A boolean flag.
        """
        self.metadata["is_upper_limit"] = bool(decision)
        return None


    @property
    def use_for_stellar_parameter_inference(self):
        """
        Return whether this spectral model should be used during the
        determination of stellar parameters.
        """
        return self.metadata["use_for_stellar_parameter_inference"]


    @use_for_stellar_parameter_inference.setter
    def use_for_stellar_parameter_inference(self, decision):
        """
        Mark whether this spectral model should be used when inferring stellar
        parameters.

        :param decision:
            A boolean flag.
        """
        self.metadata["use_for_stellar_parameter_inference"] = bool(decision)
        return None


    @property
    def use_for_stellar_composition_inference(self):
        """
        Return whether this spectral model should be used for the determination
        of stellar composition.
        """
        return self.metadata["use_for_stellar_composition_inference"]



    @use_for_stellar_composition_inference.setter
    def use_for_stellar_composition_inference(self, decision):
        """
        Mark whether this spectral model should be used when inferring the 
        stellar composition.

        :param decision:
            A boolean flag.
        """
        self.metadata["use_for_stellar_composition_inference"] = bool(decision)
        return None


    def apply_quality_constraints(self, quality_constraints):
        """
        Apply quality constraints to the this spectral model. If the model does
        not meet the quality constraints it will be marked as unacceptable.

        :param quality_constraints:
            A dictionary containing constraint names as keys and a 2-length
            tuple with the (lower, upper) bounds specified as values.

        :returns:
            Whether this model met the specified quality constraints.
        """

        is_ok = constraints(self, quality_constraints)
        self.is_acceptable = is_ok
        return is_ok


    def meets_quality_constraints(self, quality_constraints):
        """
        Check whether this spectral model meets specified quality constraints.

        :param quality_constraints:
            A dictionary containing constraint names as keys and a 2-length
            tuple with the (lower, upper) bounds specified as values.

        :returns:
            Whether this model met the specified quality constraints.
        """
        return constraints(self, quality_constraints)


    @property
    def meets_quality_constraints_in_parent_session(self):
        """
        Returns whether this spectral model meets the quality constraints set
        in the parent session.
        """

        return self.meets_quality_constraints(
            self._session.setting(("spectral_model_quality_constraints", )) or {})


    @property
    def transitions(self):
        """ Return the transitions associateed with this class. """

        # This is left as a property to prevent users from arbitrarily setting
        # the .transitions attribute.
        return self._transitions


    @property
    def elements(self):
        """ Return the elements to be measured from this class. """
        return self.metadata["elements"]


    @property
    def num_elems(self):
        return len(self.elements)


    @property
    def species(self):
        """ Return the species to be measured from this class. """
        return self.metadata["species"]


    @property
    def session(self):
        """ Return the parent session that this model is associated with. """
        return self._session

    @property
    def abundances(self):
        """ Return abundances if fit, else None """
        try:
            return self.metadata["fitted_result"][2]["abundances"]
        except KeyError:
            return None
        
    @property
    def abundances_to_solar(self):
        """ Return [X/H] if fit, else None """
        abunds = self.abundances
        if abunds is None: return None
        elems = self.elements
        return [AX - solar_composition(X) for AX,X in zip(abunds,elems)]

    @property
    def abundance_uncertainties(self):
        return None

    @property
    def expot(self):
        raise NotImplementedError
    
    @property
    def loggf(self):
        raise NotImplementedError
    
    @property
    def equivalent_width(self):
        return None

    @property
    def equivalent_width_uncertainty(self):
        return None

    @property
    def reduced_equivalent_width(self):
        return None

    @property
    def user_flag(self):
        return self.metadata["user_flag"]

    @user_flag.setter
    def user_flag(self, flag):
        self.metadata["user_flag"] = flag
        return None

    @property
    def parameters(self):
        """
        Return the model parameters.
        This must be implemented by the sub-classes.
        """
        raise NotImplementedError(
            "parameters property must be implemented by the sub-classes")


    @property
    def parameter_bounds(self):
        """ Return the fitting limits on the parameters. """
        return self._parameter_bounds


    @property
    def parameter_names(self):
        """ Return the model parameter names. """
        return self._parameter_names
    
    @property
    def snr(self):
        """ Return the SNR: median of sqrt(ivar) for masked pixels """
        spectrum = self._verify_spectrum(None)
        return np.nanmedian(spectrum.ivar[self.mask(spectrum)]**0.5)
    
    def __call__(self, dispersion, *args, **kwargs):
        """ The data-generating function. """
        raise NotImplementedError(
            "the data-generating function must be implemented by sub-classes")


    def __getstate__(self):
        """ Return a serializable state of this spectral model. """

        # LineList is a superclass of astropy.table.Table so it serializes
        # This could be a problem if e.g. astropy version changes
        state = {
            "type": self.__class__.__name__,
            "transitions": self.transitions.as_array(),
            "metadata": self.metadata
        }
        return state


    def __setstate__(self, state):
        """ Disallow the state to be instantiated from a serialised object. """
        return None


    def _verify_transitions(self):
        """
        Verify that the transitions provided are valid.
        """
        # TODO
        return True


    def _verify_spectrum(self, spectrum):
        """
        Check that the spectrum provided is valid and has data in the wavelength
        range that we are interested.

        :param spectrum:
            The observed rest-frame normalized spectrum.
        """

        spectrum = spectrum or self.session.normalized_spectrum

        # Check the transition is in the spectrum range.
        wavelength = self.transitions["wavelength"]
        try:
            wavelength = wavelength[0]
        except IndexError:
            None
        if wavelength + 1 > spectrum.dispersion[-1] \
        or wavelength - 1 < spectrum.dispersion[0]:
            raise ValueError(
                "the observed spectrum contains no data over the wavelength "
                "range we require")

        return spectrum


    def mask(self, spectrum):
        """
        Return a pixel mask based on the metadata and existing mask information
        available.

        :param spectrum:
            A spectrum to generate a mask for.
        """

        if self.metadata["antimask_flag"]:
            antimask = np.ones_like(spectrum.dispersion,dtype=bool)
            for start, end in self.metadata["mask"]:
                antimask *= ~((spectrum.dispersion >= start) \
                            * (spectrum.dispersion <= end))
            return ~antimask

        window = abs(self.metadata["window"])
        wavelengths = self.transitions["wavelength"]
        try:
            lower_wavelength = wavelengths[0]
            upper_wavelength = wavelengths[-1]
        except IndexError:
            # Single row.
            lower_wavelength, upper_wavelength = (wavelengths, wavelengths)

        mask = (spectrum.dispersion >= lower_wavelength - window) \
             * (spectrum.dispersion <= upper_wavelength + window)

        # Any masked ranges specified in the metadata?
        for start, end in self.metadata["mask"]:
            mask *= ~((spectrum.dispersion >= start) \
                     * (spectrum.dispersion <= end))

        return mask


    def fitting_function(self, dispersion, *parameters):
        """
        Generate data at the dispersion points, given the parameters, but
        respect the boundaries specified on model parameters.

        :param dispersion:
            An array of dispersion points to calculate the data for.

        :param parameters:
            Keyword arguments of the model parameters and their values.
        """

        for parameter_name, (lower, upper) in self.parameter_bounds.items():
            value = parameters[self.parameter_names.index(parameter_name)]
            if not (upper >= value and value >= lower):
                return np.nan * np.ones_like(dispersion)

        return self.__call__(dispersion, *parameters)


    def _fill_masked_arrays(self, spectrum, x, *y):
        """
        Detect masked regions and fill masked regions in y-axis arrays with
        NaNs.

        :param spectrum:
            The spectrum used in the fit.

        :param x:
            The x values that were used in the fit.

        :param *y:
            The y-axis arrays to fill.
        """

        indices = spectrum.dispersion.searchsorted(x)
        x_actual = spectrum.dispersion[indices[0]:1 + indices[-1]]

        filled_arrays = [x_actual]
        for yi in y:
            yi_actual = np.nan * np.ones_like(x_actual)
            if len(yi_actual.shape) == 2:
                yi_actual[:, indices - indices[0]] = yi
            else:
                yi_actual[indices - indices[0]] = yi
            filled_arrays.append(yi_actual)

        return tuple(filled_arrays)

    @property
    def reduced_chi2(self):
        try:
            (named_p_opt, cov, meta) = self.metadata["fitted_result"]
            return meta["chi_sq"] / meta["dof"]
        except:
            return np.nan
            
    @property
    def chi2(self):
        try:
            (named_p_opt, cov, meta) = self.metadata["fitted_result"]
            return meta["chi_sq"]
        except:
            return np.nan
            
    @property
    def dof(self):
        try:
            (named_p_opt, cov, meta) = self.metadata["fitted_result"]
            return meta["dof"]
        except:
            return np.nan
            
    @property
    def residual_slope(self):
        try:
            (named_p_opt, cov, meta) = self.metadata["fitted_result"]
            x = meta["model_x"]
            y = meta["residual"]
            coeff, scat = rpolyfit(x, y, 1)
            return coeff[0]
        except:
            return np.nan
    
    @property
    def residual_slope_and_err(self):
        try:
            (named_p_opt, cov, meta) = self.metadata["fitted_result"]
            x = meta["model_x"]
            y = meta["residual"]
            coeff, scat = rpolyfit(x, y, 1)
            yfit = np.polyval(coeff, x)
            xmean = np.mean(x)
            N = len(x)
            merr = np.sqrt(np.sum((y-yfit)**2)/((N-2)*np.sum((x-xmean)**2)))
            return coeff[0], merr
        except:
            return np.nan, np.nan

def penalized_curve_fit_lm(f, xdata, ydata,
                           penalty_function=None,
                           p0=None, sigma=None, absolute_sigma=False,
                           check_finite=True,
                           *, full_output=False, **kwargs):
    """
    Wrap scipy.optimize.curve_fit using lm method, but add penalty function on the parameters p0
    Enables adding regularization/prior to parameters p0
    This means we cannot specify bounds, as we use lm, but they can be updated with the penalty function
    lm also numerically determines the derivative itself
    Uses scipy.optimize.least_squares as the underlying minimizer (instead of leastsq as in default curve_fit for lm)
    # Uses scipy.optimize.leastsq as the underlying minimizer
    
    Assumes ``ydata = f(xdata, *params) + eps``.
    
    This means we minimize ((ydata-f(xdata, *params))/sigma)^2 + penalty_function(params)
    i.e. penalty_function(params) = scalar
    
    """
    
    ## creating a function that returns ((y-f(x,p))/sigma)^2 + penalty_function(p)
    def _wrap_func(func, xdata, ydata, transform, penalty_function):
        if transform is None:
            def func_wrapped(params):
                return (func(xdata, *params) - ydata)**2 + penalty_function(params)
        elif transform.ndim == 1:
            def func_wrapped(params):
                return (transform * (func(xdata, *params) - ydata))**2 + penalty_function(params)
        else:
            # Chisq = (y - yd)^T C^{-1} (y-yd)
            # transform = L such that C = L L^T
            # C^{-1} = L^{-T} L^{-1}
            # Chisq = (y - yd)^T L^{-T} L^{-1} (y-yd)
            # Define (y-yd)' = L^{-1} (y-yd)
            # by solving
            # L (y-yd)' = (y-yd)
            # and minimize (y-yd)'^T (y-yd)'
            def func_wrapped(params):
                return (solve_triangular(transform, func(xdata, *params) - ydata, lower=True))**2 + penalty_function(params)
        return func_wrapped
    
    ## This next code is taken from scipy.optimize.curve_fit v1.9.3, removing the bounds parts
    if p0 is None:
        from scipy._lib._util import getfullargspec_no_self as _getfullargspec
        # determine number of parameters by inspecting the function
        sig = _getfullargspec(f)
        args = sig.args
        if len(args) < 2:
            raise ValueError("Unable to determine number of fit parameters.")
        n = len(args) - 1
    else:
        p0 = np.atleast_1d(p0)
        n = p0.size

    if penalty_function is None:
        def penalty_wrapped(params):
            return 0.
    else:
        penalty_wrapped = penalty_function

    # NaNs cannot be handled
    if check_finite:
        ydata = np.asarray_chkfinite(ydata, float)
    else:
        ydata = np.asarray(ydata, float)

    if isinstance(xdata, (list, tuple, np.ndarray)):
        # `xdata` is passed straight to the user-defined `f`, so allow
        # non-array_like `xdata`.
        if check_finite:
            xdata = np.asarray_chkfinite(xdata, float)
        else:
            xdata = np.asarray(xdata, float)

    if ydata.size == 0:
        raise ValueError("`ydata` must not be empty!")

    # Determine type of sigma
    if sigma is not None:
        sigma = np.asarray(sigma)

        # if 1-D, sigma are errors, define transform = 1/sigma
        if sigma.shape == (ydata.size, ):
            transform = 1.0 / sigma
        # if 2-D, sigma is the covariance matrix,
        # define transform = L such that L L^T = C
        elif sigma.shape == (ydata.size, ydata.size):
            try:
                # scipy.linalg.cholesky requires lower=True to return L L^T = A
                transform = cholesky(sigma, lower=True)
            except LinAlgError as e:
                raise ValueError("`sigma` must be positive definite.") from e
        else:
            raise ValueError("`sigma` has incorrect shape.")
    else:
        transform = None


    ## This is where the difference comes in!
    func = _wrap_func(f, xdata, ydata, transform, penalty_wrapped)
    
    # if ydata.size == 1, this might be used for broadcast.
    if ydata.size != 1 and n > ydata.size:
        raise TypeError(f"The number of func parameters={n} must not"
                        f" exceed the number of data points={ydata.size}")
    
    ## This is the original curve_fit code
    #res = op.leastsq(func, p0, Dfun=None, full_output=1, **kwargs)
    #popt, pcov, infodict, errmsg, ier = res
    #if ier not in [1, 2, 3, 4]:
    #    raise RuntimeError("Optimal parameters not found: " + errmsg)
    
    ## This is adapting the curve_fit code to try and reproduce optimize.leastsq
    if "ftol" in kwargs: ftol = kwargs.pop("ftol")
    else: ftol = 1.49012e-8
    if "xtol" in kwargs: xtol = kwargs.pop("xtol")
    else: xtol = 1.49012e-8
    if "gtol" in kwargs: gtol = kwargs.pop("gtol")
    else: gtol = 1.0e-8
    res = op.least_squares(func, p0, method='lm',
                           ftol=ftol, xtol=xtol, gtol=gtol,
                           **kwargs)
    if not res.success:
        raise RuntimeError("Optimal parameters not found: " + res.message)

    infodict = dict(nfev=res.nfev, fvec=res.fun)
    ier = res.status
    errmsg = res.message

    ysize = len(res.fun)
    cost = 2 * res.cost  # res.cost is half sum of squares!
    popt = res.x

    # Do Moore-Penrose inverse discarding zero singular values.
    _, s, VT = svd(res.jac, full_matrices=False)
    threshold = np.finfo(float).eps * max(res.jac.shape) * s[0]
    s = s[s > threshold]
    VT = VT[:s.size]
    pcov = np.dot(VT.T / s**2, VT)

    ysize = len(infodict['fvec'])
    cost = np.sum(infodict['fvec'] ** 2)

    warn_cov = False
    if pcov is None:
        # indeterminate covariance
        pcov = zeros((len(popt), len(popt)), dtype=float)
        pcov.fill(inf)
        warn_cov = True
    elif not absolute_sigma:
        if ysize > p0.size:
            s_sq = cost / (ysize - p0.size)
            pcov = pcov * s_sq
        else:
            pcov.fill(inf)
            warn_cov = True

    if warn_cov:
        warnings.warn('Covariance of the parameters could not be estimated',
                      category=OptimizeWarning)

    if full_output:
        return popt, pcov, infodict, errmsg, ier
    else:
        return popt, pcov
