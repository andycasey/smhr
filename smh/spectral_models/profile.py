#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model for fitting an absorption profile to spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["ProfileFittingModel"]

import logging
import numpy as np
import scipy.optimize as op
import astropy.table
from astropy.table import Row
from astropy.constants import c as speed_of_light
from collections import OrderedDict
from scipy.special import wofz
from scipy import integrate

from .base import BaseSpectralModel
from ..specutils import Spectrum1D
from ..linelists import LineList

logger = logging.getLogger(__name__)


def _gaussian(x, *parameters):
    """
    Evaluate a Gaussian profile at x, given the profile parameters.

        y = amplitude * exp(-(x - position)**2 / (2.0 * sigma**2))

    :param x:
        The x-values to evaluate the Gaussian profile at.

    :param parameters:
        The position, sigma, and amplitude of the Gaussian profile.
    """
    position, sigma, amplitude = parameters
    return amplitude * np.exp(-(x - position)**2 / (2.0 * sigma**2))


def _lorentzian(x, *parameters):
    """
    Evaluate a Lorentzian profile at x, given the profile parameters:

        y = (amplitude/PI) * (width/((x - positions)**2 + width**2))

    :param x:
        The x-values to evaluate the Lorentzian profile at.

    :param parameters:
        The position, width, and amplitude of the Lorentzian profile.
    """
    position, width, amplitude = parameters
    return (amplitude/np.pi) * (width/((x - position)**2 + width**2))


def _voigt(x, *parameters):
    """
    Evaluate a Voigt profile at x, given the profile parameters.

    :param x:
        The x-values to evaluate the Voigt profile at.

    :param parameters:
        The position, fwhm, amplitude, and shape of the Voigt profile.
    """
    try:
        n = len(x)
    except TypeError:
        n = 1

    position, fwhm, amplitude, shape = parameters
    
    profile = 1 / wofz(np.zeros((n)) + 1j * np.sqrt(np.log(2.0)) * shape).real
    profile *= amplitude * wofz(2*np.sqrt(np.log(2.0)) * (x - position)/fwhm \
        + 1j * np.sqrt(np.log(2.0))*shape).real
    return profile



class ProfileFittingModel(BaseSpectralModel):

    _profiles = {
        "gaussian": (_gaussian, ("mean", "sigma", "amplitude")),
        "lorentzian": (_lorentzian, ("mean", "width", "amplitude")),
        "voigt": (_voigt, ("mean", "fwhm", "amplitude", "shape"))
    }

    def __init__(self, session, transitions, **kwargs):
        """
        Initialize a class for modelling spectra with analytic profile.

        :param session:
            The session that this spectral model will be associated with.

        :param transitions:
            A linelist containing atomic data for this model.
        """

        super(ProfileFittingModel, self).__init__(session, transitions,
            **kwargs)

        # Initialize metadata with default fitting values.
		# E. Holmbeck changed two lines. sigma from 0.5 to 0.3, 
		# and pixels from 3 to 5
                # APJ reverted, 0.3 is bad for MIKE
        self.metadata.update({
            "profile": "gaussian",
            "central_weighting": True,
            "window": 5,
            "continuum_order": 1,
            "detection_sigma": 0.5,
            "detection_pixels": 3,
            "max_iterations": 5,
            "wavelength_tolerance": 0.5,
            "velocity_tolerance": None,
            "mask": [],
            "antimask_flag": False,
            "elements": [self._verify_elements()],
            "species": [self._verify_species()]
        })

        # Set the model parameter names based on the current metadata.
        self._update_parameter_names()
        self._verify_transitions()
        self._verify_metadata()
                
        # Create a _repr_element for this.
        self._repr_element = self.transitions["element"][0]

        return None


    @property
    def abundance_uncertainties(self):
        try:
            result = self.metadata["fitted_result"][2]
            return result["abundance_uncertainties"][0]
        except KeyError:
            return None

    @property
    def expot(self):
        return self.transitions[0]["expot"]
    
    @property
    def loggf(self):
        return self.transitions[0]["loggf"]

    @property
    def equivalent_width(self):
        try:
            result = self.metadata["fitted_result"][2]
            equivalent_width = result["equivalent_width"][0]
        except KeyError:
            return None
        return 1000. * equivalent_width

    @property
    def equivalent_width_uncertainty(self):
        try:
            result = self.metadata["fitted_result"][2]
            err = 1000.*np.nanmax(np.abs(result["equivalent_width"][1:3]))
            return err
        except:
            return None
    
    @property
    def reduced_equivalent_width(self):
        eqw = self.equivalent_width
        if eqw is None: return None
        return np.log10(eqw/self.wavelength) - 3.

    @property
    def measurement_type(self):
        return "eqw"


    @property
    def fwhm(self):
        try:
            popt = self.metadata["fitted_result"][0]
            if self.metadata["profile"] == "gaussian":
                return popt["sigma"]*2.355
            elif self.metadata["profile"] == "lorentzian":
                return popt["width"] # I may be wrong about this value
            elif self.metadata["profile"] == "voigt":
                return popt["fwhm"] # I may be wrong about this value
            else:
                return None
        except (KeyError, TypeError):
            return None

    _profiles = {
        "gaussian": (_gaussian, ("mean", "sigma", "amplitude")),
        "lorentzian": (_lorentzian, ("mean", "width", "amplitude")),
        "voigt": (_voigt, ("mean", "fwhm", "amplitude", "shape"))
    }

    def _verify_transitions(self):
        """
        Verify that the atomic or molecular transitions associated with this
        class are valid.
        """

        # Check format first.
        transitions = self.transitions
        super(ProfileFittingModel, self)._verify_transitions()
        if len(transitions) > 1 and not isinstance(transitions, Row):
            raise ValueError("only a single transition can be associated with "
                             "a ProfileFittingModel")

        # Check that the transition does not have multiple element names.
        try:
            elem2 = transitions["elem2"][0]
        except IndexError:
            elem2 = transitions["elem2"]
        if elem2 != "":
            raise ValueError("only an atomic transition can be associated with "
                             "a ProfileFittingModel")
        return True

    def _verify_elements(self):
        """
        Return the element that will be measured by this model.
        """
        return self.transitions["element"][0].split()[0]
        
    def _verify_species(self):
        """
        Return the species that will be measured by this model.
        Ignore isotopes.
        """
        return np.floor(self.transitions["species"][0]*10)/10

    def _verify_metadata(self):
        """
        Verify the metadata associated with this class.
        """
        # TODO
        return True


    def _update_parameter_names(self):
        """
        Update the model parameter names based on the current metadata.
        """

        # Profile parameter names.
        func, parameter_names = self._profiles[self.metadata["profile"]]
        parameter_names = list(parameter_names)

        # Continuum coefficients.
        parameter_names += ["c{0}".format(i) \
            for i in range(self.metadata["continuum_order"] + 1)]

        # Update the bounds.
        bounds = {}
        if self.metadata["profile"] == "gaussian":
            bounds.update({
                "sigma": (-0.5, 0.5),
                "amplitude": (0, 1),
            })

        elif self.metadata["profile"] == "voigt":
            bounds["fwhm"] = (-0.5, 0.5)

        
        if self.metadata["wavelength_tolerance"] is not None \
        or self.metadata["velocity_tolerance"] is not None:

            # Convert velocity tolerance into wavelength.
            wavelength = self.transitions["wavelength"]
            try:
                wavelength = wavelength[0]
            except IndexError:
                None
            vt = abs(self.metadata.get("velocity_tolerance", None) or np.inf)
            wt = abs(self.metadata.get("wavelength_tolerance", None) or np.inf)

            bound = np.nanmin([
                wt, wavelength * vt/speed_of_light.to("km/s").value])

            bounds["mean"] = (wavelength - bound, wavelength + bound)

        else:
            # TODO: Allow the wavelength to be fixed.
            raise NotImplementedError("wavelength cannot be fixed yet; "
                                      "set a small tolerance on wavelength")

        self._parameter_names = parameter_names
        self._parameter_bounds = bounds

        return True


    def _initial_guess(self, spectrum, **kwargs):
        """
        Generate an initial guess for the model parameters.

        :param spectrum:
            The observed spectrum.
        """

        wavelength = self.transitions["wavelength"]
        try:
            wavelength = wavelength[0]
        except IndexError:
            None

        p0 = [
            wavelength,
            kwargs.pop("p0_sigma", 0.1),
        ]

        if spectrum is None:
            p0.append(0.5)
        else:
            p0.append(1.0 - \
                spectrum.flux[spectrum.dispersion.searchsorted(wavelength)])

        if self.metadata["profile"] == "voigt":
            p0.append(kwargs.pop("p0_shape", 0))

        # Continuum?
        if self.metadata["continuum_order"] > -1:
            p0.extend(([0] * self.metadata["continuum_order"]) + [1])

        return np.array(p0)


    def fit(self, spectrum=None, **kwargs):
        """
        Fit an absorption profile to the transition in the spectrum.

        :param spectrum: [optional]
            The observed spectrum to fit the profile transition model. If None
            is given, this will default to the normalized rest-frame spectrum in
            the parent session.
        """

        spectrum = self._verify_spectrum(spectrum)
        
        failure = False

        # Update internal metadata with any input parameters.
        # Ignore additional parameters because other BaseSpectralModels will
        # have different input arguments.
        for key in set(self.metadata).intersection(kwargs):
            self.metadata[key] = kwargs[key]

        # What model parameters are in the fitting process?
        # In synth these would be abundances, etc. Here they are profile/cont
        # parameters.
        self._update_parameter_names()
        
        # Get a bad initial guess.
        p0 = self._initial_guess(spectrum, **kwargs)

        # Build a mask based on the window fitting range, and prepare the data.
        mask = self.mask(spectrum)
        x, y = spectrum.dispersion[mask], spectrum.flux[mask]
        yerr, absolute_sigma = ((1.0/spectrum.ivar[mask])**0.5, True)
        if not np.all(np.isfinite(yerr)):
            yerr, absolute_sigma = (np.ones_like(x), False)

        # Central weighting?
        if self.metadata["central_weighting"]:
            yerr /= (1 + np.exp(-(x - p0[0])**2 / (4.0 * p0[1]**2)))

        # How many iterations to do?
        nearby_lines = []
        iterative_mask = np.isfinite(y * yerr)
        for iteration in range(self.metadata["max_iterations"]):
                
            if not any(iterative_mask):
                self.metadata["is_acceptable"] = False
                try:
                    del self.metadata["fitted_result"]
                except KeyError:
                    None

                return failure

            try:
                p_opt, p_cov = op.curve_fit(self.fitting_function, 
                    xdata=x[iterative_mask],
                    ydata=y[iterative_mask],
                    sigma=yerr[iterative_mask],
                    p0=p0, absolute_sigma=absolute_sigma)

            except:
                logger.exception(
                    "Exception raised in fitting atomic transition {0} "\
                    "on iteration {1}".format(self, iteration))

                if iteration == 0:
                    self.metadata["is_acceptable"] = False
                    try:
                        del self.metadata["fitted_result"]
                    except KeyError:
                        None

                    return failure

            # Look for outliers peaks.
            # TODO: use continuum or model?
            O = self.metadata["continuum_order"]
            continuum = np.ones_like(x) \
                if 0 > O else np.polyval(p_opt[-(O + 1):], x)    

            model = self(x, *p_opt)
            
            sigmas = (y - model)/np.std(y[iterative_mask])
            sigmas[~iterative_mask] = 0 # Ignore points that are already masked
            outliers = np.where(sigmas < -self.metadata["detection_sigma"])[0]

            # Look for groups of neighbouring outlier points.
            separators = np.repeat(1 + np.where(np.diff(outliers) > 1)[0], 2)
            separators = np.hstack([0, separators, None]).reshape(-1, 2)

            for start, end in separators:
                indices = outliers[start:end]

                # Require a minimum group size by number of pixels.
                if indices.size < self.metadata["detection_pixels"]: continue

                # Try and fit an absorption function to the centroid of the
                # region.

                lower_group_wl = x[indices[0]]
                upper_group_wl = x[indices[-1]]

                def model_nearby_line(x_, *p):
                    # Strict requirements on these parameters, since we don't
                    # need to know them precisely.
                    if not (lower_group_wl <= p[0] <= upper_group_wl) \
                    or abs(p[1]) > p_opt[1] \
                    or not (1 >= p[2] > 0):
                        return np.nan * np.ones_like(x_)
                    
                    return model[iterative_mask] * self(x_, *p)

                # Initial parameters for this line.
                p0_outlier = [
                    np.mean(x[indices]),
                    0.5 * p_opt[1],
                    np.max(continuum[indices] - y[indices])
                ]
                if self.metadata["profile"] == "voigt":
                    p0_outlier.append(p0[3])

                try:
                    p_out, cov = op.curve_fit(model_nearby_line,
                            xdata=x[iterative_mask],
                            ydata=y[iterative_mask],
                            sigma=yerr[iterative_mask],
                            p0=p0_outlier, absolute_sigma=absolute_sigma,
                            check_finite=True)

                except:
                    # Just take a narrow range and exclude that?
                    # TODO: The old SMH just did nothing in this scenario, but
                    #       we may want to revise that behaviour.
                    None

                else:

                    # Exclude +/- 3 sigma of the fitted line
                    l, u = (p_out[0] - 3 * p_out[1], p_out[0] + 3 * p_out[1])
                    
                    # Store this line and the region masked out by it.
                    nearby_lines.append([p_out, (u, l)])

                    # Now update the iterative mask to exclude this line.                    
                    iterative_mask *= ~((u > x) * (x > l))

            # Update p0 with the best guess from this iteration.
            p0 = p_opt.copy()

        # Finished looking for neighbouring lines
        # `max_iterations` rounds of removing nearby lines. Now do final fit:
        p_opt, p_cov = op.curve_fit(self.fitting_function, 
            xdata=x[iterative_mask],
            ydata=y[iterative_mask],
            sigma=yerr[iterative_mask],
            p0=p0, absolute_sigma=absolute_sigma)

        assert p_cov is not None

        # Make many draws from the covariance matrix.
        draws = kwargs.pop("covariance_draws", 
            self.session.setting("covariance_draws",100))
        percentiles = kwargs.pop("percentiles", \
            self.session.setting("error_percentiles",(16, 84)))
        if np.all(np.isfinite(p_cov)):
            p_alt = np.random.multivariate_normal(p_opt, p_cov, size=draws)
        else:
            p_alt = np.nan * np.ones((draws, p_opt.size))

        # Integrate the profile.
        profile, _ = self._profiles[self.metadata["profile"]]
        if profile == _gaussian:
            ew = abs(p_opt[1] * p_opt[2] * np.sqrt(2 * np.pi))
            ew_alt = np.abs(p_alt[:, 1] * p_alt[:, 2] * np.sqrt(2 * np.pi))
            ew_uncertainty = np.percentile(ew_alt, percentiles) - ew

        else:
            N, integrate_sigma = (len(_), kwargs.pop("integrate_sigma", 10))
            l, u = (
                p_opt[0] - integrate_sigma * p_opt[1],
                p_opt[0] + integrate_sigma * p_opt[1]
            )

            ew = np.abs(integrate.quad(profile, l, u, args=tuple(p_opt[:N]))[0])
            ew_alt = np.abs([integrate.quad(profile, l, u, args=tuple(_[:N]))[0] \
                for _ in p_alt])
            ew_uncertainty = np.percentile(ew_alt, percentiles) - ew
        
        # Calculate chi-square for the points that we modelled.
        ivar = spectrum.ivar[mask]
        if not np.any(np.isfinite(ivar)): ivar = 1
        residuals = y - self(x, *p_opt)
        residuals[~iterative_mask] = np.nan
        chi_sq = residuals**2 * ivar

        dof = np.isfinite(chi_sq).sum() - len(p_opt) - 1
        chi_sq = np.nansum(chi_sq)
        

        model_y = self(x, *p_opt)
        model_yerr = np.percentile(
            [self(x, *_) for _ in p_alt], percentiles, axis=0) - model_y
        model_yerr = np.max(np.abs(model_yerr), axis=0)

        ### DEBUG PLOT
        ##fig, ax = plt.subplots()
        ##ax.plot(x, y, c='k', drawstyle='steps-mid')
        ##
        ##O = self.metadata["continuum_order"]
        ##bg = np.ones_like(x) if 0 > O else np.polyval(p_opt[-(O + 1):], x)    
        ##for p, (u, l) in nearby_lines:
        ##    bg *= self(x, *p)
        ##
        ##    m = (u >= x) * (x >= l)
        ##    ax.scatter(x[m], y[m], facecolor="r")
        ##
        ##ax.plot(x, bg, c='r')
        ##ax.plot(x, model_y, c='b')
        ##
        ##model_err = np.percentile(
        ##    [self(x, *_) for _ in p_alt], percentiles, axis=0)
        ##
        ##ax.fill_between(x, model_err[0] + model_y, model_err[1] + model_y,
        ##    edgecolor="None", facecolor="b", alpha=0.5)
        
        # Convert x, model_y, etc back to real-spectrum indices.
        if self.session.setting("show_full_profiles", False):
            # HACK #152
            indices = spectrum.dispersion.searchsorted(x)
            x = spectrum.dispersion[indices[0]:1 + indices[-1]]
            y = spectrum.flux[indices[0]:1 + indices[-1]]
            model_y = self(x, *p_opt)
            model_yerr = np.percentile(
                [self(x, *_) for _ in p_alt], percentiles, axis=0) - model_y
            model_yerr = np.max(np.abs(model_yerr), axis=0)
            residuals = y - model_y
        else:
            x, model_y, model_yerr, residuals = self._fill_masked_arrays(
                spectrum, x, model_y, model_yerr, residuals)
        
        # We ignore the uncertainty in wavelength position because it only
        # affects the uncertainty in REW at the ~10^-5 level.
        rew = np.log10(ew/p_opt[0])
        rew_uncertainty = np.log10((ew + ew_uncertainty)/p_opt[0]) - rew

        fitting_metadata = {
            "equivalent_width": (ew, ew_uncertainty[0], ew_uncertainty[1]),
            "reduced_equivalent_width": np.hstack([rew, rew_uncertainty]),
            "data_indices": np.where(mask)[0][iterative_mask],
            "model_x": x,
            "model_y": model_y,
            "model_yerr": model_yerr,
            "residual": residuals,
            "nearby_lines": nearby_lines,
            "chi_sq": chi_sq,
            "dof": dof
        }

        # Update the equivalent width in the transition.
        # REMOVED: see Issue #38
        #self.transitions["equivalent_width"] = ew

        # Convert p_opt to ordered dictionary
        named_p_opt = OrderedDict(zip(self.parameter_names, p_opt))
        self.metadata["fitted_result"] = (named_p_opt, p_cov, fitting_metadata)

        # Only mark as acceptable if the model meets the quality constraints.
        self.is_acceptable = self.meets_quality_constraints_in_parent_session

        # Used normal fit
        self.metadata["used_monte_carlo_fit"] = False
        return self.metadata["fitted_result"]


    def montecarlo_fit(self, N=100, **kwargs):
        """
        Run a fit N times using spectrum ivar and report the resulting median and scatter
        """
        # Create base spectrum to add noise to
        spectrum = self._verify_spectrum(None)
        window = abs(self.metadata["window"])
        mask = (spectrum.dispersion >= self.transitions["wavelength"][0] - window) & \
               (spectrum.dispersion <= self.transitions["wavelength"][-1]+ window)
        spectrum = Spectrum1D(spectrum.dispersion[mask], spectrum.flux[mask], spectrum.ivar[mask])

        # Get rid of some kwargs
        orig_covariance_draws = kwargs.pop("covariance_draws", 
            self.session.setting("covariance_draws",100))

        # Fit multiple times
        EQWs = np.zeros(N) + np.nan
        for i in range(N):
            spec_i = spectrum.add_noise()
            try:
                meta = self.fit(spec_i, covariance_draws=1, **kwargs)
            except Exception as e:
                print(e)
                EQWs[i] = np.nan
            else:
                if meta:
                    EQWs[i] = meta[2]["equivalent_width"][0]
                else:
                    EQWs[i] = np.nan
        
        if np.any(np.isnan(EQWs)):
            logger.debug("{:.1f} {:.1f}: {} EQWs are nan".format(self.species[0], self.wavelength,
                                                             np.sum(np.isnan(EQWs))))
        
        # Run a fit to reset the plotting part to the original spectrum
        while True:
            meta = self.fit(covariance_draws=orig_covariance_draws, **kwargs)
            if "fitted_result" in self.metadata: break
            logger.debug("ERROR: default fit did not work! Trying again...")
        # Update with the montecarlo errors
        self.metadata["used_monte_carlo_fit"] = True
        percentiles = kwargs.pop("percentiles", \
            self.session.setting("error_percentiles",(16, 84)))
        ew  = np.nanmedian(EQWs)
        ew_uncertainty = np.nanpercentile(EQWs, percentiles) - ew
        rew = np.log10(ew/self.wavelength)
        rew_uncertainty = np.log10((ew + ew_uncertainty)/self.wavelength) - rew
        self.metadata["fitted_result"][2]["equivalent_width"] = (ew, ew_uncertainty[0], ew_uncertainty[1])
        self.metadata["fitted_result"][2]["reduced_equivalent_width"] = np.hstack([rew, rew_uncertainty])
        return 1000.*ew, 1000.*np.max(np.abs(ew_uncertainty))

    def _find_abund(self, ew_err=None):
        eqw = self.equivalent_width
        if eqw is None: return None
        if ew_err is not None:
            transitions = LineList.vstack([self.transitions[0], self.transitions[0]])
            transitions["equivalent_width"] = [eqw, eqw + ew_err]
            failure_message = "Failed to compute abundances for {:.1f} {:.3f}, eqws={:.2f}, {:.2f}".format(
                self.species[0], self.wavelength, eqw, eqw+ew_err)
        else:
            transitions = self.transitions.copy()
            transitions["equivalent_width"] = [eqw]
            failure_message = "Failed to compute abundances for {:.1f} {:.3f}, eqw={:.2f}".format(
                self.species[0], self.wavelength, eqw)
            
        try:
            abundances = self.session.rt.abundance_cog(
                self.session.stellar_photosphere, transitions, twd=self.session.twd)
        except Exception as e:
            logger.warn(failure_message)
            return None
        else:
            if len(abundances)==2:
                return np.abs(abundances[1] - abundances[0])
            else:
                return abundances[0]
    def find_error(self, sigma=1, ew_err_scale=1.0):
        """
        Find the sigma-th abundance uncertainty corresponding to sigma-th deviation in EQW
        Scale the equivalent width error by ew_err_scale
        """
        eqw = self.equivalent_width
        ew_err = self.equivalent_width_uncertainty
        if eqw is None or ew_err is None: return None
        ew_err *= ew_err_scale*sigma
        abund_err = self._find_abund(ew_err)
        self.metadata["{}_sigma_abundance_error".format(sigma)] = abund_err
        self.metadata["{}_sigma_eqw_error".format(sigma)] = ew_err
        return abund_err
        
    def propagate_stellar_parameter_error(self):
        e_Teff, e_logg, e_vt, e_MH = self.session.stellar_parameters_err
        Teff, logg, vt, MH = self.session.stellar_parameters
        alpha = self.session.metadata["stellar_parameters"]["alpha"]
        try:
            self.session.set_stellar_parameters(
                Teff, logg, vt, MH, alpha)
            abund0 = self._find_abund()
            self.session.set_stellar_parameters(
                Teff+e_Teff, logg, vt, MH, alpha)
            abund1 = self._find_abund()
            self.session.set_stellar_parameters(
                Teff, logg+e_logg, vt, MH, alpha)
            abund2 = self._find_abund()
            self.session.set_stellar_parameters(
                Teff, logg, vt+e_vt, MH, alpha)
            abund3 = self._find_abund()
            self.session.set_stellar_parameters(
                Teff, logg, vt, MH+e_MH, alpha)
            abund4 = self._find_abund()
            dTeff_error = abund1-abund0
            dlogg_error = abund2-abund0
            dvt_error = abund3-abund0
            dMH_error = abund4-abund0
        except:
            self.session.set_stellar_parameters(
                Teff, logg, vt, MH, alpha)
        else:
            self.session.set_stellar_parameters(
                Teff, logg, vt, MH, alpha)
            self.metadata["systematic_abundance_error"] = np.sqrt(
                dTeff_error**2 + dlogg_error**2 + dvt_error**2 + dMH_error**2)
            self.metadata["systematic_stellar_parameter_abundance_error"] = {
                "effective_temperature": dTeff_error,
                "surface_gravity": dlogg_error,
                "microturbulence": dvt_error,
                "metallicity": dMH_error
            }
        return self.metadata["systematic_abundance_error"]
        

    def __call__(self, dispersion, *parameters):
        """
        Generate data at the dispersion points, given the parameters.

        :param dispersion:
            An array of dispersion points to calculate the data for.

        :param parameters:
            Keyword arguments of the model parameters and their values.
        """

        function, profile_parameters = self._profiles[self.metadata["profile"]]

        N = len(profile_parameters)
        y = 1.0 - function(dispersion, *parameters[:N])

        # Assume rest of the parameters are continuum coefficients.
        if parameters[N:]:
            y *= np.polyval(parameters[N:], dispersion)
        
        return y


if __name__ == "__main__":


    # Load some atomic transitions
    # Measure them in a spectrum.
    # Get abundances from MOOGSILENT.


    import smh
    a = smh.Session([
        "/Users/arc/codes/smh/hd44007red_multi.fits",
        "/Users/arc/codes/smh/hd44007blue_multi.fits",
    ])

    spec = smh.specutils.Spectrum1D.read("../hd140283.fits")


    foo = ProfileFittingModel({"wavelength": [5202.336]}, a)

    raise a


    """
    with open("a", "r") as fp:
        b = fp.readlines()

    for each in b:
        

        foo = ProfileFittingModel({"wavelength": [float(each.split()[0])]}, a)

        try:
            foo.fit(spec)
        except:
            print("failed on {}".format(each))

    """


