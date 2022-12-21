
#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model for fitting synthetic spectra to spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["SpectralSynthesisModel"]

import logging
import itertools
import numpy as np
import scipy.optimize as op
import scipy.interpolate
from collections import OrderedDict
from six import string_types, iteritems
from scipy.ndimage import gaussian_filter
from scipy import stats

from .base import BaseSpectralModel, penalized_curve_fit_lm
from smh import utils
from smh.specutils import Spectrum1D
from smh.photospheres.abundances import asplund_2009 as solar_composition
# HACK TODO REVISIT SOLAR SCALE: Read from session defaults?
from astropy.table import Table

logger = logging.getLogger(__name__)


def _fix_rt_abundances(rt_abundances):
    rt_abundances = rt_abundances.copy()
    for key in rt_abundances:
        if np.isnan(rt_abundances[key]):
            rt_abundances[key] = -9.0
    return rt_abundances
    
def approximate_spectral_synthesis(model, centroids, bounds, rt_abundances={},
                                   isotopes={}):

    # Generate spectra at +/- the initial bounds. For all elements.
    species \
        = [utils.element_to_species(e) for e in model.metadata["elements"]]
    N = len(species)
    try:
        centroids[0]
    except IndexError:
        centroids = centroids * np.ones(N)

    calculated_abundances = np.array(list(
        itertools.product(*[[c-bounds, c, c+bounds] for c in centroids])))

    # Explicitly specified abundances
    # np.nan goes to -9.0
    rt_abundances = _fix_rt_abundances(rt_abundances)

    # Group syntheses into 5 where possible.
    M, K = 5, calculated_abundances.shape[0]
    fluxes = None
    for i in range(1 + int(K / M)):
        abundances = {}
        for j, specie in enumerate(species):
            abundances[specie] = calculated_abundances[M*i:M*(i + 1), j]

        # Include explicitly specified abundances.
        abundances.update(rt_abundances)
        #logger.debug(abundances)

        spectra = model.session.rt.synthesize(
            model.session.stellar_photosphere, 
            model.transitions, abundances=abundances,
            isotopes=isotopes,
            twd=model.session.twd) # TODO other kwargs?

        dispersion = spectra[0][0]
        if fluxes is None:
            fluxes = np.nan * np.ones((K, dispersion.size))

        for j, spectrum in enumerate(spectra):
            fluxes[M*i + j, :] = spectrum[1]

    def call(*parameters):
        #print(*parameters)
        N = len(model.metadata["elements"])
        if len(parameters) < N:
            raise ValueError("missing parameters")

        flux = scipy.interpolate.griddata(calculated_abundances, fluxes, 
            np.array(parameters[:N]).reshape(-1, N)).flatten()
        return (dispersion, flux)

    # Make sure the function works, otherwise fail early so we can debug.
    assert np.isfinite(call(*centroids)[1]).all()

    return call


class SpectralSynthesisModel(BaseSpectralModel):

    def __init__(self, session, transitions, elements,
                 what_species=None, what_wavelength=None, what_expot=None, what_loggf=None,
                 **kwargs):
        """
        Initialize a class for modelling spectra with synthesis.

        :param session:
            The session that this spectral model will be associated with.

        :param transitions:
            A linelist containing atomic data for this model.
        
        :param elements:
            The element(s) to be measured from the data.

        :param what_species:
            Specify this synthesis is associated with only some species of its elements
        :param what_wavelength:
            Specify this synthesis should be labeled with specific wavelength
        :param what_expot:
            Specify this synthesis should be labeled with specific expot
        :param what_loggf:
            Specify this synthesis should be labeled with specific loggf
        """
        
        super(SpectralSynthesisModel, self).__init__(session, transitions,
            **kwargs)

        # Initialize metadata with default fitting values.
        self.metadata.update({
            "mask": [],
            "window": 0, 
            "continuum_order": 1,
            "velocity_tolerance": 5,
            "smoothing_kernel": True,
            "initial_abundance_bounds": 1,
            "elements": self._verify_elements(elements),
            "species": self._verify_species(elements, what_species),
            "rt_abundances": {},
            "manual_continuum": 1.0,
            "manual_sigma_smooth":0.15,
            "manual_rv":0.0
        })

        # Set the model parameter names based on the current metadata.
        self._update_parameter_names()
        self._verify_transitions()

        # Set rt_abundances to have all the elements with nan
        unique_elements = np.unique(np.array(self.transitions["elem1"]))
        unique_elements = np.concatenate([unique_elements,np.array(np.unique(self.transitions["elem2"]))])
        #logger.debug(unique_elements)
        #logger.debug(type(unique_elements))
        unique_elements = np.unique(unique_elements)
        
        rt_abundances = {}
        try:
            unique_elements = [x.decode() for x in unique_elements]
        except:
            pass
        for elem in unique_elements:
            if elem in ["","H"]: continue
            if elem in self.elements: continue
            assert elem in utils.periodic_table, elem
            rt_abundances[elem] = np.nan
        self.metadata.update({"rt_abundances": rt_abundances})

        if len(self.elements) == 1:
            # Which of these lines is in the line list?
            matches = (self.transitions["elem1"] == self.elements[0])

            if not np.any(matches):
                self._repr_element = ", ".join(self.elements)
            else:
                unique_species = np.unique(self.transitions["element"][matches])
                if len(unique_species) == 1:
                    self._repr_element = unique_species[0]
                else:
                    self._repr_element = ", ".join(self.elements)
        else:
            # Many elements fitted simultaneously.
            self._repr_element = ", ".join(self.elements)

        if "Fe" not in elements:
            self.metadata["use_for_stellar_parameter_inference"] = False

        ## Set some display variables
        if what_wavelength is not None:
            self._wavelength = what_wavelength
        if what_expot is None: what_expot = np.nan
        if what_loggf is None: what_loggf = np.nan
        self._expot = what_expot
        self._loggf = what_loggf

        return None

    @property
    def expot(self):
        ## TODO for most syntheses this is well-defined
        return self._expot
    
    @property
    def loggf(self):
        ## TODO for most syntheses the combined loggf is well-defined
        return self._loggf

    @property
    def measurement_type(self):
        return "syn"

    @property
    def fwhm(self):
        try:
            return self.metadata["manual_sigma_smooth"]*2.355
        except (IndexError, ValueError): # Otherwise, use manual value
            return None

    def _verify_elements(self, elements):
        """
        Verify that the atomic or molecular transitions associated with this
        class are valid.

        :param elements:
            The element(s) (string or list-type of strings) that will be
            measured in this model.
        """

        # Format the elements and then check that all are real.
        if isinstance(elements, string_types):
            elements = [elements]

        elements = [str(element).title() for element in elements]
        transitions = self.transitions
        for element in elements:
            # Is the element real?
            if element not in utils.periodic_table:
                raise ValueError("element '{}' is not valid".format(element))

            # Is it in the transition list?
            if  element not in transitions["elem1"] \
            and element not in transitions["elem2"]:
                raise ValueError(
                    "element '{0}' does not appear in the transition list"\
                        .format(element))

        return elements


    def _verify_species(self, elements, what_species=None):
        # Format the elements and then check that all are real.
        if isinstance(elements, string_types):
            elements = [elements]

        elements = [str(element).title() for element in elements]
        species = []
        transitions = self.transitions
        for element in elements:
            # Get the species associated with this element
            ii = np.logical_or(
                transitions["elem1"] == element,
                transitions["elem2"] == element)

            # Note plurality/singularity of specie/species.
            # APJ modified to remove isotopes in species
            specie = transitions[ii]["species"]
            specie = (specie*10).astype(int)/10.0
            specie = list(np.unique(specie))
            if what_species is not None:
                for s in specie[::-1]: # go backwards to remove properly
                    if s not in what_species: specie.remove(s)
            species.append(specie)

        return species
        

    def _initial_guess(self, spectrum, **kwargs):
        """
        Return an initial guess about the model parameters.

        :param spectrum:
            The observed spectrum.
        """

        # Potential parameters:
        # elemental abundances, continuum coefficients, smoothing kernel,
        # velocity offset

        defaults = {
            "sigma_smooth": self.metadata["manual_sigma_smooth"], #0.1,
            "vrad": self.metadata["manual_rv"]
        }
        p0 = []
        for parameter in self.parameter_names:
            default = defaults.get(parameter, None)
            if default is not None:
                p0.append(default)

            elif parameter.startswith("log_eps"):
                # The elemental abundances are in log_epsilon format.
                # If the value was already fit, we use that as the initial guess
                # Otherwise, we assume a scaled-solar initial value based on the stellar [M/H]

                # Elemental abundance.
                element = parameter.split("(")[1].rstrip(")")

                # Assume scaled-solar composition.
                default_value = solar_composition(element) + \
                        self.session.metadata["stellar_parameters"]["metallicity"]
                try:
                    fitted_result = self.metadata["fitted_result"]
                except KeyError:
                    p0.append(default_value)
                else:
                    value = fitted_result[0][parameter]
                    if np.isnan(value):
                        p0.append(default_value)
                    else:
                        p0.append(value)
                
            elif parameter.startswith("c"):
                # Continuum coefficient.
                p0.append((0, 1)[parameter == "c0"])

            else:
                raise ParallelUniverse("this should never happen")

        return np.array(p0)


    def _update_parameter_names(self):
        """
        Update the model parameter names based on the current metadata.
        """

        bounds = {}
        parameter_names = []

        # Abundances of different elements to fit simultaneously.
        parameter_names.extend(
            ["log_eps({0})".format(e) for e in self.metadata["elements"]])

        # Radial velocity?
        vt = abs(self.metadata["velocity_tolerance"] or 0)
        if vt > 0:
            rv_mid = self.metadata["manual_rv"]
            parameter_names.append("vrad")
            bounds["vrad"] = [rv_mid-vt, rv_mid+vt]

        # Continuum coefficients?
        parameter_names += ["c{0}".format(i) \
            for i in range(self.metadata["continuum_order"] + 1)]

        # Gaussian smoothing kernel?
        if self.metadata["smoothing_kernel"]:
            # TODO: Better init of this
            smooth_mid = self.metadata["manual_sigma_smooth"]
            bounds["sigma_smooth"] = (-5, +5)
            parameter_names.append("sigma_smooth")

        self._parameter_bounds = bounds
        self._parameter_names = parameter_names
        return True


    def fit(self, spectrum=None, penalty_function=None, **kwargs):
        """
        Fit a synthesised model spectrum to the observed spectrum.

        :param spectrum: [optional]
            The observed spectrum to fit the synthesis spectral model. If None
            is given, this will default to the normalized rest-frame spectrum in
            the parent session.

        :param penalty_function: [optional]
            A function of the form penalty_function(params) [not *params]
            to penalize the parameters. Can be used to make a prior on some parameter.
        """

        # Check the observed spectrum for validity.
        spectrum = self._verify_spectrum(spectrum)

        # Update internal metadata with any input parameters.
        # Ignore additional parameters because other BaseSpectralModels will
        # have different input arguments.
        for key in set(self.metadata).intersection(kwargs):
            self.metadata[key] = kwargs.pop(key)

        # Update the parameter names in case they have been updated due to the
        # input kwargs.
        self._update_parameter_names()
        
        # Build a mask based on the window fitting range, and prepare the data.
        mask = self.mask(spectrum)
        x, y = spectrum.dispersion[mask], spectrum.flux[mask]
        yerr, absolute_sigma = ((1.0/spectrum.ivar[mask])**0.5, True)
        if not np.all(np.isfinite(yerr)):
            yerr, absolute_sigma = (np.ones_like(x), False)

        # Get a bad initial guess.
        p0 = self._initial_guess(spectrum, **kwargs)

        cov, bounds = None, self.metadata["initial_abundance_bounds"]

        # We allow the user to specify their own function to decide if the
        # approximation is good.
        tolerance_metric = kwargs.pop("tolerance_metric",
            lambda t, a, y, yerr, abs_sigma: np.nansum(np.abs(t - a)) < 0.05)

        # Here we set a hard bound limit where linear interpolation must be OK.
        while bounds > 0.01: # dex
            central = p0[:len(self.metadata["elements"])]
            approximater \
                = approximate_spectral_synthesis(self, central, bounds, 
                                                 rt_abundances=self.metadata["rt_abundances"],
                                                 isotopes=self.session.metadata["isotopes"])

            def objective_function(x, *parameters):
                synth_dispersion, intensities = approximater(*parameters)
                return self._nuisance_methods(
                    x, synth_dispersion, intensities, *parameters)
                
            if penalty_function is None:
                p_opt, p_cov = op.curve_fit(objective_function, xdata=x, ydata=y,
                        sigma=yerr, p0=p0, absolute_sigma=absolute_sigma)
            else:
                p_opt, p_cov = penalized_curve_fit_lm(
                    objective_function, xdata=x, ydata=y, penalty_function=penalty_function,
                    sigma=yerr, p0=p0, absolute_sigma=absolute_sigma)

            # At small bounds it can be difficult to estimate the Jacobian.
            # So if it is correctly approximated once then we keep that in case
            # things screw up later, since it provides a conservative estimate.
            if cov is None or np.isfinite(p_cov).any():
                cov = p_cov.copy()

            # Is the approximation good enough?
            model_y = self(x, *p_opt)
            close_enough = tolerance_metric(
                model_y,                        # Synthesised
                objective_function(x, *p_opt),  # Approxsised
                y, yerr, absolute_sigma)

            if close_enough: break
            else:
                bounds /= 2.0
                p0 = p_opt.copy()

                logger.debug("Dropping bounds to {} because approximation was "
                             "not good enough".format(bounds))

        # Use (p_opt, cov) not (p_opt, p_cov)

        # Make many draws from the covariance matrix.
        draws = kwargs.pop("covariance_draws", 
            self.session.setting("covariance_draws",100))
        percentiles = kwargs.pop("percentiles", \
            self.session.setting("error_percentiles",(16, 84)))
        if np.all(np.isfinite(cov)):
            p_alt = np.random.multivariate_normal(p_opt, cov, size=draws)
            model_yerr = model_y - np.percentile(
                [objective_function(x, *_) for _ in p_alt], percentiles, axis=0)
        else:
            p_alt = np.nan * np.ones((draws, p_opt.size))
            model_yerr = np.nan * np.ones((2, x.size))
        model_yerr = np.max(np.abs(model_yerr), axis=0)
        
        # Calculate chi-square for the points that we modelled.
        ivar = spectrum.ivar[mask]
        if not np.any(np.isfinite(ivar)): ivar = 1
        residuals = y - model_y
        chi_sq = residuals**2 * ivar

        dof = np.isfinite(chi_sq).sum() - len(p_opt) - 1
        chi_sq = np.nansum(chi_sq)

        x, model_y, model_yerr, residuals = self._fill_masked_arrays(
            spectrum, x, model_y, model_yerr, residuals)

        # Convert result to ordered dict.
        named_p_opt = OrderedDict(zip(self.parameter_names, p_opt))
        
        # Save plotting
        if self.session.setting("full_synth_resolution",False):
            plot_x = self.metadata["raw_synthetic_dispersion"]
            plot_y = self._nuisance_methods(plot_x, plot_x,
                                            self.metadata["raw_synthetic_flux"],
                                            *named_p_opt.values())
        else:
            plot_x = x
            plot_y = model_y

        fitting_metadata = {
            "model_x": x,
            "model_y": model_y,
            "model_yerr": model_yerr,
            "residual": residuals,
            "plot_x": plot_x,
            "plot_y": plot_y,
            "chi_sq": chi_sq,
            "dof": dof,
            "abundances": p_opt[:len(self.elements)]
        }

        self.metadata["fitted_result"] = (named_p_opt, cov, fitting_metadata)
        self.metadata["is_acceptable"] = True
        if "vrad" in self.parameter_names:
            self.metadata["manual_rv"] = named_p_opt["vrad"]
        if "sigma_smooth" in self.parameter_names:
            self.metadata["manual_sigma_smooth"] = named_p_opt["sigma_smooth"]
        
        return self.metadata["fitted_result"]


    def iterfit(self, maxiter=10, tol=.01, init_abund = np.nan, **kwargs):
        """
        Iteratively run fit() until abundance is converged.
        Only works for syntheses with single abundance.
        Uses current fitting parameter settings.
        """
        ## Only do this with a single value right now
        assert len(self.elements) == 1, self.elements

        abund = init_abund
        for i in range(maxiter):
            lastabund = abund
            self.fit(**kwargs)
            abund = self.abundances[0]
            if np.abs(abund - lastabund) < tol: break
        else:
            logger.warn("iterfit: Reached {}/{} iter without converging. Now {:.3f}, last {:.3f}, tol={:.3f}".format(
                    i, maxiter, abund, lastabund, tol))
        return abund

    def export_fit(self, synth_fname, data_fname, parameters_fname):
        self._update_parameter_names()
        spectrum = self._verify_spectrum(None)
        try:
            (named_p_opt, cov, meta) = self.metadata["fitted_result"]
        except KeyError:
            logger.info("Please run a fit first!")
            return None
        ## Write synthetic spectrum
        # take the mean of the two errors for simplicity
        if len(meta["model_yerr"].shape) == 1:
            ivar = (meta["model_yerr"])**-2.
        elif len(meta["model_yerr"].shape) == 2:
            assert meta["model_yerr"].shape[0] == 2, meta["model_yerr"].shape
            ivar = (np.nanmean(meta["model_yerr"], axis=0))**-2.
        #synth_spec = Spectrum1D(meta["model_x"], meta["model_y"], ivar)
        try:
            plot_x = meta["plot_x"]
            plot_y = meta["plot_y"]
        except KeyError:
            plot_x = meta["model_x"]
            plot_y = meta["model_y"]
            if self.session.setting("full_synth_resolution", False):
                logger.warn("Note: exporting synthesis sampled at data points, not full synth resolution!")

        if len(meta["plot_x"]) != len(ivar):
            synth_spec = Spectrum1D(meta["plot_x"], meta["plot_y"], np.ones(len(meta["plot_x"]))*1e6)
        else:
            synth_spec = Spectrum1D(meta["plot_x"], meta["plot_y"], ivar)
        synth_spec.write(synth_fname)
        
        ## Write data only in the mask range
        mask = self.mask(spectrum)
        x, y, ivar = spectrum.dispersion[mask], spectrum.flux[mask], spectrum.ivar[mask]
        data_spec = Spectrum1D(x, y, ivar)
        data_spec.write(data_fname)
        
        ## Write out all fitting parameters
        with open(parameters_fname, "w") as fp:
            ## Fitted abundances
            fp.write("----------Fit Abundances----------\n")
            for elem, abund in zip(self.elements, meta["abundances"]):
                fp.write("{} {:.3f}\n".format(elem, abund))
            fp.write("\n")
            fp.write("----------Fixed Abundances----------\n")
            for elem, abund in iteritems(self.metadata["rt_abundances"]):
                fp.write("{} {:.3f}\n".format(elem, abund))
            fp.write("\n")
            fp.write("p_opt: {}\n".format(named_p_opt))
            fp.write("cov: {}\n".format(cov))
            fp.write("\n")
            fp.write("----------Parameters----------\n")
            for key, value in iteritems(self.metadata):
                if key=="rt_abundances": continue
                if key=="fitted_result": continue
                fp.write("{}: {}\n".format(key, value))
            for key, value in iteritems(meta):
                if key in ['chi_sq','dof']:
                    fp.write("{}: {}\n".format(key, value))
            fp.write("\n")
            fp.write("----------Isotopes----------\n")
            fp.write(str(self.session.metadata["isotopes"]))
            fp.write("\n")
        return None
    

    def export_plot_data(self, logeps_err, wlmin=None, wlmax=None,
                         normalize_data=False):
        """
        Exports data and synthesized models to be used for a synthesis plot.
        Returns two tables:
        (1) Model table: wl, flux(A(X)=-9), flux(fit-err), flux(fit), flux(fit+err)
        (2) Data table: wl, flux, err, and model residuals
        If normalize_data=True, applies the continuum fit to the data rather than the model
          (so normalized flux = 1 is the continuum)
        """
        assert len(self.elements) == 1, self.elements
        assert self.session.setting("full_synth_resolution",True)
        try:
            named_p, cov, meta = self.metadata["fitted_result"]
        except KeyError:
            logger.warn("No fit run yet!")
            return None
        elem = self.elements[0]
        assert "log_eps({})".format(elem) in named_p
        
        ## Assemble data
        self._update_parameter_names()
        spectrum = self._verify_spectrum(None)
        if (wlmin is None) or (wlmax is None):
            wavelengths = self.transitions["wavelength"]
            try:
                wlmin = wavelengths[0]
                wlmax = wavelengths[-1]
            except IndexError: # Single row.
                wlmin, wlmax = (wavelengths, wavelengths)
            window = abs(self.metadata["window"])
            wlmin -= window
            wlmax += window
        mask = (spectrum.dispersion >= wlmin) \
             * (spectrum.dispersion <= wlmax)
        data_disp = spectrum.dispersion[mask]
        data_flux = spectrum.flux[mask]
        data_errs = 1./np.sqrt(spectrum.ivar[mask])
        wlrange = (data_disp[0], data_disp[-1])
        
        model_output = [self.metadata["raw_synthetic_dispersion"]]
        data_output = [data_disp, data_flux, data_errs]

        # Run synths
        abund0 = self.abundances[0]
        abunds_to_synth = [-9, abund0 - logeps_err, abund0, abund0 + logeps_err]
        orig_p = named_p.copy()
        for abund in abunds_to_synth:
            abund = float(abund)
            named_p["log_eps({})".format(elem)] = abund
            self.update_fit_after_parameter_change(synthesize=True, wlrange=wlrange, with_mask=False)
            data_output.append(meta["residual"].copy())
            model_output.append(meta["plot_y"].copy())
            
        # Reset state
        named_p = orig_p
        self.update_fit_after_parameter_change(synthesize=True, wlrange=wlrange)
        
        data_output = Table(data_output, names=["wl","flux","err","r_none","r_-err","r_fit","r_+err"])
        model_output = Table(model_output, names=["wl","f_none","f_-err","f_fit","f_+err"])
        
        if normalize_data:
            ## remove the continuum from model and data
            modeldisp = model_output["wl"]
            datadisp  = data_output["wl"]
            parameters = self.metadata["fitted_result"][0].values()

            names = self.parameter_names
            O = self.metadata["continuum_order"]
            if 0 > O:
                continuum1 = 1 * self.metadata["manual_continuum"]
                continuum2 = 1 * self.metadata["manual_continuum"]
            else:
                continuum1 = np.polyval([parameters[names.index("c{}".format(i))] \
                    for i in range(O + 1)][::-1], modeldisp)
                continuum2 = np.polyval([parameters[names.index("c{}".format(i))] \
                    for i in range(O + 1)][::-1], datadisp)
            for col in ["f_none", "f_-err", "f_fit", "f_+err"]:
                model_output[col] = model_output[col]/continuum1
            for col in ["flux", "r_none", "r_-err", "r_fit", "r_+err"]:
                data_output[col] = data_output[col]/continuum2

        return elem, orig_p, logeps_err, model_output, data_output
        
    def update_fit_after_parameter_change(self, synthesize=True, wlrange=None, with_mask=True):
        """
        After manual parameter change, update fitting metadata 
        for plotting purposes.
        This somewhat messes up the state, but if you are trying to
        fiddle manually you probably don't care about those anyway.
        
        :param synthesize:
            If True (default), resynthesizes the spectrum.
            If False, use existing raw synthetic spectrum and reapply nuisance methods.
        """
        self._update_parameter_names()
        spectrum = self._verify_spectrum(None)
        try:
            (named_p_opt, cov, meta) = self.metadata["fitted_result"]
        except KeyError:
            logger.info("Please run a fit first!")
            return None
        if wlrange is None:
            model_x = meta["model_x"]
            data_y = meta["residual"] + meta["model_y"]
        else:
            if with_mask:
                mask = self.mask(spectrum)
            else:
                mask = (spectrum.dispersion >= wlrange[0]) & (spectrum.dispersion <= wlrange[-1])
            model_x = spectrum.dispersion[mask]
            data_y = spectrum.flux[mask]
        model_yerr = np.nan * np.ones((1, model_x.size))
        
        # recompute model
        if synthesize:
            model_y = self(model_x, *named_p_opt.values())
        else:
            try:
                model_y = self.metadata["raw_synthetic_flux"]
            except KeyError:
                logger.debug("No raw spectrum, synthesizing...")
                model_y = self(model_x, *named_p_opt.values())
            else:
                model_y = self._nuisance_methods(model_x,
                                                 self.metadata["raw_synthetic_dispersion"],
                                                 self.metadata["raw_synthetic_flux"],
                                                 *named_p_opt.values())
        residuals = data_y - model_y
        
        model_x, model_y, model_yerr, residuals = self._fill_masked_arrays(
            spectrum, model_x, model_y, model_yerr, residuals)
        # recompute plot
        if self.session.setting("full_synth_resolution",False):
            plot_x = self.metadata["raw_synthetic_dispersion"]
            plot_y = self._nuisance_methods(plot_x, plot_x,
                                            self.metadata["raw_synthetic_flux"],
                                            *named_p_opt.values())
        else:
            plot_x = model_x
            plot_y = model_y

        meta["model_y"] = model_y
        meta["model_yerr"] = model_yerr
        meta["residual"] = residuals
        meta["plot_x"] = plot_x
        meta["plot_y"] = plot_y
        
        return None

    def get_synth(self, abundances):
        try:
            (named_p_opt, cov, meta) = self.metadata["fitted_result"]
        except KeyError:
            logger.info("Please run a fit first!")
            return None
        abundances = _fix_rt_abundances(abundances)
        synth_dispersion, intensities, meta = self.session.rt.synthesize(
            self.session.stellar_photosphere, self.transitions,
            abundances=abundances, 
            isotopes=self.session.metadata["isotopes"],
            twd=self.session.twd)[0] # TODO: Other RT kwargs......
        return synth_dispersion, self._nuisance_methods(
            synth_dispersion, synth_dispersion, intensities, *named_p_opt.values())

    def __call__(self, dispersion, *parameters):
        """
        Generate data at the dispersion points, given the parameters.

        :param dispersion:
            An array of dispersion points to calculate the data for.

        :param *parameters:
            Keyword arguments of the model parameters and their values.
        """

        # Parse the elemental abundances, because they need to be passed to
        # the synthesis routine.
        abundances = {}
        names =  self.parameter_names
        for name, parameter in zip(names, parameters):
            if name.startswith("log_eps"):
                element = name.split("(")[1].rstrip(")")
                abundances[utils.element_to_species(element)] = parameter
            else: break # The log_eps abundances are always first.

        rt_abundances = _fix_rt_abundances(self.metadata["rt_abundances"])
        abundances.update(rt_abundances)

        # Produce a synthetic spectrum.
        synth_dispersion, intensities, meta = self.session.rt.synthesize(
            self.session.stellar_photosphere, self.transitions,
            abundances=abundances, 
            isotopes=self.session.metadata["isotopes"],
            twd=self.session.twd)[0] # TODO: Other RT kwargs......

        # Save raw synthetic spectrum
        self.metadata["raw_synthetic_dispersion"] = synth_dispersion
        self.metadata["raw_synthetic_flux"] = intensities

        return self._nuisance_methods(
            dispersion, synth_dispersion, intensities, *parameters)


    def _nuisance_methods(self, dispersion, synth_dispersion, intensities,
        *parameters):
        """
        Apply nuisance operations (convolution, continuum, radial velocity, etc)
        to a model spectrum.

        Attempts to use model parameters for continuum, smoothing, and radial velocity.
        If those are not there, uses self.metadata["manual_<x>"] where
        <x> = "continuum", "sigma_smooth", and "rv" respectively.

        :param dispersion:
            The dispersion points to evaluate the model at.

        :param synth_dispersion:
            The model dispersion points.

        :param intensities:
            The intensities at the model dispersion points.

        :param *parameters:
            The model parameter values.

        :returns:
            A convolved, redshifted model with multiplicatively-entered
            continuum.
        """

        # Continuum.
        names = self.parameter_names
        O = self.metadata["continuum_order"]
        if 0 > O:
            continuum = 1 * self.metadata["manual_continuum"]
        else:
            continuum = np.polyval([parameters[names.index("c{}".format(i))] \
                for i in range(O + 1)][::-1], synth_dispersion)

        model = intensities * continuum

        # Smoothing.
        try: # If in parameters to vary, use that
            sigma_smooth = parameters[names.index("sigma_smooth")]
        except (IndexError, ValueError): # Otherwise, use manual value
            sigma_smooth = self.metadata["manual_sigma_smooth"]
        # Scale value by pixel diff.
        kernel = abs(sigma_smooth)/np.mean(np.diff(synth_dispersion))
        if kernel > 0:
            model = gaussian_filter(model, kernel)

        try:
            v = parameters[names.index("vrad")]
        except ValueError:
            v = self.metadata["manual_rv"]

        #logger.debug("smooth: {:.4f}, rv: {:.4f}".format(sigma_smooth,v))
        # Interpolate the model spectrum onto the requested dispersion points.
        return np.interp(dispersion, synth_dispersion * (1 + v/299792458e-3), 
            model, left=1, right=1)


    def find_error(self, sigma=1, max_elem_diff=2.0, abundtol=.001, pix_per_element=1.0):
        """
        Increase abundance of elem until chi2 is <sigma> higher.
        Return the difference in abundance
        """
        assert len(self.metadata["elements"]) == 1, self.metadata["elements"]
        assert sigma > 0, sigma
        elem = self.metadata["elements"][0]
        self._update_parameter_names()

        # Load current fit
        try:
            (orig_p_opt, cov, meta) = self.metadata["fitted_result"]
        except KeyError:
            logger.info("Please run a fit first!")
            return None
        elem_p_opt_name = "log_eps({})".format(elem)
        assert elem_p_opt_name in orig_p_opt, (elem, orig_p_opt)
        abund0 = orig_p_opt[elem_p_opt_name]

        # Set up data again from scratch
        spectrum = self._verify_spectrum(None)
        mask = self.mask(spectrum)
        x = spectrum.dispersion[mask]
        data_y = spectrum.flux[mask]
        data_ivar = spectrum.ivar[mask]
        
        # recompute chi2 and get target
        model_y = self(x, *orig_p_opt.values())
        chi2_best = np.nansum((data_y - model_y)**2 * data_ivar)
        frac = stats.norm.cdf(sigma) - stats.norm.cdf(-sigma)
        target_chi2 = chi2_best + stats.chi2.ppf(frac, 1)*pix_per_element
        logger.debug("chi2={:.2f}, target chi2={:.2f}".format(chi2_best, target_chi2))
        
        # Find abundance where chi2 matches
        p_opt = orig_p_opt.copy() # temporary thing
        def minfn(abund):
            logger.debug("abund = {:.3f}".format(abund))
            # calculate model_y: run synthesis with new abund
            p_opt[elem_p_opt_name] = abund
            #try:
            model_y = self(x, *p_opt.values())
            #except Exception as e:
            #    print(e)
            #    #logger.debug(e)
            #    return np.nan
            
            # return difference in chi2
            _chi2 = np.nansum((data_y - model_y)**2 * data_ivar)
            logger.debug("chi2 diff = {:.3f}".format(_chi2 - target_chi2))
            return _chi2 - target_chi2
        
        abund1 = op.brentq(minfn, abund0, abund0+max_elem_diff, xtol = abundtol)
        logger.debug("orig abund={:.2f} error abund={:.2f}".format(abund0, abund1))
        # Reset the raw synthetic spectrum back to what it was and return
        _ = self(x, *orig_p_opt.values())

        self.metadata["{}_sigma_abundance_error".format(sigma)] = np.abs(abund1 - abund0)
        return np.abs(abund1 - abund0)
        

    def check_line_detection(self, sigma=3):
        """
        Compare the current fit to a synthetic fit with 10^-9 of the current abundance.
        Returns (is_detection, delta_chi2).
        Detection assumes a gaussian sigma of 3 (by default), or tail probability of 100-99.7%,
          assuming 1 degree of freedom
        """
        # Load current fit
        try:
            (orig_p_opt, cov, meta) = self.metadata["fitted_result"]
        except KeyError:
            logger.info("Please run a fit first!")
            return None, np.nan
        elem = self.metadata["elements"][0]
        elem_p_opt_name = "log_eps({})".format(elem)
        assert elem_p_opt_name in orig_p_opt, (elem, orig_p_opt)
        abund0 = orig_p_opt[elem_p_opt_name]

        # Set up data again from scratch
        spectrum = self._verify_spectrum(None)
        mask = self.mask(spectrum)
        x = spectrum.dispersion[mask]
        data_y = spectrum.flux[mask]
        data_ivar = spectrum.ivar[mask]
        
        # recompute chi2 for this and for -9
        model_y = self(x, *orig_p_opt.values())
        chi2_current = np.nansum((data_y - model_y)**2 * data_ivar)
        p_opt = orig_p_opt.copy()
        p_opt[elem_p_opt_name] = abund0 - 9
        model_y2 = self(x, *p_opt.values())
        chi2_none = np.nansum((data_y - model_y2)**2 * data_ivar)
        
        delta_chi2 = chi2_none - chi2_current # should be > 0, since chi2_none is worse than the fit
        frac = stats.norm.cdf(sigma) - stats.norm.cdf(-sigma)
        threshold = stats.chi2.ppf(frac, 1) # 1 df for one abundance
        return delta_chi2 > threshold, delta_chi2

    def find_upper_limit(self, sigma=3, start_at_current=True, max_elem_diff=12.0, tol=.01, pix_per_element=1.0):
        """
        Does a simple chi2 check to find the upper limit
        
        (1) fit abundance until converged
        (2) increase abundance until delta chi2 matches the difference from sigma
        
        If start_at_current is True, use the current best-fit abundance as 
        the starting point to find the upper limit.
        Otherwise, run iterfit and start there.
        
        Start with none of the element, then increase abundance
        until it is (x) sigma above.
        """
        assert len(self.metadata["elements"]) == 1, self.metadata["elements"]
        
        if start_at_current:
            abund = self.abundances[0]
            if abund is None or np.isnan(abund):
                logger.debug("find_upper_limit: No fit yet! Running fit...")
                abund = self.iterfit(tol=tol)
        else:
            abund = self.iterfit(tol=tol)
            
        err = self.find_error(sigma=sigma, max_elem_diff=max_elem_diff, abundtol=tol, pix_per_element=pix_per_element)
        upper_limit = abund + err
        
        # Save it as the current abundance and mark as upper limit
        elem = self.metadata["elements"][0]
        key = "log_eps({})".format(elem)
        fitted_result = self.metadata["fitted_result"]
        fitted_result[0][key] = upper_limit
        fitted_result[2]["abundances"][0] = upper_limit
        self.metadata["is_upper_limit"] = True
        
        # Update synthesized figure data
        self.update_fit_after_parameter_change()

        return upper_limit

    def export_line_list(self, fname):
        self.transitions.write(fname, format="moog")
        return None

    def propagate_stellar_parameter_error(self, **kwargs):
        e_Teff, e_logg, e_vt, e_MH = self.session.stellar_parameters_err
        Teff, logg, vt, MH = self.session.stellar_parameters
        alpha = self.session.metadata["stellar_parameters"]["alpha"]
        # Save a copy to reset later
        current_metadata = self.metadata.copy()
        try:
            self.session.set_stellar_parameters(
                Teff, logg, vt, MH, alpha)
            abund0 = self.iterfit(**kwargs)
            self.session.set_stellar_parameters(
                Teff+e_Teff, logg, vt, MH, alpha)
            abund1 = self.iterfit(**kwargs)
            self.session.set_stellar_parameters(
                Teff, logg+e_logg, vt, MH, alpha)
            abund2 = self.iterfit(**kwargs)
            self.session.set_stellar_parameters(
                Teff, logg, vt+e_vt, MH, alpha)
            abund3 = self.iterfit(**kwargs)
            self.session.set_stellar_parameters(
                Teff, logg, vt, MH+e_MH, alpha)
            abund4 = self.iterfit(**kwargs)
            dTeff_error = abund1-abund0
            dlogg_error = abund2-abund0
            dvt_error = abund3-abund0
            dMH_error = abund4-abund0
        except:
            self.metadata = current_metadata
            self.session.set_stellar_parameters(
                Teff, logg, vt, MH, alpha)
            raise
        else:
            self.metadata = current_metadata
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
