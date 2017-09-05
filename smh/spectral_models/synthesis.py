
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

from .base import BaseSpectralModel
from smh import utils
from smh.specutils import Spectrum1D
from smh.photospheres.abundances import asplund_2009 as solar_composition
# HACK TODO REVISIT SOLAR SCALE: Read from session defaults?


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
        logger.debug(abundances)

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
        print(*parameters)
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
            "window": 1, 
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
        logger.debug(unique_elements)
        logger.debug(type(unique_elements))
        unique_elements = np.unique(unique_elements)
        
        rt_abundances = {}
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

                try:
                    fitted_result = self.metadata["fitted_result"]
                except KeyError:
                    # Assume scaled-solar composition.
                    p0.append(solar_composition(element) + \
                        self.session.metadata["stellar_parameters"]["metallicity"])
                else:
                    p0.append(fitted_result[0][parameter])
                
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


    def fit(self, spectrum=None, **kwargs):
        """
        Fit a synthesised model spectrum to the observed spectrum.

        :param spectrum: [optional]
            The observed spectrum to fit the synthesis spectral model. If None
            is given, this will default to the normalized rest-frame spectrum in
            the parent session.
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
                
            p_opt, p_cov = op.curve_fit(objective_function, xdata=x, ydata=y,
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

        fitting_metadata = {
            "model_x": x,
            "model_y": model_y,
            "model_yerr": model_yerr,
            "residual": residuals,
            "chi_sq": chi_sq,
            "dof": dof,
            "abundances": p_opt[:len(self.elements)]
        }

        # Convert result to ordered dict.
        named_p_opt = OrderedDict(zip(self.parameter_names, p_opt))
        self.metadata["fitted_result"] = (named_p_opt, cov, fitting_metadata)
        self.metadata["is_acceptable"] = True
        if "vrad" in self.parameter_names:
            self.metadata["manual_rv"] = named_p_opt["vrad"]
        if "sigma_smooth" in self.parameter_names:
            self.metadata["manual_sigma_smooth"] = named_p_opt["sigma_smooth"]
        
        return self.metadata["fitted_result"]


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
        synth_spec = Spectrum1D(meta["model_x"], meta["model_y"], ivar)
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
    
    def update_fit_after_parameter_change(self, synthesize=False):
        """
        After manual parameter change, update fitting metadata 
        for plotting purposes.
        This somewhat messes up the state, but if you are trying to
        fiddle manually you probably don't care about those anyway.
        
        :param synthesize:
            Not implemented. Eventually should not recalc synth every time.
            But have to figure out when you need to recalc synth,
            and where to store the synth.
        """
        self._update_parameter_names()
        spectrum = self._verify_spectrum(None)
        try:
            (named_p_opt, cov, meta) = self.metadata["fitted_result"]
        except KeyError:
            logger.info("Please run a fit first!")
            return None
        model_x = meta["model_x"]
        model_y = self(model_x, *named_p_opt.values())
        #model_yerr = np.nan * np.ones((2, model_x.size))
        model_yerr = np.nan * np.ones((1, model_x.size))
        #model_yerr = np.max(np.abs(model_yerr), axis=0)
        residuals = (meta["residual"] + meta["model_y"]) - model_y
        model_x, model_y, model_yerr, residuals = self._fill_masked_arrays(
            spectrum, model_x, model_y, model_yerr, residuals)
        meta["model_y"] = model_y
        meta["model_yerr"] = model_yerr
        meta["residual"] = residuals
        return None

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


    def find_upper_limit(self, elem, sigma=3):
        """
        Does a simple chi2 check to find the upper limit 
        """
        assert len(self.metadata["elements"]) == 1, self.metadata["elements"]
        assert self.metadata["elements"][0] == elem, (self.metadata["elements"], elem)
        
        

    """
    def abundances(self):
        "
        Calculate the abundances (model parameters) given the current stellar
        parameters in the parent session.
        "

        #_ self.fit(self.session.normalized_spectrum)

        # Parse the abundances into the right format.
        raise NotImplementedError("nope")
    """
        
