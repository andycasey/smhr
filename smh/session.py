#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" An object to manage SMH sessions. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Session"]

import logging
import numpy as np
import os
import yaml
from six import string_types

from .linelists import LineList
from . import (photospheres, radiative_transfer, specutils, isoutils, utils)

logger = logging.getLogger(__name__)


class BaseSession(object):
    """
    An abstract class for a SMH session.
    """
    pass


class Session(BaseSession):
    """
    An object to manage a session in SMH.
    """

    # The default settings path is only defined (hard-coded) here.
    _default_settings_path = os.path.expanduser("~/.smh_session.defaults")

    def __init__(self, spectrum_paths, **kwargs):
        """
        Create a new session from input spectra.

        :param spectrum_paths:
            Filename of a single spectrum, or an iterable of spectrum filenames.
        """

        if isinstance(spectrum_paths, string_types):
            spectrum_paths = (spectrum_paths, )

        # Load the spectra and flatten all orders into a single list.
        input_spectra = []
        for path in spectrum_paths:
            s = specutils.Spectrum1D.read(path)
            if isinstance(s, list):
                input_spectra.extend(s)
            else:
                input_spectra.append(s)
        
        # Sort orders from blue to red.
        input_spectra.sort(key=lambda order: order.dispersion.mean())

        # TODO: Store the path names internally for provenance?

        # Extract basic metadata information from the spectrum headers if
        # possible: RA, DEC, OBJECT
        # TODO: Include UTDATE, etc to calculate helio/bary-centric corrections.
        self.metadata = { "NOTES": ""}
        common_metadata_keys = ["RA", "DEC", "OBJECT"] \
            + kwargs.pop("common_metadata_keys", [])

        for key in common_metadata_keys:
            for order in input_spectra:
                if key in order.metadata:
                    self.metadata[key] = order.metadata[key]
                    break

        self.input_spectra = input_spectra
        self.input_spectra_paths = spectrum_paths
        
        # Initialize metadata dictionary.
        N = len(self.input_spectra)
        self.metadata.update({
            "rv": {},
            "normalization": {
                "continuum": [None] * N,
                "normalization_kwargs": [{}] * N
            },
            "isotopes": {},
            "stellar_parameters": {
                "effective_temperature": 5777, # K
                "surface_gravity": 4.4,
                "metallicity": 0.0, # Solar-scaled
                "microturbulence": 1.06, # km/s
            }
        })

        return None


    @property
    def rt(self):
        """
        Access radiative transfer functions.
        """
        
        # What RT do we prefer?
        # TODO: Respect user preferences about which rt they want.
        return radiative_transfer.moog


    def setting(self, key_tree):
        """
        Return a setting value either from the session metadata, or the default.

        :param key_tree:
            A tuple containing a tree of dictionary keys.
        """

        value = self.metadata.copy()
        try:
            for key in key_tree:
                value = value[key]

        except KeyError:
            # Check in defaults.
            with open(self._default_settings_path, "rb") as fp:
                default = yaml.load(fp)

            try:
                for key in key_tree:
                    default = default[key]
            except KeyError:
                return None

            else:
                return default

        else:
            return value


    """
    def _default(self, input_value, default_key_tree):
        ""
        Return the input value if it is valid (i.e., not `None`), or return the
        default session value.

        :param input_value:
            The value provided by the user.

        :param default_key_tree:
            A tuple containing a tree of dictionary keys.
        ""

        if input_value is not None:
            return input_value

        with open(self._default_settings_path, "rb") as fp:
            default = yaml.load(fp)

        for key in default_key_tree:
            try:
                default = default[key]
            except KeyError:
                raise KeyError("no default session value found for {0}".format(
                    default_key_tree))
                
        return default
    """


    @classmethod
    def from_filename(cls, session_path, **kwargs):
        """
        Create a Session from a path saved to disk.
        """

        raise NotImplementedError


    def _get_overlap_order(self, wavelength_regions, template_spectrum=None):
        """
        Find the order (and order index) that most overlaps with the template
        spectrum in any of the wavelength_regions provided.

        :param wavelength_regions:
            A list of the wavelength regions to search for overlap between the 
            template spectrum and the input spectra.

        :param template_spectrum: [optional]
            An optional template spectrum that should also overlap with the
            requested ranges. The should be a `specutils.Spectrum1D` object.

        :returns:
            The overlapping order, the overlap index, and the wavelength_region
            that matched.
        """

        # Check to see if wavelength region is a list of entries.
        try:
            int(wavelength_regions[0])
        except (TypeError, ValueError):
            # It is (probably) a list of 2-length tuples.
            None
        else:
            wavelength_regions = [wavelength_regions]

        # Find the order best suitable for the preferred wavelength region.
        for wl_start, wl_end in wavelength_regions:
            # Does the template cover this range?
            if template_spectrum is not None and \
                not (wl_start > template_spectrum.dispersion[0] \
                and  wl_end   < template_spectrum.dispersion[-1]):
                continue

            # Do any observed orders cover any part of this range?
            overlaps, indices = specutils.find_overlaps(
                self.input_spectra, (wl_start, wl_end), return_indices=True)
            if not overlaps:
                continue

            # The first spectral index has the most overlap with the range.
            overlap_index = indices[0]
            overlap_order = overlaps[0]
            break

        else:
            raise ValueError("no wavelength regions are common to the template "
                             "and the observed spectra")

        return (overlap_order, overlap_index, (wl_start, wl_end))


    def rv_measure(self, template_spectrum=None, wavelength_region=None,
        resample=None, apodize=None, normalization_kwargs=None, **kwargs):
        """
        Measure the observed radial velocity by cross-correlating an individual
        echelle order with a normalized rest-frame template spectrum. The most
        suitable order is determined by the `wavelength_region` given.

        :param template_spectrum: [optional]
            The rest-frame template spectrum to cross-correlate against. This
            should be a `specutils.Spectrum1D` object or a `str`-type pointing
            to a spectrum path.

        :param wavelength_region: [optional]
            The preferred wavelength region(s) to use for cross-correlation. The
            most suitable input spectrum will be determined from this supplied
            range.

        :param resample: [optional]
            Re-sample to the 'template' or 'observed' spectrum.

        :param apodize: [optional]
            What fraction of the pixels to apodize (on both ends) before
            performing the cross-correlation.

        :param normalization_kwargs: [optional]
            Keyword arguments that are passed directly to the 
            `Spectrum1D.fit_continuum` function.

        :param kwargs:
            Dummy variable to take all extra keywords

        Note
        ----
        If these parameters are not specified, then defaults are read from the
        session defaults file.
        """

        # Read in everything from defaults as necessary.
        if template_spectrum is None:
            template_spectrum = self.setting(("rv", "template_spectrum"))
        if wavelength_region is None:
            wavelength_region = self.setting(("rv", "wavelength_regions"))
        if resample is None:
            resample = self.setting(("rv", "resample"))
        if apodize is None:
            apodize = self.setting(("rv", "apodize"))
        if normalization_kwargs is None:
            normalization_kwargs = self.setting(("rv", "normalization"))

        # Is the template spectrum actually a filename?
        if isinstance(template_spectrum, string_types):
            self.metadata["rv"]["template_spectrum_name"] = template_spectrum
            template_spectrum = specutils.Spectrum1D.read(template_spectrum,
                debug=True)

        overlap_order, overlap_index, wavelength_region = \
            self._get_overlap_order(wavelength_region, template_spectrum) 

        # Normalize that order using the normalization settings supplied.
        observed_spectrum, continuum, _, __ = overlap_order.fit_continuum(
            full_output=True, **normalization_kwargs)

        # Perform cross-correlation with the template spectrum.
        rv, rv_uncertainty, ccf = specutils.cross_correlate(
            observed_spectrum, template_spectrum, wavelength_region, 
            apodize=apodize, resample=resample)

        # Store the measured information as part of the session.
        # TODO: Should we store these as a NamedTuple instead?

        try:
            v_helio, v_bary = specutils.motions.corrections_from_headers(
                overlap_order.metadata)
        except Exception as e:
            # TODO not raising an exception for testing purposes, but may want to
            #logger.exception(
            #    "Exception in calculating heliocentric and barycentric motions")
            logger.error(
                "Exception in calculating heliocentric and barycentric motions")
            logger.error(e)
            v_helio, v_bary = (np.nan, np.nan)

        else:
            v_helio = v_helio.to("km/s").value
            v_bary = v_bary.to("km/s").value

        self.metadata["rv"].update({
            # Measurements
            "rv_measured": rv,
            "rv_uncertainty": rv_uncertainty,
            "order_index": overlap_index,
            "normalized_order": observed_spectrum,
            "continuum": continuum,
            "ccf": ccf,
            "heliocentric_correction": v_helio,
            "barycentric_correction": v_bary,

            # Input settings
            "template_spectrum": template_spectrum,
            "wavelength_region": wavelength_region,
            "resample": resample,
            "apodize": apodize,
            "normalization": normalization_kwargs.copy()
        })

        return (rv, rv_uncertainty)


    def rv_correct(self, rv):
        """
        Apply a radial velocity correction to the input spectra.
        
        :param rv:
            The radial velocity correction (in km/s) to apply.
        """

        self.metadata["rv"]["rv_applied"] = -float(rv)
        return None


    def stitch_and_stack(self, **kwargs):

        normalized_orders = []
        for i, (spectrum, continuum) \
        in enumerate(zip(self.input_spectra,
        self.metadata["normalization"]["continuum"])):
            normalized_orders.append(specutils.Spectrum1D(spectrum.dispersion,
                spectrum.flux / continuum,
                continuum * spectrum.ivar * continuum))

        self.normalized_spectrum = specutils.spectrum.stitch(normalized_orders)

        # Ensure the radial velocity is accounted for.
        self.normalized_spectrum.redshift(v=self.metadata["rv"]["rv_applied"])

        return self.normalized_spectrum


    def normalize_input_spectra(self, **kwargs):
        """
        Continuum-normalize all orders in the input spectra.
        """

        # Get the relevant masks.
        mask_name = self.setting(("normalization", "default_mask"))
        try:
            mask = self.setting(("normalization", "masks"))[mask_name]
        except KeyError:
            mask = {}

        # Get the default input parameters, and overwrite them with any inputs
        # provided to this function.



        # Fit & store continuum for all input spectra.
        raise NotImplementedError


    @property
    def stellar_photosphere(self):
        """
        Return a photosphere model of the current stellar parameters.
        """

        # TODO: HACK -- how to specify different models?
        #               or optional arguments, e.g. [alpha/fe] for CK 2004

        try:
            self._photosphere_interpolator
        except AttributeError:
            self._photosphere_interpolator = photospheres.interpolator()

        meta = self.metadata["stellar_parameters"]
        photosphere = self._photosphere_interpolator(*[meta[k] for k in \
            ("effective_temperature", "surface_gravity", "metallicity")])
        # Update other metadata (e.g., microturbulence)
        # TODO: Convert session.metadata to session.meta to be consistent
        photosphere.meta["stellar_parameters"].update(meta)

        return photosphere


    def stellar_parameter_state(self, full_output=False, **kwargs):
        """
        Calculate the abundances of all spectral models that are used in the
        determination of stellar parameters.
        """

        # Get the transitions & EWs together from spectral models.
        equivalent_widths = []
        transition_indices = []
        spectral_model_indices = []
        for i, model in enumerate(self.metadata["spectral_models"]):
            if model.use_for_stellar_parameter_inference:

                # TODO assert it is a profile model.
                spectral_model_indices.append(i)
                transition_indices.extend(model._transition_indices)
                if model.is_acceptable:
                    equivalent_widths.append(1e3 * \
                        model.metadata["fitted_result"][-1]["equivalent_width"][0])
                else:
                    equivalent_widths.append(np.nan)


        if len(equivalent_widths) == 0 \
        or np.isfinite(equivalent_widths).sum() == 0:
            raise ValueError("no measured transitions to calculate abundances")

        # Construct a copy of the line list table.
        transition_indices = np.array(transition_indices)
        spectral_model_indices = np.array(spectral_model_indices)
        transitions = self.metadata["line_list"][transition_indices].copy()
        transitions["equivalent_width"] = equivalent_widths

        finite = np.isfinite(transitions["equivalent_width"])

        # Calculate abundances and put them back into the spectral models stored
        # in the session metadata.
        abundances = self.rt.abundance_cog(
            self.stellar_photosphere, transitions[finite])


        for index, abundance in zip(spectral_model_indices[finite], abundances):
            self.metadata["spectral_models"][index]\
                .metadata["fitted_result"][-1]["abundances"] = [abundance]

        transitions["abundance"] = np.nan * np.ones(len(transitions))
        transitions["abundance"][finite] = abundances

        # By default just return a transitions table for convenience.
        if not full_output:
            return transitions

        transitions["reduced_equivalent_width"] = np.log10(1e-3 * \
            transitions["equivalent_width"] / transitions["wavelength"])

        slopes = None
        #slopes = utils.equilibrium_state(transitions,
        #    ("expot", "reduced_equivalent_width", "wavelength"))

        # Otherwise return full state information.
        return (transitions, slopes, spectral_model_indices)



    def optimize_stellar_parameters(self, **kwargs):
        """
        Optimize the stellar parameters for this star using the spectral models
        associated with stellar parameter inference.
        """

        # Get the list of relevant spectral models.

        # Any synth?
        # raise NotImplementedError yet.

        # interpolator, do obj. function

        raise NotImplementedError
