#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" An object to manage SMH sessions. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Session"]

import logging
import numpy as np
import os
import sys
import tarfile
import yaml
import time
from six import string_types
from six.moves import cPickle as pickle
from shutil import copyfile, rmtree
from tempfile import mkdtemp

from .linelists import LineList
from . import (photospheres, radiative_transfer, specutils, isoutils, utils)
from .spectral_models import ProfileFittingModel, SpectralSynthesisModel

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

        # We only store this so that we can retain the original file format for
        # when we save the session.
        self._input_spectra_paths = spectrum_paths
        
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


    def save(self, session_path, overwrite=False, **kwargs):
        """
        Save the Session to disk.

        :param session_path:
            The disk location where to save the session.

        :param overwrite: [optional]
            Overwrite the file if it already exists.

        :returns:
            True if the file was saved successfully.

        :raise IOError:
            If the `session_path` already exists and `overwrite` was set to
            `False`.
        """

        if os.path.exists(session_path) and not overwrite:
            raise IOError("path '{}' already exists".format(session_path))

        # What files do we need:
        #   --> input spectra from class.
        #   --> template spectrum from metadata ("rv", "template_spectrum_path")
        #   --> line list from metadata ("line_list")

        # What internals do we need:
        # - the metadata
        #   --> spectral_models

        # Create a temporary directory and copy the required files there.
        twd_paths = []
        metadata = self.metadata.copy()
        protocol = kwargs.pop("protocol", 2)

        # Store the input spectrum paths before we copy them.
        metadata["reconstruct_paths"] = {
            "input_spectra": [os.path.basename(path) \
                for path in self._input_spectra_paths],
        }

        # Create a temporary working directory and copy files over.
        twd = mkdtemp(**kwargs)
        for path in self._input_spectra_paths:
            copyfile(path, os.path.join(twd, os.path.basename(path)))
            twd_paths.append(path)

        # Save the template spectrum.
        if "template_spectrum_path" in metadata.get("rv", {}):

            # TODO: Give a random (unused) path name to the template spectrum?
            path = metadata["rv"].pop("template_spectrum_path")
            new_path = os.path.join(twd, ".{}".format(os.path.basename(path)))
            copyfile(path, new_path)
            twd_paths.append(new_path)
            metadata["reconstruct_paths"]["template_spectrum_path"] \
                = ".{}".format(os.path.basename(path))

        # The line list and session file will always have the same name, since
        # the line list could be constructed from many paths.

        # Save the line list.
        if "line_list" in metadata:
            # TODO: Give a random (unused) path name to the line list?
            metadata.pop("line_list").write(os.path.join(twd, "line_list.fits"),
                format="fits")
            twd_paths.append(os.path.join(twd, "line_list.fits"))

        # The spectral models must be treated with care.
        metadata["spectral_models"] \
            = [_.__getstate__() for _ in metadata.get("spectral_models", [])]

        # Pickle the metadata.
        twd_paths.append(os.path.join(twd, "session.pkl"))
        try:
            with open(twd_paths[-1], "wb") as fp:
                # I dreamt Python 2 was dead. It was great.
                pickle.dump(metadata, fp, protocol)

        except Exception:
            logger.exception("Exception in serializing session:")

            exc_info = sys.exc_info()

            # Go through the key/value pairs and see what can/cannot be
            # saved.
            failed_on = []
            for key, value in metadata.iteritems():
                try:
                    with open(os.path.join("twd", ".damaged", "wb")) as dfp:
                        pickle.dump([key, value], dfp, protocol)
                except:
                    failed_on.append(key)

            logger.warning("Could not pickle the following keys (and their "
                "value pairs): {}".format(", ".join(failed_on)))

            # Now re-raise the original exception.
            raise (exc_info[1], None, exc_info[2])

        # Tar it up.
        if not session_path.lower().endswith(".smh"):
            session_path = "{}.smh".format(session_path)
        tarball = tarfile.open(name=session_path, mode="w:gz")
        for path in twd_paths:
            tarball.add(path, arcname=os.path.basename(path))
        tarball.close()

        # Remove the temporary working directory.
        rmtree(twd)

        return True


    @classmethod
    def load(cls, session_path, **kwargs):
        """
        Create a Session from a path saved to disk.

        :param session_path:
            The disk location where to load the session from.
        """

        # Extract all.
        tarball = tarfile.open(name=session_path, mode="r:gz")
        twd = mkdtemp(**kwargs)
        tarball.extractall(path=twd)

        # Reconstruct the session, starting with the initial paths.
        with open(os.path.join(twd, "session.pkl"), "rb") as fp:
            metadata = pickle.load(fp)

        # Load in the line list.
        if os.path.exists(os.path.join(twd, "line_list.fits")):
            metadata["line_list"] = LineList.read(
                os.path.join(twd, "line_list.fits"), format="fits")


        # Load in the template spectrum.
        # TODO: This means we actually need to keep the template spectrum on
        #       disk until it's loaded. Let's store it now but don't load it in.
        #if "template_spectrum_path" in metadata["reconstruct_paths"]:
        #    metadata["rv"]["template_spectrum_path"] = 

        # Create the object using the temporary working directory input spectra.
        session = cls([os.path.join(twd, basename) \
            for basename in metadata["reconstruct_paths"]["input_spectra"]])

        # Remove any reconstruction paths.
        metadata.pop("reconstruct_paths")

        # Update the new session with the metadata.
        session.metadata = metadata

        # Reconstruct any spectral models.
        reconstructed_spectral_models = []
        for state in session.metadata.get("spectral_models", []):

            args = [session, state["transition_hashes"]]
            if state["type"] == "SpectralSynthesisModel":
                klass = SpectralSynthesisModel
                args.append(state["metadata"]["elements"])

            elif state["type"] == "ProfileFittingModel":
                klass = ProfileFittingModel

            else:
                raise ValueError("unrecognized spectral model class '{}'"\
                    .format(state["type"]))

            model = klass(*args)
            model.metadata = state["metadata"]
            reconstructed_spectral_models.append(model)

        # Update the session with the spectral models.
        session.metadata["spectral_models"] = reconstructed_spectral_models

        # TODO: We need to clean up!
        #       The line list and input spectra are stored in a TWD, and we at
        #       least need to keep the input spectra until the model is saved
        #       later on so that the input_spectra can be copied into the new
        #       'save as' temporary working directory.

        return session


    def index_spectral_models(self):
        """
        (Re-)Index the spectral models so that they are linked correctly
        against the session.
        """

        for spectral_model in self.metadata.get("spectral_models", []):
            spectral_model.index_transitions()
        return None



    def apply_quality_constraints(self, constraints, only=None):
        """
        Apply quality constraints to some or all of the spectral models in this 
        session. Spectral models that do not match the specified constraints are
        marked as unacceptable.

        :param constraints:
            A dictionary specifying the constraints to apply.

        :param only: [optional]
            A function that takes a single argument (a spectral model) and
            returns True or False whether that spectral model should be subject
            to the `constraints`.

        :returns:
            The number of spectral models affected by the quality constraints.
        """

        N = 0
        only = only or (lambda model: True)
        for spectral_model in self.metadata.get("spectral_models", []):
            if not only(spectral_model): continue

            for key, (lower, upper) in constraints.items():

                # Get data.

                if (lower is not None and value < lower) \
                or (upper is not None and value > upper):
                    spectral_model.is_acceptable = False
                    N += 1

        return N



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

        value = self.metadata
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


    def update_default_setting(self, key_tree, value):
        """
        Update a default value in the local settings file.

        :param key_tree:
            A tuple containing a tree of dictionary keys.

        :param value:
            The value for the setting.
        """

        # Open the defaults.
        with open(self._default_settings_path, "rb") as fp:
            defaults = yaml.load(fp)

        branch = defaults
        for key in key_tree[:-1]:
            branch = branch[key]
        branch[key_tree[-1]] = value

        with open(self._default_settings_path, "w") as fp:
            fp.write(yaml.dump(defaults))

        return True





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
            self.metadata["rv"]["template_spectrum_path"] = template_spectrum
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

        filtering = kwargs.pop("filtering",
            lambda model: model.use_for_stellar_parameter_inference)
        for i, model in enumerate(self.metadata["spectral_models"]):
            if filtering(model):

                # TODO assert it is a profile model.
                spectral_model_indices.append(i)
                transition_indices.append(model._transition_indices[0])
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
        return (transitions, slopes)



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

    def measure_abundances(self, spectral_models=None, 
                           save_abundances=True,
                           calculate_uncertainties=True):
        """
        Sort through list of spectral models (default is self.metadata["spectral_models"]).
        Measure all ProfileFittingModels with self.rt.abundance_cog() at the end
        Calculate abundance uncertainties too.
        Measure synthesis on the fly (TODO).
        save_abundances: if True, save measured values into the session
        """
        start = time.time()
        if spectral_models is None:
            spectral_models = self.metadata["spectral_models"]

        equivalent_widths = []
        equivalent_width_errs = []
        transition_indices = []
        spectral_model_indices = []
        
        num_profile = 0
        num_synth = 0
        for i,spectral_model in enumerate(spectral_models):
            if isinstance(spectral_model, ProfileFittingModel):
                spectral_model_indices.append(i)
                transition_indices.extend(spectral_model._transition_indices)
                if spectral_model.is_acceptable:
                    equivalent_widths.append(1000.* \
                        spectral_model.metadata["fitted_result"][-1]["equivalent_width"][0])
                    equivalent_width_errs.append(1000.* \
                        np.nanmax(spectral_model.metadata["fitted_result"][-1]["equivalent_width"][1:3]))
                else:
                    equivalent_widths.append(np.nan)
                    equivalent_width_errs.append(np.nan)
                num_profile += 1
            elif isinstance(spectral_model, SpectralSynthesisModel):
                print("Ignoring synthesis",spectral_model)
                num_synth += 1
            else:
                raise RuntimeError("Unknown model type: {}".format(type(spectral_model)))
            
        if num_profile > 0 and \
        (len(equivalent_widths) == 0 \
        or np.isfinite(equivalent_widths).sum() == 0):
            raise ValueError("no measured transitions to calculate abundances")
        
        transition_indices = np.array(transition_indices)
        spectral_model_indices = np.array(spectral_model_indices)
        transitions = self.metadata["line_list"][transition_indices].copy()
        transitions["equivalent_width"] = equivalent_widths
        finite = np.isfinite(transitions["equivalent_width"])
        
        abundances = self.rt.abundance_cog(
            self.stellar_photosphere, transitions[finite])

        if calculate_uncertainties:
            # Increase EW by uncertainty and measure again
            equivalent_width_errs = np.array(equivalent_width_errs)
            transitions["equivalent_width"] += equivalent_width_errs
            finite_uncertainty = np.isfinite(transitions["equivalent_width"])
            # Some EW uncertainties are HUGE. 
            # Set a maximum EW of 9999, and later max abund uncertainty of 9
            transitions["equivalent_width"][transitions["equivalent_width"] > 9999] = 9999.

            uncertainties = self.rt.abundance_cog(
                self.stellar_photosphere, transitions[finite_uncertainty])
            
            # These are not the same size. Make them the same size by filling with nan
            # Inelegant but works...
            _all = np.zeros(len(finite))*np.nan
            _all[finite_uncertainty] = uncertainties
            uncertainties = _all[finite] - abundances
            # Set a maximum abund uncertainty of 9
            uncertainties[uncertainties > 9] = 9
        else:
            uncertainties = np.nan*np.ones_like(abundances)
        assert len(uncertainties) == len(abundances)
        
        # This appears to save into the session just fine
        if save_abundances:
            for index, abundance, uncertainty in \
                    zip(spectral_model_indices[finite], abundances, uncertainties):
                spectral_models[index].metadata["fitted_result"][-1]["abundances"] = [abundance]
                spectral_models[index].metadata["fitted_result"][-1]["abundance_uncertainties"] = [uncertainty]
        print("Time to measure abundances: {:.1f}".format(time.time()-start))
        return abundances, uncertainties if calculate_uncertainties else abundances

