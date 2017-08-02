#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" An object to manage SMH sessions. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Session"]

import atexit
import logging
import numpy as np
import os
import sys
import tarfile
import yaml
import time
from six import string_types, iteritems
from six.moves import cPickle as pickle
from shutil import copyfile, rmtree
#from tempfile import mkdtemp

import astropy.table
from .linelists import LineList
from .utils import mkdtemp
from . import (photospheres, radiative_transfer, specutils, isoutils, utils)
from .spectral_models import ProfileFittingModel, SpectralSynthesisModel
from smh.photospheres.abundances import asplund_2009 as solar_composition
from . import (smh_plotting)

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

    def __init__(self, spectrum_paths, twd=None, **kwargs):
        """
        Create a new session from input spectra.

        :param spectrum_paths:
            Filename of a single spectrum, or an iterable of spectrum filenames.
        """

        if isinstance(spectrum_paths, string_types):
            spectrum_paths = (spectrum_paths, )

        if twd is None:
            twd = mkdtemp()
        self.twd = twd
        logger.info("Working directory: {}".format(twd))

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
                "alpha":0.4,
            }
        })

        # Load any line list?
        line_list_filename = self.setting(("line_list_filename",))
        if line_list_filename is not None and os.path.exists(line_list_filename):
            self.metadata["line_list"] = LineList.read(
                line_list_filename, format="fits")

        # Construct any default spectral models.
        deconstructed_spectral_models = self.setting(("default_spectral_models", ))
        if deconstructed_spectral_models is not None:

            # Reconstruct them.
            self.metadata.setdefault("spectral_models", [])

            for state in deconstructed_spectral_models:

                args = [self, state["transition_hashes"]]
                if state["type"] == "SpectralSynthesisModel":
                    klass = SpectralSynthesisModel
                    args.append(state["metadata"]["elements"])

                elif state["type"] == "ProfileFittingModel":
                    klass = ProfileFittingModel

                else:
                    raise ValueError("unrecognized spectral model class '{}'"\
                        .format(state["type"]))

                model = klass(*args)
                model.metadata.update(state["metadata"])
                self.metadata["spectral_models"].append(model)

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

        if not session_path.lower().endswith(".smh"):
            session_path = "{}.smh".format(session_path)

        if os.path.exists(session_path) and not overwrite:
            raise IOError("path '{}' already exists".format(session_path))

        metadata = self.metadata.copy()
        protocol = kwargs.pop("protocol", 2)

        # Create a temporary working directory and copy files over.
        twd = mkdtemp(**kwargs)
        twd_paths = [] + self._input_spectra_paths

        # Input spectra.
        metadata["reconstruct_paths"] = {
            "input_spectra": list(map(os.path.basename, twd_paths))
        }


        def safe_path(path, twd, metadata):
            """
            Function to ensure that we don't have filename clashes in the
            temporary working directory.
            """

            basename = os.path.basename(path)
            existing_paths = np.hstack(metadata["reconstruct_paths"].values())

            if basename in existing_paths:

                while basename in existing_paths:
                    basename = ".".join([
                        utils.random_string(), os.path.basename(path)])

                new_path = os.path.join(twd, basename)

                if os.path.exists(path):
                    copyfile(path, new_path)

                return new_path

            return path


        # Radial velocity template.
        if "template_spectrum_path" in metadata.get("rv", {}) and \
        os.path.exists(metadata["rv"]["template_spectrum_path"]):
            twd_paths.append(safe_path(
                metadata["rv"].pop("template_spectrum_path"), twd, metadata))
            metadata["reconstruct_paths"]["template_spectrum_path"] \
                = os.path.basename(twd_paths[-1])

        # Line list.
        if "line_list" in metadata:
            twd_paths.append(safe_path(os.path.join(twd, "line_list.fits"),
                twd, metadata))
            metadata.pop("line_list").write(twd_paths[-1], format="fits")
            metadata["reconstruct_paths"]["line_list"] \
                = os.path.basename(twd_paths[-1])
            

        # normalized spectrum.
        #twd_path = safe_path(os.path.join(twd, "normalized_spectrum.fits"),
        #    twd, metadata)
        ## HACK while fits is broken, it should be called .txt
        twd_path = safe_path(os.path.join(twd, "normalized_spectrum.txt"),
            twd, metadata)
        try:
            self.normalized_spectrum.write(twd_path)

        except AttributeError:
            None

        else:
            twd_paths.append(twd_path)
            metadata["reconstruct_paths"]["normalized_spectrum"] \
                = os.path.basename(twd_path)

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
            for key, value in iteritems(metadata):
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
        exception_occurred = False
        tarball = tarfile.open(name=session_path, mode="w:gz")
        for path in twd_paths:
            try:
                tarball.add(path, arcname=os.path.basename(path))

            except:
                logger.exception(
                    "Cannot save path '{}' to session:".format(path))

                logger.warn(
                    "Continuing to save session without this path "
                    "before raising the issue")
                exception_occurred = True
                continue

        tarball.close()

        # Remove the temporary working directory.
        rmtree(twd)

        if exception_occurred:
            raise

        logger.info("Saved file to {}".format(session_path))

        return True


    def import_transitions(self, path):
        """
        Import transitions (line list data and spectral models) from disk.

        :param path:
            The disk location of the serialized transitions.
        """

        with open(path, "rb") as fp:
            line_list, spectral_model_states = pickle.load(fp)

        # Integrate line list with the existing list.
        if "line_list" in self.metadata:
            self.metadata["line_list"] = self.metadata["line_list"].merge(
                line_list, in_place=False)
        else:
            self.metadata["line_list"] = line_list

        # Add the spectral models.
        self.metadata.setdefault("spectral_models", [])

        reconstructed_spectral_models = []
        for state in spectral_model_states:

            args = [self, state["transition_hashes"]]
            if state["type"] == "SpectralSynthesisModel":
                klass = SpectralSynthesisModel
                args.append(state["metadata"]["elements"])

            elif state["type"] == "ProfileFittingModel":
                klass = ProfileFittingModel

            else:
                raise ValueError("unrecognized spectral model class '{}'"\
                    .format(state["type"]))

            model = klass(*args)
            model.metadata.update(state["metadata"])
            reconstructed_spectral_models.append(model)

        self.metadata["spectral_models"].extend(reconstructed_spectral_models)
        return len(reconstructed_spectral_models)


    # TODO: put export_transitions here from the spectral model GUI too?


    @classmethod
    def load(cls, session_path, skip_spectral_models=False, **kwargs):
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
        line_list = metadata["reconstruct_paths"].get("line_list", None)
        if line_list is not None:
            metadata["line_list"] = LineList.read(
                os.path.join(twd, line_list), format="fits")
            metadata["line_list_argsort_hashes"] = np.argsort(
                metadata["line_list"]["hash"])

        # Load in the template spectrum.
        template_spectrum_path \
            = metadata["reconstruct_paths"].get("template_spectrum_path", None)
        if template_spectrum_path is not None:
            metadata["rv"]["template_spectrum_path"] \
                = os.path.join(twd, template_spectrum_path)

        # Create the object using the temporary working directory input spectra.
        # TODO #225 and other issues have had major saving/loading problems because
        # the temporary directories get deleted.
        session = cls([os.path.join(twd, basename) \
            for basename in metadata["reconstruct_paths"]["input_spectra"]],
            twd=twd)
        
        # Load in any normalized spectrum.
        normalized_spectrum \
            = metadata["reconstruct_paths"].get("normalized_spectrum", None)
        if normalized_spectrum is not None:
            session.normalized_spectrum = specutils.Spectrum1D.read(
                os.path.join(twd, normalized_spectrum))

        # Remove any reconstruction paths.
        metadata.pop("reconstruct_paths")

        # Update the new session with the metadata.
        session.metadata = metadata
        # A hack to maintain backwards compatibility
        session.metadata["stellar_parameters"].setdefault("alpha", 0.4)

        if skip_spectral_models:
            atexit.register(rmtree, twd)
            return session

        # Reconstruct any spectral models.
        reconstructed_spectral_models = []
        start = time.time()
        for state in session.metadata.get("spectral_models", []):
            start2 = time.time()

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
            model.metadata.update(state["metadata"])
            reconstructed_spectral_models.append(model)
            t2 = time.time()-start2
            if t2 > 1.0:
                logger.debug("  Long time to load model {:.3f} {} {} {}\n".format(t2, len(model._transitions), model.elements, model.wavelength))
        logger.debug("Time to reconstruct spectral models: {:.3f}".format(time.time()-start))
        
        # Update the session with the spectral models.
        session.metadata["spectral_models"] = reconstructed_spectral_models

        # Clean up the TWD when Python exits.
        atexit.register(rmtree, twd)

        logger.info("Loaded file {}".format(session_path))
        logger.debug("Input spectra paths: {}".format(session._input_spectra_paths))

        return session


    def index_spectral_models(self):
        """
        (Re-)Index the spectral models so that they are linked correctly
        against the session.
        """

        self.metadata["line_list_argsort_hashes"] = np.argsort(
            self.metadata["line_list"]["hash"])
        for spectral_model in self.metadata.get("spectral_models", []):
            spectral_model.index_transitions()
        return None


    def apply_spectral_model_quality_constraints(self, constraints, only=None,
        full_output=False):
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

        :param full_output: [optional]
            Optionally return the number of models affected and their indices.

        :returns:
            The number of spectral models affected by the quality constraints,
            and optionally, the indices of those affected spectral models.
        """

        indices = []
        only = only or (lambda model: True)
        self.metadata.setdefault("spectral_models", [])
        for i, spectral_model in enumerate(self.metadata["spectral_models"]):
            if not only(spectral_model) \
            or not spectral_model.is_acceptable \
            or spectral_model.is_upper_limit: continue

            if not spectral_model.apply_quality_constraints(constraints):
                indices.append(i)

        N = len(indices)
        return N if not full_output else (N, indices)


    @property
    def rt(self):
        """
        Access radiative transfer functions.
        """
        
        # What RT do we prefer?
        # TODO: Respect user preferences about which rt they want.
        return radiative_transfer.moog


    def setting(self, key_tree, default_return_value=None):
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
                return default_return_value

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
            branch.setdefault(key, {})
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

        logger.info(
            "Heliocentric velocity correction: {0:.2f} km/s".format(v_helio))
        logger.info(
            "Barycentric velocity correction: {0:.2f} km/s".format(v_bary))
        
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
            ("effective_temperature", "surface_gravity", "metallicity", "alpha")])
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
        ews = []
        rews = []
        ew_uncertainties = []

        transition_indices = []
        spectral_model_indices = []

        filtering = kwargs.pop("filtering",
            lambda model: model.use_for_stellar_parameter_inference)
        for i, model in enumerate(self.metadata["spectral_models"]):

            transition_indices.append(model._transition_indices[0])
            
            # TODO assert it is a profile model.
            if filtering(model):
                spectral_model_indices.append(i)

                if model.is_acceptable and not model.is_upper_limit:

                    meta = model.metadata["fitted_result"][-1]
                    model_ew = meta["equivalent_width"]

                    ews.append(1e3 * model_ew[0])
                    rews.append(meta["reduced_equivalent_width"][0])

                    # Get the largest absolute uncertainty.
                    ew_uncertainties.append(1e3 * np.abs(model_ew[1:]).max())

                else:
                    ews.append(np.nan)
                    rews.append(np.nan)
                    ew_uncertainties.append(np.nan)
            else:
                spectral_model_indices.append(np.nan)
                ews.append(np.nan)
                rews.append(np.nan)
                ew_uncertainties.append(np.nan)


        if len(ews) == 0 \
        or np.isfinite(ews).sum() == 0:
            raise ValueError("no measured transitions to calculate abundances")


        # Construct a copy of the line list table.
        transition_indices = np.array(transition_indices)
        spectral_model_indices = np.array(spectral_model_indices)
        transitions = self.metadata["line_list"][transition_indices].copy()
        transitions["equivalent_width"] = ews
        transitions["reduced_equivalent_width"] = rews

        min_eqw = kwargs.pop("minimum_equivalent_width", 0.01) # mA
        finite = np.logical_and(np.isfinite(transitions["equivalent_width"]),
                                transitions["equivalent_width"] > min_eqw)

        # Calculate abundances and put them back into the spectral models stored
        # in the session metadata.
        abundances = self.rt.abundance_cog(
            self.stellar_photosphere, transitions[finite], twd=self.twd)

        for index, abundance in zip(spectral_model_indices[finite], abundances):
            self.metadata["spectral_models"][int(index)]\
                .metadata["fitted_result"][-1]["abundances"] = [abundance]

        transitions["abundance"] = np.nan * np.ones(len(transitions))
        transitions["abundance"][finite] = abundances

        # Calculate abundance uncertainties by propagating the positive
        # uncertainty in equivalent width.

        _transitions = transitions.copy()
        _transitions["equivalent_width"] += ew_uncertainties
        finite = np.isfinite(_transitions["equivalent_width"])

        propagated_abundances = self.rt.abundance_cog(
            self.stellar_photosphere, _transitions[finite], twd=self.twd)

        for index, abundance, propagated_abundance \
        in zip(spectral_model_indices[finite], transitions["abundance"][finite],
            propagated_abundances):

            self.metadata["spectral_models"][int(index)]\
                .metadata["fitted_result"][-1]["abundance_uncertainties"] \
                    = [propagated_abundance - abundance]

        transitions["abundance_uncertainty"] = np.nan * np.ones(len(transitions))
        transitions["abundance_uncertainty"][finite] \
            = propagated_abundances - transitions["abundance"][finite]

        # Include minimum uncertainty of 0.005 dex.
        transitions["abundance_uncertainty"][finite] \
            = np.clip(transitions["abundance_uncertainty"][finite], 0.005, np.inf)

        # By default just return a transitions table for convenience.
        if not full_output:
            return transitions

        # TODO TAKE THIS FROM THE SPECTRAL MODEL METADATA.
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
        min_eqw = .01
        finite = np.logical_and(np.isfinite(transitions["equivalent_width"]),
                                transitions["equivalent_width"] > min_eqw)
        
        abundances = self.rt.abundance_cog(
            self.stellar_photosphere, transitions[finite], twd=self.twd)

        if calculate_uncertainties:
            # Increase EW by uncertainty and measure again
            equivalent_width_errs = np.array(equivalent_width_errs)
            transitions["equivalent_width"] += equivalent_width_errs
            finite_uncertainty = np.logical_and(np.isfinite(transitions["equivalent_width"]),
                                                transitions["equivalent_width"] > min_eqw)
            # Some EW uncertainties are HUGE. 
            # Set a maximum EW of 9999, and later max abund uncertainty of 9
            transitions["equivalent_width"][transitions["equivalent_width"] > 9999] = 9999.

            uncertainties = self.rt.abundance_cog(
                self.stellar_photosphere, transitions[finite_uncertainty], twd=self.twd)
            
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
        print("Time to measure {} abundances: {:.1f}".format(np.sum(finite), time.time()-start))
        return abundances, uncertainties if calculate_uncertainties else abundances

    def summarize_spectral_models(self, spectral_models=None, organize_by_element=False,
                                  use_weights = None, use_finite = True):
        """
        Loop through all spectral_models and return a summary dict

        Returns:
        summary_dict[key] = [num_models, logeps, stdev, stderr, XH, XFe]

        :param spectral_models:
            List of spectral models. Defaults to self.metadata["spectral_models"]
        :param organize_by_element:
            If False (default), key is species (without isotopes)
            If True, key is element (sum all species together)
        :param use_weights:
            TODO Not implemented yet!
            If True, use line-by-line weights
            If False, weight all lines equally
            Defaults to session settings
        :param use_finite:
            If True (default), only use finite abundances
            If False, use any acceptable abundances
            I cannot imagine why you'd set it to False unless debugging
        """
        what_key_type = "element" if organize_by_element else "species"

        start = time.time()
        if spectral_models is None:
            spectral_models = self.metadata.get("spectral_models", [])

        all_logeps = {}
        # TODO abundance uncertainties too
        def update_one_key(key,logeps):
            if key not in all_logeps:
                all_logeps[key] = []
            all_logeps[key].append(logeps)
            return
        for spectral_model in spectral_models:
            if not spectral_model.is_acceptable or spectral_model.is_upper_limit: continue
            
            abundances = spectral_model.abundances
            if abundances is None: continue

            for elem,species,logeps in zip(spectral_model.elements, \
                                          spectral_model.species, abundances):
                # Elements doesn't have the ionization
                # Species do
                if organize_by_element:
                    key = elem
                else:
                    key = species

                if isinstance(key,list):
                    # In synthesis, species may be a list
                    for species in key:
                        if isinstance(species,float):
                            update_one_key(species, logeps)
                        elif isinstance(species,list):
                            for _species in species: update_one_key(_species, logeps)
                        else:
                            raise TypeError("key {} is of type {} (organized by {})".format(\
                                    species,type(species),what_key_type))
                elif isinstance(key,float):
                    assert not organize_by_element
                    update_one_key(key, logeps)
                elif isinstance(key,str):
                    assert organize_by_element
                    update_one_key(key, logeps)
                else:
                    raise TypeError("key {} is of type {} (organized by {})".format(\
                            key,type(key),what_key_type))

        summary_dict = {}
        for key in all_logeps:
            logepss = np.array(all_logeps[key])
            if use_finite:
                finite = np.isfinite(logepss)
                num_models = np.sum(finite)
                logepss = logepss[finite]
            else:
                num_models = len(logepss)

            logeps = np.mean(logepss)
            # TODO weight
            stdev = np.std(logepss)
            stderr= stdev/np.sqrt(num_models)
            XH = logeps - solar_composition(key)
            summary_dict[key] = [num_models, logeps, stdev, stderr, XH, np.nan]

        try:
            if organize_by_element:
                FeH = summary_dict['Fe'][4]
            else:
                # TODO: using Fe I for now, should make this configurable
                FeH = summary_dict[26.0][4]
        except KeyError:
            # Fe not measured yet
            FeH = np.nan
        for key in all_logeps:
            summary = summary_dict[key]
            summary[5] = summary[4] - FeH
            summary_dict[key] = summary

        total_num_models_summarized = np.sum([len(x) for x in all_logeps.values()])
        ## This is basically instantaneous, which is good!
        #print("Time to summarize {} measurements (organized by {}): {:.1f}".format(\
        #        total_num_models_summarized, what_key_type, time.time()-start))
        return summary_dict
    
    def export_abundance_table(self, filepath):
        ## TODO: put in upper limits too.
        summary_dict = self.summarize_spectral_models()
        if filepath.endswith(".tex"):
            self._export_latex_abundance_table(filepath, summary_dict)
        else:
            self._export_ascii_abundance_table(filepath, summary_dict)
        logger.info("Exported to {}".format(filepath))
        return None
    def _export_latex_abundance_table(self, filepath, summary_dict):
        raise NotImplementedError
    def _export_ascii_abundance_table(self, filepath, summary_dict):
        out = np.zeros((len(summary_dict), 7))
        for i,(species, (N, logeps, stdev, stderr, XH, XFe)) in \
                enumerate(iteritems(summary_dict)):
            out[i,:] = [species, N, logeps, stdev, stderr, XH, XFe]
        names = ["species", "N", "logeps", "stdev", "stderr", "[X/H]", "[X/Fe]"]
        tab = astropy.table.Table(out, names=names)
        tab["N"].format = ".0f"
        tab["logeps"].format = "5.2f"
        tab["stdev"].format = "5.2f"
        tab["stderr"].format = "5.2f"
        tab["[X/H]"].format = "5.2f"
        tab["[X/Fe]"].format = "5.2f"
        tab.write(filepath, format="ascii.fixed_width_two_line")
        return True #raise NotImplementedError

    def export_spectral_model_measurements(self, filepath):
        ## TODO: 
        ## Make sure to include synthesis measurements.
        ## We'll eventually put in upper limits too.
        spectral_models = self.metadata.get("spectral_models", [])
        linedata = np.zeros((len(spectral_models), 6)) + np.nan
        for i,spectral_model in enumerate(spectral_models):
            # TODO include upper limits
            if not spectral_model.is_acceptable or spectral_model.is_upper_limit: continue
            # TODO make this work with syntheses as well
            if isinstance(spectral_model, SpectralSynthesisModel):
                raise NotImplementedError
            elif isinstance(spectral_model, ProfileFittingModel):
                line = spectral_model.transitions[0]
                wavelength = line['wavelength']
                species = line['species']
                expot = line['expot']
                loggf = line['loggf']

                try:
                    EW = 1000.*spectral_model.metadata["fitted_result"][2]["equivalent_width"][0]
                    logeps = spectral_model.abundances[0]
                except Exception as e:
                    print(e)
                    EW = np.nan
                    logeps = np.nan
                if EW is None: EW = np.nan
                if logeps is None: logeps = np.nan
            else:
                raise NotImplementedError
            linedata[i,:] = [species, wavelength, expot, loggf, EW, logeps]
        ii_bad = np.logical_or(np.isnan(linedata[:,5]), np.isnan(linedata[:,4]))
        linedata = linedata[~ii_bad,:]

        if filepath.endswith(".tex"):
            self._export_latex_measurement_table(filepath, linedata)
        else:
            self._export_ascii_measurement_table(filepath, linedata)
        logger.info("Exported to {}".format(filepath))
        return None
    def _export_latex_measurement_table(self, filepath, linedata):
        raise NotImplementedError
    def _export_ascii_measurement_table(self, filepath, linedata):
        names = ["species", "wavelength", "expot", "loggf", "EW", "logeps"]
        tab = astropy.table.Table(linedata, names=names)
        tab.sort(["species","wavelength","expot"])
        tab["wavelength"].format = ".3f"
        tab["expot"].format = "5.3f"
        tab["loggf"].format = "6.3f"
        tab["EW"].format = "6.2f"
        tab["logeps"].format = "6.3f"
        tab.write(filepath, format="ascii.fixed_width_two_line")
        return True

    def make_summary_plot(self, figure=None):
        with open(self._default_settings_path, "rb") as fp:
            defaults = yaml.load(fp)
        if "summary_figure" not in defaults:
            raise RuntimeError("Defaults file ({}) must have summary_figure".format(\
                    self._default_settings_path))
        if not isinstance(self.normalized_spectrum, specutils.Spectrum1D):
            print("Must have normalized spectrum to make summary plot")
            return None
        return smh_plotting.make_summary_plot(defaults["summary_figure"],
                                       self.normalized_spectrum, figure)
    def make_ncap_summary_plot(self, figure=None):
        with open(self._default_settings_path, "rb") as fp:
            defaults = yaml.load(fp)
        if "summary_figure_ncap" not in defaults:
            raise RuntimeError("Defaults file ({}) must have summary_figure".format(\
                    self._default_settings_path))
        if not isinstance(self.normalized_spectrum, specutils.Spectrum1D):
            print("Must have normalized spectrum to make summary plot")
            return None
        return smh_plotting.make_summary_plot(defaults["summary_figure_ncap"],
                                       self.normalized_spectrum, figure)

    def import_linelists(self, filenames, ignore_conflicts=False, full_output=False):
        if isinstance(filenames, string_types):
            filenames = [filenames]
        
        line_list = LineList.read(filenames[0], verbose=True)

        filename_transitions = [line_list]
        for filename in filenames[1:]:
            new_lines = LineList.read(filename)
            # Use extremely intolerant to force hashes to be the same
            line_list = line_list.merge(new_lines, in_place=False,
                                        skip_exactly_equal_lines=True,
                                        ignore_conflicts=ignore_conflicts)
            filename_transitions.append(new_lines)

        # Merge the line list with any existing line list in the session.
        if self.metadata.get("line_list", None) is None:
            self.metadata["line_list"] = line_list
            N = len(line_list)
        else:
            N = len(self.metadata["line_list"]) - len(line_list)
            self.metadata["line_list"] \
                = self.metadata["line_list"].merge(
                    line_list, in_place=False, skip_exactly_equal_lines=True,
                    ignore_conflicts=ignore_conflicts)

        # Must update hash sorting after any modification to line list
        self.metadata["line_list_argsort_hashes"] = np.argsort(
            self.metadata["line_list"]["hash"])
        
        if full_output:
            return line_list, filename_transitions
        return line_list
    
    def import_transitions_with_measured_equivalent_widths(self, filenames=None, ignore_conflicts=False):
        line_list = self.import_linelists(filenames, ignore_conflicts=ignore_conflicts)
        try:
            line_list["equivalent_width"]
        except KeyError:
            raise KeyError("no equivalent widths found in imported line lists")
        
        spectral_models_to_add = []
        for idx in range(len(line_list)):
            model = ProfileFittingModel(self, line_list["hash"][[idx]])
            model.metadata.update({
                "is_acceptable": True,
                "fitted_result": [None, None, {
                    # We assume supplied equivalent widths are in milliAngstroms
                    "equivalent_width": \
                    (1e-3 * line_list["equivalent_width"][idx], 0.0, 0.0),
                    "reduced_equivalent_width": \
                    (-3+np.log10(line_list["equivalent_width"][idx]/line_list["wavelength"][idx]),
                      0.0, 0.0)
                }]
            })
            spectral_models_to_add.append(model)
        self.metadata.setdefault("spectral_models", [])
        self.metadata["spectral_models"].extend(spectral_models_to_add)
        self._spectral_model_conflicts = utils.spectral_model_conflicts(
            self.metadata["spectral_models"],
            self.metadata["line_list"])
        
