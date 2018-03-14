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
from astropy.io import ascii
from .linelists import LineList
from .utils import mkdtemp
from . import (photospheres, radiative_transfer, specutils, isoutils, utils)
from .spectral_models import ProfileFittingModel, SpectralSynthesisModel
from smh.photospheres.abundances import asplund_2009 as solar_composition
from . import (smh_plotting, __version__)

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

    def __init__(self, spectrum_paths, twd=None, from_load=False, **kwargs):
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
        self.metadata = {
            "VERSION": __version__,
            "NOTES": "Spectrum paths: "+",".join(spectrum_paths)+"\n\n"
        }
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
                "effective_temperature":
                    self.setting("default_Teff",5777), # K
                "surface_gravity":
                    self.setting("default_logg",4.4),
                "metallicity":
                    self.setting("default_MH",0.0), # Solar-scaled
                "microturbulence":
                    self.setting("default_vt",1.06), # km/s
                "alpha":
                    self.setting("default_aFe",0.4),
            }
        })

        # Set defaults for metadata dictionary
        self.metadata.setdefault("spectral_models", [])
        self.metadata.setdefault("reconstruct_copied_paths", [])

        # Only do these things if creating session for the first time
        if not from_load:
            # Construct default profile models
            line_list_filename = self.setting(("line_list_filename",))
            if line_list_filename is not None and os.path.exists(line_list_filename):
                self.import_linelist_as_profile_models(line_list_filename)

            # Construct any default spectral models.
            deconstructed_spectral_models = self.setting(("default_spectral_models", ))
            if deconstructed_spectral_models is not None:
                logger.warn("default_spectral_models is not implemented (skipping)")

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

        ### Some comments for clarity
        ## twd: a scratch directory for creatintg files that will be tar'd as the smh save file
        ## twd_paths: paths of files to save in tarball. Does not have to be actually point to twd.
        ## metadata["reconstruct_paths"]: names of files in twd needed to load the .smh file
        ## metadata["reconstruct_copied_paths"]: names of files kept with .smh file but not needed to load

        start0 = time.time()

        if not session_path.lower().endswith(".smh"):
            session_path = "{}.smh".format(session_path)

        if os.path.exists(session_path) and not overwrite:
            raise IOError("path '{}' already exists".format(session_path))

        metadata = self.metadata.copy()
        protocol = kwargs.pop("protocol", 2)

        # Create a temporary working directory and copy files over.
        twd = mkdtemp(**kwargs)
        twd_paths = [] + list(self._input_spectra_paths)

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
            raise IOError("This is an old session (version<0.2) with line_list (NOT SAVING)!"
                          "Running a conversion is required to save..")
        
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

        # Keep files copied to the temporary working directory (e.g. original linelists).
        twd_paths.extend(self.metadata["reconstruct_copied_paths"])
        
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
        logger.info("Saved file to {} ({:.1f}s)".format(session_path, time.time()-start0))

        # Remove the temporary working directory.
        rmtree(twd)

        if exception_occurred:
            raise


        return True


    def import_spectral_model_states(self, path):
        """
        Import list of spectral models from disk and append to current spectral models

        :param path:
            The disk location of the serialized transitions.
        """

        with open(path, 'rb') as fp:
            spectral_model_states = pickle.load(fp)
        spectral_models = self.reconstruct_spectral_models(spectral_model_states)
        self.metadata["spectral_models"].extend(spectral_models)
        return len(spectral_models)

    def export_spectral_model_states(self, path):
        # TODO implement mask saving etc.
        states = [_.__getstate__() for _ in self.spectral_models]
        with open(path, 'w') as fp:
            pickle.dump(states)
        return True

    def reconstruct_spectral_models(self, spectral_model_states):
        """
        When saving an SMH file or exporting its spectral models, we serialize 
        the spectral model classes into a state. 
        This function reconstructs the spectral models from that serialized state.
        """
        reconstructed_spectral_models = []
        for state in spectral_model_states:
            start = time.time()
            if "transitions" in state.keys():
                args = [self, LineList(state["transitions"])]
            else:
                assert "transition_hashes" in state.keys()
                raise IOError("Old spectral model format! (v<0.2)"
                              "(hashes instead of linelist) Cannot load")
            
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
            t2 = time.time()-start
            if t2 > 1.0:
                logger.debug("  Long time to load model {:.3f} {} {} {}\n".format(t2, len(model.transitions), model.elements, model.wavelength))
            
        return reconstructed_spectral_models

    def export_spectral_models(self, path, overwrite=False, keep_measurements=True):
        """
        Export list of spectral models to disk.
        
        :param path:
            The disk location to serialize transitions.
        :param overwrite:
            (default False) If True, overwrite path if it exists
        :param keep_measurements:
            (default True)
            If True, keep measurements specific to this star (equivalent widths, abundances).
            If False, remove that information and only keep atomic data, masks,
              fitting parameters, etc. [TODO make sure to include automasks]
        """
        if os.path.exists(path) and not overwrite:
            raise IOError("path '{}' already exists".format(path))
        raise NotImplementedError
        
    @classmethod
    def load(cls, session_path, skip_spectral_models=False, **kwargs):
        """
        Create a Session from a path saved to disk.

        :param session_path:
            The disk location where to load the session from.
        """

        start0 = time.time()

        # Extract all.
        tarball = tarfile.open(name=session_path, mode="r:gz")
        twd = mkdtemp(**kwargs)
        tarball.extractall(path=twd)

        # Reconstruct the session, starting with the initial paths.
        with open(os.path.join(twd, "session.pkl"), "rb") as fp:
            metadata = pickle.load(fp)

        # Load in the template spectrum.
        template_spectrum_path \
            = metadata["reconstruct_paths"].get("template_spectrum_path", None)
        if template_spectrum_path is not None:
            metadata["rv"]["template_spectrum_path"] \
                = os.path.join(twd, template_spectrum_path)

        # Create the object using the temporary working directory input spectra.
        session = cls([os.path.join(twd, basename) \
            for basename in metadata["reconstruct_paths"]["input_spectra"]],
            twd=twd, from_load=True)
        
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
        # Note that to serialize metadata, we converted the spectral models
        # into a serializable state
        start = time.time()
        spectral_model_states = session.metadata.get("spectral_models", [])
        reconstructed_spectral_models = session.reconstruct_spectral_models(spectral_model_states)
        logger.debug("Time to reconstruct {} spectral models: {:.3f}".format(
                len(reconstructed_spectral_models), time.time()-start))
        session.metadata["spectral_models"] = reconstructed_spectral_models

        # Clean up the TWD when Python exits.
        atexit.register(rmtree, twd)

        logger.info("Loaded file {} ({:.1f}s)".format(session_path, time.time()-start0))
        logger.debug("Input spectra paths: {}".format(session._input_spectra_paths))

        return session


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
    def spectral_models(self):
        """
        Shortcut for accessing spectral models
        """
        return self.metadata.get("spectral_models", [])

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

        if isinstance(key_tree, string_types):
            key_tree = [key_tree]
        
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

        transitions = []
        spectral_model_indices = []

        filtering = kwargs.pop("filtering",
            lambda model: model.use_for_stellar_parameter_inference)
        for i, model in enumerate(self.metadata["spectral_models"]):

            transitions.append(model.transitions[0])
            
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
        transitions = LineList.vstack(transitions)
        spectral_model_indices = np.array(spectral_model_indices)
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
        
        NOTE: there may be a bug where the wrong EQW gets passed if you specify spectral_models
        It hasn't seemed to be a huge problem.
        Not going to fix because it is made obsolete by SMHR v0.2
        """
        start = time.time()
        if spectral_models is None:
            spectral_models = self.metadata["spectral_models"]

        equivalent_widths = []
        equivalent_width_errs = []
        transitions = []
        spectral_model_indices = []
        
        num_profile = 0
        num_synth = 0
        for i,spectral_model in enumerate(spectral_models):
            if isinstance(spectral_model, ProfileFittingModel):
                spectral_model_indices.append(i)
                transitions.append(spectral_model.transitions[0])
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
                logger.info("Ignoring synthesis",spectral_model)
                num_synth += 1
            else:
                raise RuntimeError("Unknown model type: {}".format(type(spectral_model)))
            
        if num_profile > 0 and \
        (len(equivalent_widths) == 0 \
        or np.isfinite(equivalent_widths).sum() == 0):
            raise ValueError("no measured transitions to calculate abundances")
        
        transitions = LineList.vstack(transitions)
        spectral_model_indices = np.array(spectral_model_indices)
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
        logger.info("Time to measure {} abundances: {:.1f}".format(np.sum(finite), time.time()-start))
        return abundances, uncertainties if calculate_uncertainties else abundances

    def summarize_spectral_models(self, spectral_models=None, organize_by_element=False,
                                  use_weights = None, use_finite = True, what_fe = 1):
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
        :param what_fe:
            1 or 2 depending on Fe I or Fe II
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
                if what_fe == 1:
                    FeH = summary_dict[26.0][4]
                elif what_fe == 2:
                    FeH = summary_dict[26.1][4]
                else:
                    raise ValueError(str(what_fe))
        except KeyError:
            # Fe not measured yet
            FeH = np.nan
        for key in all_logeps:
            summary = summary_dict[key]
            summary[5] = summary[4] - FeH
            summary_dict[key] = summary

        total_num_models_summarized = np.sum([len(x) for x in all_logeps.values()])
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
            logger.warn("Must have normalized spectrum to make summary plot")
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
            logger.warn("Must have normalized spectrum to make summary plot")
            return None
        return smh_plotting.make_summary_plot(defaults["summary_figure_ncap"],
                                       self.normalized_spectrum, figure)

    def make_snr_plot(self, figure=None):
        return smh_plotting.make_snr_plot(self.normalized_spectrum, figure)

    def copy_file_to_working_directory(self, filename, twd=None):
        """
        Generic utility to copy files safely to twd
        :param filename:
            the file to copy to working directory
        :param twd:
            if None (default) use self.twd, otherwise use specified twd
        """
        if twd is None:
            twd = self.twd
        else:
            assert os.path.isdir(twd)
        basename = os.path.basename(filename)
        path_to_copy = os.path.join(twd, basename)
        while os.path.exists(path_to_copy):
            new_basename = ".".join([
                    utils.random_string(), basename])
            path_to_copy = os.path.join(twd, new_basename)
        copyfile(filename, path_to_copy)
        self.metadata["reconstruct_copied_paths"].append(path_to_copy)
        logger.info("Made copy of {} in {}".format(filename, path_to_copy))

    def import_linelist_as_profile_models(self, filename, import_equivalent_widths=False,
                                          copy_to_working_dir=True):
        """
        :param filename:
            path to a valid line list to be converted into (many) profile models
        :param import_equivalent_widths:
            (default False)
            if True, add equivalent width in line list to the profile model 
            and mark as acceptable (used for importing literature values)
        :param copy_to_working_dir:
            if True [default], copy the linelist to the session's working directory
        """
        assert os.path.exists(filename), filename

        if copy_to_working_dir:
            self.copy_file_to_working_directory(filename)

        start = time.time()
        line_list = LineList.read(filename, verbose=True)
        logger.debug("Time to load linelist {:.1f}".format(time.time()-start))
        
        if import_equivalent_widths:
            try:
                assert np.any(np.isfinite(line_list["equivalent_width"]))
            except:
                raise KeyError("no equivalent widths found in imported line list")
        
        start = time.time()
        spectral_models_to_add = []
        for i in range(len(line_list)):
            line = line_list[i]
            model = ProfileFittingModel(self, line)
            if import_equivalent_widths and np.isfinite(line["equivalent_width"]):
                model.metadata.update({
                        "is_acceptable": True,
                        "fitted_result": [None, None, {
                                # We assume supplied equivalent widths are in milliAngstroms
                                "equivalent_width": \
                                    (1e-3 * line["equivalent_width"], 0.0, 0.0),
                                "reduced_equivalent_width": \
                                    (-3+np.log10(line["equivalent_width"]/line["wavelength"]),
                                      0.0, 0.0)
                        }]
                })
            spectral_models_to_add.append(model)
        self.metadata["spectral_models"].extend(spectral_models_to_add)
        logger.debug("Created {} profile models in {:.1f}s".format(len(line_list),
                                                                   time.time()-start))
        return
        
    def import_linelist_as_synthesis_model(self, filename, elements,
                                           copy_to_working_dir=True, **kwargs):
        """
        :param filename:
            path to a valid line list to be converted into one synthesis model
        :param elements:
            elements to measure in this synthesis model
        :param copy_to_working_dir:
            if True [default], copy the linelist to the session's working directory
        kwargs are passed to SpectralSynthesisModel.__init__
        """
        assert os.path.exists(filename), filename

        if copy_to_working_dir:
            self.copy_file_to_working_directory(filename)

        start = time.time()
        line_list = LineList.read(filename, verbose=True)
        logger.debug("Time to load linelist {:.1f}".format(time.time()-start))
        
        start = time.time()
        spectral_model = SpectralSynthesisModel(self, line_list, elements, **kwargs)
        self.metadata["spectral_models"].append(spectral_model)
        logger.debug("Created synthesis model with {} lines in {:.1f}s".format(len(line_list),
                                                                               time.time()-start))
        return

    def import_master_list(self, filename,
                           copy_to_working_dir=True, **kwargs):
        """
        Use a "master" list to create a bunch of measurements
        wavelength, species, expot, loggf, type, filename(for syn or list)
        """
        assert os.path.exists(filename), filename

        if copy_to_working_dir:
            self.copy_file_to_working_directory(filename)
            
        master_list = ascii.read(filename, **kwargs).filled()
        logger.debug(master_list)
        types = np.array(map(lambda x: x.lower(), np.array(master_list["type"])))
        assert np.all(map(lambda x: (x=="eqw") or (x=="syn") or (x=="list"), types)), types

        num_added = 0

        ## Add EQW
        eqw = master_list[types=="eqw"]
        if len(eqw) > 0:
            line_list = LineList.create_basic_linelist(eqw["wavelength"],
                                                       eqw["species"],
                                                       eqw["expot"],
                                                       eqw["loggf"])
            spectral_models_to_add = []
            for i in range(len(line_list)):
                line = line_list[i]
                model = ProfileFittingModel(self, line)
                spectral_models_to_add.append(model)
            self.metadata["spectral_models"].extend(spectral_models_to_add)
            num_added += len(spectral_models_to_add)
        
        ## Add LIST
        lists = master_list[types=="list"]
        for row in lists:
            _filename = row["filename"]
            try:
                if _filename.endswith(".fits"):
                    ll = LineList.read(_filename, format='fits')
                else:
                    ll = LineList.read(_filename)
            except Exception as e:
                logger.warn("Could not import {}".format(_filename))
                logger.warn(e)
            else:
                spectral_models_to_add = []
                for i in range(len(ll)):
                    line = ll[i]
                    model = ProfileFittingModel(self, line)
                    spectral_models_to_add.append(model)
                self.metadata["spectral_models"].extend(spectral_models_to_add)
                num_added += len(ll)

        ## Add SYN
        syn = master_list[types=="syn"]
        for row in syn:
            elem1, elem2, isotope1, isotope2, ion = \
                utils.species_to_elems_isotopes_ion(row['species'])
            logger.debug("{} -> {} {} {} {} {}".format(row['species'],elem1,elem2,isotope1,isotope2,ion))
            if elem2=="":
                element = [elem1]
            else:
                # Just hardcoding for now, have to do something about this later
                if elem1 == "H" and elem2 == "C": element = ["C"]
                logger.debug("Hardcoded element: C")
            _filename = row["filename"]
            what_wavelength = row['wavelength']
            what_species = [row['species']]
            what_expot = row['expot']
            what_loggf = row['loggf']
            kwargs = {"what_wavelength":what_wavelength,
                      "what_species":what_species,
                      "what_expot":what_expot,
                      "what_loggf":what_loggf}
            try:
                self.import_linelist_as_synthesis_model(_filename, element,
                                                        copy_to_working_dir=copy_to_working_dir,
                                                        **kwargs)
            except Exception as e:
                logger.warn("Could not import {}".format(_filename))
                logger.warn(e)
            else:
                num_added += 1
        
        if num_added != len(master_list):
            logger.warn("Created {} models out of {} in the master list".format(
                    num_added, len(master_list)))
        return None
    
    def export_normalized_spectrum(self, path):
        """ Write out the normalized spectrum """
        self.normalized_spectrum.write(path)
        return None
    
    def export_unnormalized_spectrum(self, path):
        """ Coadd input orders into one spectrum and write out """
        new_spectrum = specutils.spectrum.coadd(self.input_spectra)
        new_spectrum.redshift(self.metadata["rv"]["rv_applied"])
        new_spectrum.write(path)
        return None
    
    def export_stitched_continuum(self, path):
        """ Coadd fitted continuums into one spectrum and write out """

        assert len(self.input_spectra) == len(self.metadata["normalization"]["continuum"])
        continuums = []
        for i, (spectrum, continuum) \
        in enumerate(zip(self.input_spectra,
        self.metadata["normalization"]["continuum"])):
            continuums.append(specutils.Spectrum1D(
                    spectrum.dispersion,continuum, 1./continuum))

        new_spectrum = specutils.spectrum.coadd(continuums)
        new_spectrum.redshift(self.metadata["rv"]["rv_applied"])
        new_spectrum.write(path)
        return None
    
