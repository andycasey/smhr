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
from copy import deepcopy
import warnings

import astropy.table
from astropy.io import ascii
from .linelists import LineList
from .utils import mkdtemp
from . import (photospheres, radiative_transfer, specutils, isoutils, utils)
from .spectral_models import ProfileFittingModel, SpectralSynthesisModel
from smh.photospheres.abundances import asplund_2009 as solar_composition
from . import (smh_plotting, __version__)
from .optimize_stellar_params import optimize_stellar_parameters as run_optimize_stellar_parameters
from .optimize_stellar_params import optimize_feh as run_optimize_feh

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
                    self.setting("default_Teff",4500), # K
                "surface_gravity":
                    self.setting("default_logg",0.85),
                "metallicity":
                    self.setting("default_MH",-2.0), # Solar-scaled
                "microturbulence":
                    self.setting("default_vt",1.3), # km/s
                "alpha":
                    self.setting("default_aFe",0.4),
                "syserr_effective_temperature":
                    self.setting("default_Teff_syserr",150),
                "syserr_surface_gravity":
                    self.setting("default_logg_syserr",0.3),
                "syserr_metallicity":
                    self.setting("default_MH_syserr",0.0),
                "syserr_microturbulence":
                    self.setting("default_vt_syserr",0.2),
                "staterr_effective_temperature":
                    self.setting("default_Teff_staterr",10),
                "staterr_surface_gravity":
                    self.setting("default_logg_staterr",0.1),
                "staterr_metallicity":
                    self.setting("default_MH_staterr",0.1),
                "staterr_microturbulence":
                    self.setting("default_vt_staterr",0.1),
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
                    with open(os.path.join(twd, ".damaged"), "wb") as dfp:
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
        failed_paths = []
        for path in twd_paths:
            try:
                tarball.add(path, arcname=os.path.basename(path))

            except:
                logger.warn(
                    "Cannot save path '{}' to session:".format(path))

                logger.warn(
                    "Continuing to save session without this path "
                    "before raising the issue")
                failed_paths.append(path)
                exception_occurred = True
                continue

        tarball.close()
        logger.info("Saved file to {} ({:.1f}s)".format(session_path, time.time()-start0))

        # Remove the temporary working directory.
        rmtree(twd)

        if exception_occurred:
            logger.exception("Paths unable to save to session: {}".format(failed_paths))
            raise


        return True


    def import_spectral_model_states(self, path):
        """
        Import list of spectral models from disk and append to current spectral models

        :param path:
            The disk location of the serialized transitions.
        """

        with open(path, 'rb') as fp:
            spectral_model_states = pickle.load(fp,encoding="latin1")
        spectral_models = self.reconstruct_spectral_models(spectral_model_states)
        self.metadata["spectral_models"].extend(spectral_models)
        return len(spectral_models)

    def export_spectral_model_states(self, path):
        # TODO implement mask saving etc.
        states = [_.__getstate__() for _ in self.spectral_models]
        with open(path, 'wb') as fp:
            pickle.dump(states, fp)
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
            ## python 2/3 issue
            if "rt_abundances" in state["metadata"].keys():
                state["metadata"]["rt_abundances"] = utils._fix_bytes_dict(state["metadata"]["rt_abundances"])
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
            metadata = pickle.load(fp,encoding="latin1")

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

        # Python 2/3
        if "isotopes" in metadata:
            metadata["isotopes"] = utils._fix_bytes_dict(metadata["isotopes"])
        
        # Update the new session with the metadata.
        session.metadata = metadata
        # A hack to maintain backwards compatibility
        session.metadata["stellar_parameters"].setdefault("alpha", 0.4)
        session.metadata["stellar_parameters"].setdefault("syserr_effective_temperature", 150)
        session.metadata["stellar_parameters"].setdefault("syserr_surface_gravity", 0.3)
        session.metadata["stellar_parameters"].setdefault("syserr_metallicity", 0.0)
        session.metadata["stellar_parameters"].setdefault("syserr_microturbulence", 0.2)
        session.metadata["stellar_parameters"].setdefault("staterr_effective_temperature", 10)
        session.metadata["stellar_parameters"].setdefault("staterr_surface_gravity", 0.1)
        session.metadata["stellar_parameters"].setdefault("staterr_metallicity", 0.1)
        session.metadata["stellar_parameters"].setdefault("staterr_microturbulence", 0.1)
        
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
                try:
                    default = yaml.load(fp, yaml.FullLoader)
                except AttributeError:
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
            try:
                default = yaml.load(fp, yaml.FullLoader)
            except AttributeError:
                default = yaml.load(fp)

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
            v_helio, v_bary = specutils.motions.corrections_from_headers(\
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
            try:
                v_helio = v_helio.to("km/s").value
                v_bary = v_bary.to("km/s").value
            except:
                pass
                
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
        
        # -----------------------------------------------------------------
        # E. Holmbeck: calculate the bcv if it doesn't exist
        if "barycentric_correction" in self.metadata["rv"]:
            return
        
        names = self._input_spectra_paths
        spectrum = names[0]
        if len(names) > 1:
            for s in names:
                if 'red' in s:
                    spectrum = s
                    break

        try:
            from astropy.io import fits
            _, headers = fits.getdata(spectrum, header=True)
        except OSError as e:
            print("Failure to read FITS headers, will not have heliocentric/barycentric corrections")

        try:
            v_helio, v_bary = specutils.motions.corrections_from_headers(\
                headers)
        
        except Exception as e:
            logger.error(
                "Exception in calculating heliocentric and barycentric motions")
            logger.error(e)
            v_helio, v_bary = (np.nan, np.nan)

        else:
            try:
                v_helio = v_helio.to("km/s").value
                v_bary = v_bary.to("km/s").value
            except:
                pass

        self.metadata["rv"].update({
            # Measurements
            "rv_measured": rv,
            "heliocentric_correction": v_helio,
            "barycentric_correction": v_bary,
        })

        logging.info(
            "Heliocentric velocity correction: {0:.2f} km/s".format(v_helio))
        logging.info(
            "Barycentric velocity correction: {0:.2f} km/s".format(v_bary))
        # -----------------------------------------------------------------

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

    @property
    def stellar_parameters(self):
        """
        Return Teff, logg, vt, MH
        """
        return self.metadata["stellar_parameters"]["effective_temperature"], \
               self.metadata["stellar_parameters"]["surface_gravity"], \
               self.metadata["stellar_parameters"]["microturbulence"], \
               self.metadata["stellar_parameters"]["metallicity"]
    @property
    def stellar_parameters_staterr(self):
        """
        Return Teff, logg, vt, MH
        """
        return self.metadata["stellar_parameters"]["staterr_effective_temperature"], \
               self.metadata["stellar_parameters"]["staterr_surface_gravity"], \
               self.metadata["stellar_parameters"]["staterr_microturbulence"], \
               self.metadata["stellar_parameters"]["staterr_metallicity"]
    @property
    def stellar_parameters_syserr(self):
        """
        Return Teff, logg, vt, MH
        """
        return self.metadata["stellar_parameters"]["syserr_effective_temperature"], \
               self.metadata["stellar_parameters"]["syserr_surface_gravity"], \
               self.metadata["stellar_parameters"]["syserr_microturbulence"], \
               self.metadata["stellar_parameters"]["syserr_metallicity"]
    @property
    def stellar_parameters_err(self):
        """
        Return Teff, logg, vt, MH
        """
        e1,e2,e3,e4 = self.stellar_parameters_staterr
        s1,s2,s3,s4 = self.stellar_parameters_syserr
        return np.sqrt(e1**2 + s1**2), np.sqrt(e2**2 + s2**2), \
               np.sqrt(e3**2 + s3**2), np.sqrt(e4**2 + s4**2)

    def set_stellar_parameters(self, Teff, logg, vt, MH, alpha=None):
        """ 
        Set stellar parameters (Teff, logg, vt, MH[, alpha])
        """
        self.metadata["stellar_parameters"]["effective_temperature"] = Teff
        self.metadata["stellar_parameters"]["surface_gravity"] = logg
        self.metadata["stellar_parameters"]["microturbulence"] = vt
        self.metadata["stellar_parameters"]["metallicity"] = MH
        if alpha is not None:
            self.metadata["stellar_parameters"]["alpha"] = alpha
        
        return None
        
    def set_stellar_parameters_errors(self, sysstat, dTeff, dlogg, dvt, dMH):
        """ 
        Set stellar parameter errors (sys/stat, Teff, logg, vt, MH)
        """
        assert sysstat in ["sys","stat"], sysstat
        self.metadata["stellar_parameters"][sysstat+"err_effective_temperature"] = dTeff
        self.metadata["stellar_parameters"][sysstat+"err_surface_gravity"] = dlogg
        self.metadata["stellar_parameters"][sysstat+"err_microturbulence"] = dvt
        self.metadata["stellar_parameters"][sysstat+"err_metallicity"] = dMH
        return None
    
    def stellar_parameter_uncertainty_analysis(self, transitions=None,
                                               tolerances=[5,0.01,0.01],
                                               systematic_errors=None,
                                               expot_balance_species=26.0,
                                               ionization_balance_species_1=26.0,
                                               ionization_balance_species_2=26.1,
                                               rew_balance_species=26.0):
        """
        Performs an uncertainty analysis on the stellar parameters.
        Teff from varying Teff until the slope is +1 sigma off from its current value
          (sigma is slope fitting standard error on mean)
        logg from varying logg until difference between Fe 1 and Fe 2 is +1 sigma off from its current value
          (sigma is standard error added in quadrature)
        microturbulence is changing Teff until slope is +1 sigma off from its current value
        MH is maximum standard deviation of Fe 1 and Fe 2
        alpha has no error
        
        systematic_errors: if a 4-length value, then the order is
        [dTeff, dlogg, dvt, dMH]
        
        tolerances for accuracy in Teff, logg, vt can be specified with tolerances keyword
        (default 5, .01, .01)
        
        expot_balance_species, ionization_balance_species_1, ionization_balance_species_2,
        rew_balance_species: set by default to be 26.0, 26.0, 26.1, 26.0.
        Most likely you'll only change rew_balance_species to 26.1 as needed.
        """
        
        from scipy.optimize import fmin
        from scipy.stats import linregress
        start = time.time()
        saved_stellar_params = deepcopy(self.metadata["stellar_parameters"])
        
        # Use abundance errors in fit? Should just not use this by default for now
        if self.setting(("stellar_parameter_inference", 
            "use_abundance_uncertainties_in_line_fits"), False):
            yerr_column = "abundance_uncertainty"
        else:
            yerr_column = None
        
        if transitions is None:
            # TODO right now just assuming all acceptable Fe EQW lines
            transitions = []
            spectral_models = []
            eqws = []
            rews = []
            for model in self.metadata["spectral_models"]:
                if model.is_acceptable and \
                   isinstance(model, ProfileFittingModel) and \
                   model.elements[0] == "Fe":
                    transitions.append(model.transitions[0])
                    spectral_models.append(model)
                    eqws.append(model.equivalent_width)
                    rews.append(model.reduced_equivalent_width)
            transitions = LineList.vstack(transitions)
            eqws = np.array(eqws); rews = np.array(rews)
            transitions["equivalent_width"] = eqws
            transitions["reduced_equivalent_width"] = rews
            assert np.all(np.isfinite(eqws)), np.sum(~np.isfinite(eqws))
            assert np.all(eqws > .01), np.sum(eqws <= .01)
            abundances = self.rt.abundance_cog(self.stellar_photosphere, transitions, twd=self.twd)
            transitions["abundance"] = abundances
        else:
            transitions = transitions.copy()
            eqws = transitions["equivalent_width"]
            assert np.all(np.isfinite(eqws)), np.sum(~np.isfinite(eqws))
            assert np.all(eqws > .01), np.sum(eqws <= .01)
            for needed_col in ["abundance", "expot", "reduced_equivalent_width"]:
                assert needed_col in transitions.colnames, needed_col
        initial_slopes = utils.equilibrium_state(transitions, columns=("expot", "reduced_equivalent_width"),
                                                 ycolumn="abundance", yerr_column=yerr_column)
        

        # Put everything in a big try/except block so you don't accidentally overwrite current SPs
        # (Might be unnecessary)
        try:
            # MH: stdev of all Fe lines (including Fe I and II)
            self.metadata["stellar_parameters"] = deepcopy(saved_stellar_params)
            abundances = self.rt.abundance_cog(self.stellar_photosphere, transitions, twd=self.twd)
            mh_error = np.nanstd(abundances)
            
            # Teff
            try:
                self.metadata["stellar_parameters"] = deepcopy(saved_stellar_params)
                m, b, median, sigma, N = initial_slopes[expot_balance_species].get("expot", (np.nan,np.nan,np.nan,np.nan,0))
                logger.info("Finding error in Teff slope: {:.3f} +/- {:.3f}".format(m, sigma[1]))
                Teff = saved_stellar_params["effective_temperature"]
                Teff_error = np.nan
                m_target = m + sigma[1]
                def _calculate_teff_slope(teff):
                    self.metadata["stellar_parameters"]["effective_temperature"] = teff
                    abundances = self.rt.abundance_cog(self.stellar_photosphere, transitions, twd=self.twd)
                    transitions["abundance"] = abundances
                    m, b, median, sigma, N = utils.equilibrium_state(transitions, columns=("expot",), ycolumn="abundance",
                                                                     yerr_column=yerr_column)[expot_balance_species]["expot"]
                    logger.debug("Teff={:.0f} m={:.3f} m_target={:.3f}".format(teff,m,m_target))
                    return m
                minfn = lambda teff: (m_target - _calculate_teff_slope(teff[0]))**2
                Teff_positive_slope = fmin(minfn, [Teff], xtol=tolerances[0], ftol=0.00001)[0]
                Teff_error = int(round(np.abs(Teff_positive_slope - Teff)))
                self.metadata["stellar_parameters"] = deepcopy(saved_stellar_params)
                logger.info("Teff: {} + {}".format(Teff, Teff_error))
            except:
                logger.warn("Error calculating uncertainties in Teff")
                raise

            # logg
            try:
                self.metadata["stellar_parameters"] = deepcopy(saved_stellar_params)
                def _get_fe_values(abundances, transitions):
                    ii1 = transitions["species"] == ionization_balance_species_1
                    ii2 = transitions["species"] == ionization_balance_species_2
                    N1 = np.sum(ii1); N2 = np.sum(ii2)
                    ab1 = np.mean(abundances[ii1]); ab2 = np.mean(abundances[ii2])
                    std1 = np.std(abundances[ii1])/np.sqrt(N1)
                    std2 = np.std(abundances[ii2])/np.sqrt(N2)
                    return ab1, ab2, std1, std2
                abundances = self.rt.abundance_cog(self.stellar_photosphere, transitions, twd=self.twd)
                abFe1, abFe2, semFe1, semFe2 = _get_fe_values(abundances, transitions)
                dFe0 = abFe1 - abFe2
                logger.debug("Finding error in logg: Fe1={:.2f}+/-{:.2f}, Fe2={:.2f}+/-{:.2f}".format(
                        abFe1, semFe1, abFe2, semFe2))
                Fe_offset = np.sqrt(semFe1**2 + semFe2**2)
                logg = saved_stellar_params["surface_gravity"]
                logg_error = np.nan
                dFe_target = dFe0 + Fe_offset
                
                def _calculate_logg_dFe(logg):
                    self.metadata["stellar_parameters"]["surface_gravity"] = logg
                    abundances = self.rt.abundance_cog(self.stellar_photosphere, transitions, twd=self.twd)
                    ab1,ab2,sem1,sem2 = _get_fe_values(abundances, transitions)
                    logger.debug("logg={:.2f} dFe={:.3f} dFe_target={:.3f}".format(logg,ab1-ab2,dFe_target))
                    return ab1 - ab2
                minfn = lambda logg: (dFe_target - _calculate_logg_dFe(logg[0]))**2
                logg_positive_error = fmin(minfn, [logg], xtol=tolerances[1], ftol=0.00001)[0]
                logg_error = np.abs(logg_positive_error - logg)
                self.metadata["stellar_parameters"] = deepcopy(saved_stellar_params)
                logger.info("logg: {:.2f} + {:.2f}".format(logg, logg_error))
            except:
                logger.warn("Error calculating uncertainties in logg")
                raise

            # vt
            try:
                self.metadata["stellar_parameters"] = deepcopy(saved_stellar_params)
                m, b, median, sigma, N = initial_slopes[rew_balance_species].get("reduced_equivalent_width", (np.nan,np.nan,np.nan,np.nan,0))
                logger.info("Finding error in vt slope: {:.3f} +/- {:.3f}".format(m, sigma[1]))
                vt = saved_stellar_params["microturbulence"]
                vt_error = np.nan
                m_target = m + sigma[1]
                def _calculate_vt_slope(vt):
                    self.metadata["stellar_parameters"]["microturbulence"] = vt
                    abundances = self.rt.abundance_cog(self.stellar_photosphere, transitions, twd=self.twd)
                    transitions["abundance"] = abundances
                    m, b, median, sigma, N = utils.equilibrium_state(transitions, columns=("reduced_equivalent_width",),ycolumn="abundance",
                                                                     yerr_column=yerr_column)[rew_balance_species]["reduced_equivalent_width"]
                    logger.debug("vt={:.2f} m={:.3f} m_target={:.3f}".format(vt,m,m_target))
                    return m
                minfn = lambda vt: (m_target - _calculate_vt_slope(vt[0]))**2
                vt_positive_slope = fmin(minfn, [vt], xtol=tolerances[2], ftol=0.00001)[0]
                vt_error = round(np.abs(vt_positive_slope - vt),2)
                self.metadata["stellar_parameters"] = deepcopy(saved_stellar_params)
                logger.info("vt: {:.2f} + {:.2f}".format(vt, vt_error))
            except:
                logger.warn("Error calculating uncertainties in vt")
                raise

        except Exception as e:
            self.metadata["stellar_parameters"] = deepcopy(saved_stellar_params)
            raise
        else:
            self.metadata["stellar_parameters"] = deepcopy(saved_stellar_params)
            self.set_stellar_parameters_errors("stat", Teff_error, logg_error, mh_error, vt_error)
            logger.info("Uncertainties: dTeff={:.0f} dlogg={:.2f} dvt={:.2f} dmh={:.2f}, took {:.1f}s".format(
                    Teff_error, logg_error, vt_error, mh_error, time.time()-start))
            if systematic_errors is not None:
                if len(systematic_errors) != 4:
                    logger.warn("Systematic errors are not length 4! {}".format(systematic_errors))
                    logger.warn("Skipping...")
                else:
                    self.set_stellar_parameters_errors("sys",
                                                       systematic_errors[0], systematic_errors[1],
                                                       systematic_errors[2], systematic_errors[3])
        return Teff_error, logg_error, vt_error, mh_error

    """
    Abundance uncertainty analysis.
    For every measurement, there is sigma_r and sigma_s (random and systematic error)
    Equivalent widths (fast, but not uniform):
        sigma_r = spectral_model.abundance_uncertainties[0]
                  Stored in every spectral model individually
        sigma_s = session.propagate_stellar_parameter_error_to_equivalent_widths
                  Stored in session.metadata["abundance_uncertainies_EQWsys"]
    Equivalent widths (slow, but same interface and a bit more flexibility):
        sigma_r = spectral_model.find_error(sigma=1)
                  Stored in spectral_model.metadata["1_sigma_abundance_error"]
        sigma_s = spectral_model.propagate_stellar_parameter_error()
                  Stored in spectral_model.metadata["systematic_abundance_error"]
    Synthesis:
        sigma_r = spectral_model.find_error(sigma=1)
                  stored in spectral_model.metadata["1_sigma_abundance_error"]
        sigma_s = spectral_model.propagate_stellar_parameter_error()
                  stored in spectral_model.metadata["systematic_abundance_error"]
    """
    def compute_all_abundance_uncertainties(self, print_memory_usage=False):
        """
        Call model.find_error() and model.propagate_stellar_parameter_error() for every
        acceptable spectral model that is not an upper limit.
        Tabulate and return all the results.
        """
        self.measure_abundances()
        data = np.zeros((len(self.spectral_models), 8)) + np.nan
        data[:,0] = -1
        logger.info("Starting abundance uncertainty loop")
        start = time.time()
        saved_stellar_params = self.metadata["stellar_parameters"].copy()
        for i, model in enumerate(self.spectral_models):
            if not model.is_acceptable: continue
            if model.is_upper_limit: continue
            # Save the data
            wavelength = model.wavelength
            species = np.ravel(model.species)[0]
            try:
                logeps = model.abundances[0]
            except:
                logeps = np.nan
            whattype = 0 if isinstance(model, ProfileFittingModel) else 1
            # Compute uncertainty
            try:
                staterr = model.find_error()
            except:
                staterr = np.nan
            try:
                syserr = model.propagate_stellar_parameter_error()
            except:
                syserr = np.nan
            
            data[i,0] = i
            data[i,1] = wavelength
            data[i,2] = species
            data[i,3] = whattype
            data[i,4] = logeps
            data[i,5] = staterr
            data[i,6] = syserr
            data[i,7] = np.sqrt(staterr**2 + syserr**2)

            if print_memory_usage:
                import psutil
                import resource
                process = psutil.Process(os.getpid())
                print(species, wavelength, process.memory_info().rss)  # in bytes 
                print("    ",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)  # in bytes 

        self.metadata["all_line_data"] = data
        logger.info("Finished abundance uncertainty loop in {:.1f}s".format(time.time()-start))
        return data
    def propagate_stellar_parameter_error_to_equivalent_widths(self, Teff_error, logg_error, vt_error, mh_error,
                                                               output_array=True, save_to_session=True):
        """
        Measure EQW abundances at multiple different stellar parameters
        Output:
        if output_array==True (default):
          An N x 10 array, where N is the number of species, and the columns are:
            species, +Teff, -Teff, +logg, -logg, +vt, -vt, +MH, -MH, total
        else:
          Dict[stellar_param] -> [+err, -err]
          Each err dict has two types of keys:
              species -> abundance difference for that species (logeps)
              (species1, species2) -> abundance difference for the ratio between those species
        """
        spnames = ["effective_temperature","surface_gravity","microturbulence","metallicity"]
        eqw_models = []
        for model in self.spectral_models:
            if model.is_acceptable and isinstance(model, ProfileFittingModel) and (not model.is_upper_limit):
                eqw_models.append(model)
        
        all_species = np.unique(list(map(lambda m: m.species[0], eqw_models)))
        all_ratios = []
        for _s1 in all_species:
            for _s2 in all_species:
                all_ratios.append((_s1, _s2))

        abund0 = self.measure_abundances(eqw_models, save_abundances=True, calculate_uncertainties=False)
        sp0 = self.metadata["stellar_parameters"].copy()
        summary0 = self.summarize_spectral_models()
        ratios0 = {}
        for ratio in all_ratios:
            ratios0[ratio] = summary0[ratio[0]][1] - summary0[ratio[1]][1]
        
        error_output = {}
        error_output["stellar_parameter_errors"] = [Teff_error, logg_error, vt_error, mh_error]
        error_output["orig_abund"] = abund0
        error_output["orig_ratios"] = ratios0
        for spname, error in zip(spnames,[Teff_error, logg_error, vt_error, mh_error]):
            sperr_list = []
            for sign in [+1., -1.]:
                abund_diffs = {}
                self.metadata["stellar_parameters"].update(sp0)
                self.metadata["stellar_parameters"][spname] += sign * error
                logger.debug(spname+" "+str(self.metadata["stellar_parameters"][spname]))
                self.measure_abundances(eqw_models, save_abundances=True, calculate_uncertainties=False)
                summary = self.summarize_spectral_models()
                for species in all_species:
                    logger.debug("{:4.1f} {:5.2f} {:5.2f}".format(species, summary[species][1], summary0[species][1]))
                    abund_diffs[species] = summary[species][1] - summary0[species][1]
                for ratio in all_ratios:
                    abund_diffs[ratio] = (summary[ratio[0]][1] - summary[ratio[1]][1]) - ratios0[ratio]
                sperr_list.append(abund_diffs)
            error_output[spname] = sperr_list
        # Reset abundance measurements in session
        self.measure_abundances(eqw_models, save_abundances=True, calculate_uncertainties=True)
        if output_array:
            def simplify_sperr_output(error_output):
                all_species = []
                for x in error_output["effective_temperature"][0].keys():
                    if isinstance(x, float): all_species.append(x)
                all_species = np.sort(all_species)
                spnames = ["effective_temperature","surface_gravity","microturbulence","metallicity"]
                fmt = "{:5.1f},{:5.2f},{:5.2f},{:5.2f},{:5.2f},{:5.2f},{:5.2f},{:5.2f},{:5.2f},{:5.2f}"
                data = np.zeros((len(all_species), 10))
                for i, species in enumerate(all_species):
                    data[i,0] = species
                    for j0, spname in enumerate(spnames):
                        sperr_list = error_output[spname]
                        for k0, (sign, sperr) in enumerate(zip([+1,-1], sperr_list)):
                            icol = 1 + 2*j0 + k0
                            data[i,icol] = sperr[species]
                        data[i,9] += (np.max(np.abs([data[i,1+2*j0],data[i,2*j0+1]])))**2
                    data[i,9] = np.sqrt(data[i,9])
                for row in data:
                    print(fmt.format(*row))
                return data
            error_output = simplify_sperr_output(error_output)
        if save_to_session:
            self.metadata["abundance_uncertainies_EQWsys"] = error_output
        return error_output

    def stellar_parameter_state(self, full_output=False, **kwargs):
        """
        Calculate the abundances of all spectral models that are used in the
        determination of stellar parameters.

        TODO this should be removed in favor of the measure_abundances method.
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



    def _spectral_models_to_transitions(self, spectral_models):
        equivalent_widths = []
        equivalent_width_errs = []
        transitions = []
        spectral_model_indices = []
        abundances = []
        for i,spectral_model in enumerate(spectral_models):
            if isinstance(spectral_model, ProfileFittingModel):
                spectral_model_indices.append(i)
                transitions.append(spectral_model.transitions[0])
                if spectral_model.is_acceptable:
                    try:
                        equivalent_widths.append(1000.* \
                            spectral_model.metadata["fitted_result"][-1]["equivalent_width"][0])
                    except:
                        equivalent_widths.append(np.nan)
                    try:
                        equivalent_width_errs.append(1000.* \
                            np.nanmax(spectral_model.metadata["fitted_result"][-1]["equivalent_width"][1:3]))
                    except:
                        equivalent_width_errs.append(np.nan)
                    try:
                        abundances.append(spectral_model.metadata["fitted_result"][-1]["abundances"][0])
                    except:
                        abundances.append(np.nan)
                else:
                    equivalent_widths.append(np.nan)
                    equivalent_width_errs.append(np.nan)
                    abundances.append(np.nan)
            elif isinstance(spectral_model, SpectralSynthesisModel):
                #logger.info("Ignoring synthesis",spectral_model)
                pass
            else:
                raise RuntimeError("Unknown model type: {}".format(type(spectral_model)))
        transitions = LineList.vstack(transitions)
        spectral_model_indices = np.array(spectral_model_indices)
        transitions["equivalent_width"] = equivalent_widths
        transitions["abundance"] = abundances
        return transitions

    def optimize_stellar_parameters(self, **kwargs):
        """
        Optimize the stellar parameters for this star using the spectral models
        associated with stellar parameter inference.
        """
        
        Teff, logg, vt, MH = self.stellar_parameters
        alpha = self.metadata["stellar_parameters"]["alpha"]
        initial_guess = [Teff, vt, logg, MH] # stupid me did not change the ordering to match
        logger.info("Initializing optimization at Teff={:.0f} logg={:.2f} vt={:.2f} MH={:.2f}".format(
            Teff, logg, vt, MH))
        
        # get the list of relevant spectral models.
        stellar_parameter_spectral_models = []
        for spectral_model in self.spectral_models:
            if spectral_model.use_for_stellar_parameter_inference and spectral_model.is_acceptable:
                if isinstance(spectral_model, SpectralSynthesisModel):
                    raise NotImplementedError("Syntheses cannot be used for stellar parameter state (check the Transitions manager)")
                stellar_parameter_spectral_models.append(spectral_model)
        transitions = self._spectral_models_to_transitions(stellar_parameter_spectral_models)
        finite = np.logical_and(np.isfinite(transitions["equivalent_width"]), transitions["equivalent_width"] > 0.01)
        if finite.sum() != len(finite):
            logger.warn("Number of finite transitions ({}) != number of transitions ({})".format(
                finite.sum(), len(finite)))
        transitions = transitions[finite]
        logger.info("Optimizing with {} transitions".format(len(transitions)))
        
        # interpolator, do obj. function
        out = run_optimize_stellar_parameters(initial_guess, transitions, **kwargs)
        final_parameters = out[1]
        new_Teff, new_vt, new_logg, new_MH = final_parameters
        
        if not out[0]:
            logger.warn("Optimization did not converge, stopped at Teff={:.0f} logg={:.2f} vt={:.2f} MH={:.2f}".format(
                new_Teff, new_logg, new_vt, new_MH))
            
        
        self.set_stellar_parameters(new_Teff, new_logg, new_vt, new_MH)
        return None

    
    # E. Holmbeck added this function to call the other "Solve" function
    def optimize_feh(self, params_to_optimize, **kwargs):
        """
        Optimize the stellar parameters for this star using the spectral models
        associated with stellar parameter inference.
        """
        
        Teff, logg, vt, MH = self.stellar_parameters
        initial_guess = np.array([Teff, vt, logg, MH]) # stupid me did not change the ordering to match
        logger.info("Initializing optimization at Teff={:.0f} logg={:.2f} vt={:.2f} MH={:.2f}".format(
            Teff, logg, vt, MH))
        
        # get the list of relevant spectral models.
        stellar_parameter_spectral_models = []
        for spectral_model in self.spectral_models:
            if spectral_model.use_for_stellar_parameter_inference and spectral_model.is_acceptable:
                if isinstance(spectral_model, SpectralSynthesisModel):
                    raise NotImplementedError("Syntheses cannot be used for stellar parameter state (check the Transitions manager)")
                stellar_parameter_spectral_models.append(spectral_model)
        transitions = self._spectral_models_to_transitions(stellar_parameter_spectral_models)
        finite = np.logical_and(np.isfinite(transitions["equivalent_width"]), transitions["equivalent_width"] > 0.01)
        if finite.sum() != len(finite):
            logger.warn("Number of finite transitions ({}) != number of transitions ({})".format(
                finite.sum(), len(finite)))
        transitions = transitions[finite]
        logger.info("Optimizing with {} transitions".format(len(transitions)))
        
        # interpolator, do obj. function
        out = run_optimize_feh(initial_guess, transitions, params_to_optimize, **kwargs)
        final_parameters = out[1]
        new_Teff, new_vt, new_logg, new_MH = final_parameters
        
        if not out[0]:
            logger.warn("Optimization did not converge, stopped at Teff={:.0f} logg={:.2f} vt={:.2f} MH={:.2f}".format(
                new_Teff, new_logg, new_vt, new_MH))
            
        
        self.set_stellar_parameters(new_Teff, new_logg, new_vt, new_MH)
        return None

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
                #logger.info("Ignoring synthesis",spectral_model)
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
                                  use_weights = None, use_finite = True, what_fe = 1,
                                  default_error = 0.1):
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
        def update_one_key(key,logeps):
            if key not in all_logeps:
                all_logeps[key] = []
            all_logeps[key].append(logeps)
            return
        # TODO abundance uncertainties too
        all_weights = {}
        def update_one_weight(key,err):
            if key not in all_weights:
                all_weights[key] = []
            all_weights[key].append(err**-2.)
            return
        for spectral_model in spectral_models:
            if not spectral_model.is_acceptable or spectral_model.is_upper_limit: continue
            
            abundances = spectral_model.abundances
            if abundances is None: continue
            # Try to get the abundance uncertainties
            try:
                errval = spectral_model.abundance_uncertainties or default_error
            except:
                errval = default_error
            errors = [errval for _ in abundances]
            
            for elem,species,logeps,logepserr in zip(spectral_model.elements, \
                                                     spectral_model.species, \
                                                     abundances, errors):
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
                            update_one_weight(species, logepserr)
                        elif isinstance(species,list):
                            for _species, _logeps, _logepserr in zip(species,logeps,logepserr):
                                update_one_key(_species, _logeps)
                                update_one_weight(_species, _logepserr)
                        else:
                            raise TypeError("key {} is of type {} (organized by {})".format(\
                                    species,type(species),what_key_type))
                elif isinstance(key,float):
                    assert not organize_by_element
                    update_one_key(key, logeps)
                    update_one_weight(key, logepserr)
                elif isinstance(key,str):
                    assert organize_by_element
                    update_one_key(key, logeps)
                    update_one_weight(key, logepserr)
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

            if use_weights:
                weights = np.array(all_weights[key])
                logeps = np.sum(weights*logepss)/np.sum(weights)
                stdev = np.sqrt(np.sum(weights*(logepss-logeps)**2)/np.sum(weights))
                stderr= stdev/np.sqrt(num_models)
            else:
                logeps = np.mean(logepss)
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
    
    def export_abundance_table(self, filepath, use_weights=False):
        ## TODO: put in upper limits too.
        summary_dict = self.summarize_spectral_models(use_weights=use_weights)
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
        # Erika added EW sigma to output
        linedata = np.zeros((len(spectral_models), 8)) + np.nan
        for i,spectral_model in enumerate(spectral_models):
            # TODO include upper limits
            if not spectral_model.is_acceptable or spectral_model.is_upper_limit: continue
            if isinstance(spectral_model, SpectralSynthesisModel):
                assert len(spectral_model.elements) == 1, spectral_model.elements
                wavelength = spectral_model.wavelength
                species = spectral_model.species[0][0]
                expot = spectral_model.expot
                loggf = spectral_model.loggf
                EW = np.nan
                e_EW = np.nan
                logeps = spectral_model.abundances[0]
                try:
                    logeps_err = spectral_model.metadata["2_sigma_abundance_error"]/2.0
                except:
                    logeps_err = np.nan
            elif isinstance(spectral_model, ProfileFittingModel):
                line = spectral_model.transitions[0]
                wavelength = line['wavelength']
                species = line['species']
                expot = line['expot']
                loggf = line['loggf']

                try:
                    EW = 1000.*spectral_model.metadata["fitted_result"][2]["equivalent_width"][0]
                    # Erika added EW sigma to output
                    e_EW = max(1000.*np.abs(spectral_model.metadata["fitted_result"][2]["equivalent_width"][1:]))
                    logeps = spectral_model.abundances[0]
                    logeps_err = spectral_model.abundance_uncertainties or np.nan
                except Exception as e:
                    print(e)
                    EW = np.nan
                    e_EW = np.nan
                    logeps = np.nan
                    logeps_err = np.nan
                if EW is None: EW = np.nan
                if logeps is None: logeps = np.nan
            else:
                raise NotImplementedError
            # Erika added EW sigma to output
            linedata[i,:] = [species, wavelength, expot, loggf, EW, e_EW, logeps, logeps_err]
        #ii_bad = np.logical_or(np.isnan(linedata[:,5]), np.isnan(linedata[:,4]))
        ii_bad = np.isnan(linedata[:,5])
        linedata = linedata[~ii_bad,:]
        if len(linedata) == 0:
            raise RuntimeError("No lines have abundances measured!")

        if filepath.endswith(".tex"):
            self._export_latex_measurement_table(filepath, linedata)
        else:
            self._export_ascii_measurement_table(filepath, linedata)
        logger.info("Exported to {}".format(filepath))
        return None
    def _export_latex_measurement_table(self, filepath, linedata):
        raise NotImplementedError
    def _export_ascii_measurement_table(self, filepath, linedata):
        # Erika added EW sigma to output
        names = ["species", "wavelength", "expot", "loggf", "EW", "e_EW", "logeps", "e_logeps"]
        tab = astropy.table.Table(linedata, names=names)
        tab.sort(["species","wavelength","expot"])
        tab["wavelength"].format = ".3f"
        tab["expot"].format = "5.3f"
        tab["loggf"].format = "6.3f"
        tab["EW"].format = "6.2f"
        tab["e_EW"].format = "6.2f"
        tab["logeps"].format = "6.3f"
        tab["e_logeps"].format = "6.3f"
        tab.write(filepath, format="ascii.fixed_width_two_line")
        return True

    def make_summary_plot(self, figure=None):
        with open(self._default_settings_path, "rb") as fp:
            try:
                defaults = yaml.load(fp, yaml.FullLoader)
            except AttributeError:
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
            try:
                defaults = yaml.load(fp, yaml.FullLoader)
            except AttributeError:
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
                                          copy_to_working_dir=False):
        """
        :param filename:
            path to a valid line list to be converted into (many) profile models
        :param import_equivalent_widths:
            (default False)
            if True, add equivalent width in line list to the profile model 
            and mark as acceptable (used for importing literature values)
        :param copy_to_working_dir:
            if True, copy the linelist to the session's working directory
        """
        assert os.path.exists(filename), filename

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
        
        if copy_to_working_dir:
            self.copy_file_to_working_directory(filename)

        return
        
    def import_linelist_as_synthesis_model(self, filename, elements,
                                           copy_to_working_dir=False, **kwargs):
        """
        :param filename:
            path to a valid line list to be converted into one synthesis model
        :param elements:
            elements to measure in this synthesis model
        :param copy_to_working_dir:
            if True, copy the linelist to the session's working directory
        kwargs are passed to SpectralSynthesisModel.__init__
        """
        assert os.path.exists(filename), filename

        start = time.time()
        line_list = LineList.read(filename, verbose=True)
        logger.debug("Time to load linelist {:.1f}".format(time.time()-start))
        
        start = time.time()
        spectral_model = SpectralSynthesisModel(self, line_list, elements, **kwargs)
        self.metadata["spectral_models"].append(spectral_model)
        logger.debug("Created synthesis model with {} lines in {:.1f}s".format(len(line_list),
                                                                               time.time()-start))

        if copy_to_working_dir:
            self.copy_file_to_working_directory(filename)

        return

    def import_master_list(self, filename,
                           copy_to_working_dir=False, **kwargs):
        """
        Use a "master" list to create a bunch of measurements
        wavelength, species, expot, loggf, type, filename(for syn or list)
        """
        assert os.path.exists(filename), filename

        master_list = ascii.read(filename, **kwargs).filled()
        logger.debug(master_list)
        print(master_list["type"])
        types = np.array(list(map(lambda x: x.lower(), np.array(master_list["type"]))))
        print(types)
        assert np.all([(x=="eqw") or (x=="syn") or (x=="list") for x in types]), types

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
                if copy_to_working_dir:
                    self.copy_file_to_working_directory(_filename)

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
                if elem1 == "H" and elem2 == "C":
                    element = ["C"]
                    logger.debug("Hardcoded element: CH->H")
                if elem1 == "C" and elem2 == "N":
                    element = ["N"]
                    logger.debug("Hardcoded element: CN->N")
                if elem1 == "H" and elem2 == "N":
                    element = ["N"]
                    logger.debug("Hardcoded element: NH->N")
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
                if copy_to_working_dir:
                    self.copy_file_to_working_directory(_filename)
        
        if num_added != len(master_list):
            logger.warn("Created {} models out of {} in the master list".format(
                    num_added, len(master_list)))

        if copy_to_working_dir:
            self.copy_file_to_working_directory(filename)
            
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
    
    def get_spectral_models_species_dict(self):
        all_models = {}
        for model in self.spectral_models:
            if isinstance(model, ProfileFittingModel): species = round(model.species[0],1)
            if isinstance(model, SpectralSynthesisModel):
                if len(model.species) > 1:
                    warnings.warn("Spectral model has multiple species: {}".format(model.species))
                species = round(model.species[0][0],1)
            # Use setdefault instead of get, not sure why it has to be this way
            species_models = all_models.setdefault(species, [])
            species_models.append(model)
        return all_models
    
    def initialize_normalization(self):
        N = len(self.input_spectra)
        self.metadata["normalization"] = {
            "continuum": [None] * N,
            "normalization_kwargs": [{}] * N
        }

    def initialize_rv(self):
        """
        Set things up so the RV GUI doesn't error out
        """
        wavelength_region = self.setting(("rv", "wavelength_regions"))
        resample = self.setting(("rv", "resample"))
        apodize = self.setting(("rv", "apodize"))
        normalization_kwargs = self.setting(("rv", "normalization"))
        
        template_spectrum_path = self.setting(("rv", "template_spectrum"))
        template_spectrum = specutils.Spectrum1D.read(template_spectrum_path)
        
        self.metadata["rv"].update({
            # Input settings
            "template_spectrum_path": template_spectrum_path,
            "template_spectrum": template_spectrum,
            "wavelength_region": wavelength_region,
            "resample": resample,
            "apodize": apodize,
            "normalization": normalization_kwargs.copy()
        })
        
