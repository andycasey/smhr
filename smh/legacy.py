#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import tarfile
import atexit
import os, sys, glob, time
from shutil import copyfile, rmtree
import pickle

from .utils import mkdtemp
from . import (Session, LineList, specutils, __version__)
from spectral_models import ProfileFittingModel, SpectralSynthesisModel

""" Functions for converting and loading older versions of SMHR """

def extract_filename_to_twd(filename, twd=None):
    """
    Extract the contents of filename
    """
    if twd is None:
        twd = mkdtemp()
        atexit.register(rmtree, twd)
    
    with tarfile.open(name=filename, mode="r:gz") as tarball:
        tarball.extractall(path=twd)
    
    return twd

def identify_smh_version(filename):
    """
    Attempts to identify the version that an SMH file is.
    """
    assert os.path.exists(filename)
    
    try:
        twd = extract_filename_to_twd(filename)
    except:
        raise RuntimeError("Cannot extract tarball, probably not an SMH file")
    
    # Try to find session.pkl, if not its' probably an old SMH file
    metadata_path = os.path.join(twd,"session.pkl")
    if not os.path.exists(metadata_path):
        raise RuntimeError("Old SMH?")

    with open(metadata_path, "rb") as fp:
        metadata = pickle.load(fp,encoding="latin1")
    if "VERSION" in metadata.keys():
        return session.metadata["VERSION"]

    line_list = metadata["reconstruct_paths"].get("line_list", None)
    if line_list is not None:
        linelist_path = os.path.join(twd, line_list)
        if os.path.exists(linelist_path):
            return "0.1"
    
    raise RuntimeError("Unknown SMHR type")

def convert_v0_1_to_v0_2(fname_in, fname_out, overwrite=False):
    """
    v0.1 is the last version to have session.metadata["line_list"].
    This converts to v0.2, which removes "line_list" and instead stores
      the atomic data with the spectral models.
    Note that if session.load() is changed, you will likely have to modify this code too.
    """
    assert os.path.exists(fname_in)
    if os.path.exists(fname_out) and not overwrite:
        raise IOError("{} already exists (set overwrite=True if needed)".format(fname_out))
    
    # Check file type
    version = identify_smh_version(fname_in)
    assert version=="0.1", version

    # Manually extract and get information from the file
    twd = extract_filename_to_twd(fname_in)
    metadata_path = os.path.join(twd, "session.pkl")
    with open(metadata_path, "rb") as fp:
        metadata = pickle.load(fp,encoding="latin1")

    # Create the object using the temporary working directory input spectra.
    session = Session([os.path.join(twd, basename) \
        for basename in metadata["reconstruct_paths"]["input_spectra"]])
    session.metadata.update(metadata)

    # Load in the linelist (but not into the session!)
    linelist_path = os.path.join(twd, metadata["reconstruct_paths"].get("line_list", None))
    linelist = LineList.read(linelist_path, format='fits')
    
    # Load in the template spectrum.
    template_spectrum_path \
        = metadata["reconstruct_paths"].get("template_spectrum_path", None)
    if template_spectrum_path is not None:
        metadata["rv"]["template_spectrum_path"] \
            = os.path.join(twd, template_spectrum_path)

    # Load in any normalized spectrum.
    normalized_spectrum \
        = metadata["reconstruct_paths"].get("normalized_spectrum", None)
    if normalized_spectrum is not None:
        session.normalized_spectrum = specutils.Spectrum1D.read(
            os.path.join(twd, normalized_spectrum))

    # Remove any reconstruction paths.
    metadata.pop("reconstruct_paths")

    
    # Now load as before: use hashes to index linelist
    try:
        hashes = linelist["hash"]
    except KeyError:
        print("LineList does not have hashes; computing again")
        hashes = linelist.compute_hashes()
    iisort = np.argsort(hashes)
    # turning hashes into transitions and add to the state
    spectral_model_states = metadata.get("spectral_models", [])
    for state in spectral_model_states:
        this_hashes = state["transition_hashes"]
        sorted = np.searchsorted(hashes, this_hashes, sorter=iisort)
        indices = iisort[sorted]
        transitions = linelist[indices]
        state["transitions"] = transitions
    # Then reconstruct as normal
    reconstructed_spectral_models = session.reconstruct_spectral_models(spectral_model_states)
    session.metadata["spectral_models"] = reconstructed_spectral_models

    # HACK: line_list_argsort_hashes has a bug and is really dumb not needed legacy
    if "line_list_argsort_hashes" in metadata.keys():
        _ = session.metadata.pop("line_list_argsort_hashes")

    # Create and save a new session
    session.save(fname_out, overwrite=True)
