#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Synthesis functionality using the Turbospectrum radiative transfer code. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import numpy as np
import os
import subprocess
import yaml
from pkg_resources import resource_stream

from . import utils

logger = logging.getLogger(__name__)

# Load defaults.
with resource_stream(__name__, "defaults.yaml") as fp:
    try:
        _turbospectrum_defaults = yaml.load(fp, yaml.FullLoader)
    except AttributeError:
        _turbospectrum_defaults = yaml.load(fp)


def _is_supported():

    # Must have these things on path:
    required_executables_on_path = ("babsma_lu", "bsyn_lu")
    for each in required_executables_on_path:
        try:
            subprocess.check_output("which {}".format(each), shell=True)

        except subprocess.CalledProcessError:
            raise OSError("cannot find executable {} on path".format(each))

    # and environment variable $TURBODATA
    if os.environ.get("TURBODATA") is None \
    or not os.path.exists(os.environ.get("TURBODATA")):
        raise OSError("the TURBODATA environment variable must point to the "\
                      "DATA directory for turbospectrum")

    return True



def synthesize(photosphere, transitions, abundances=None, isotopes=None,
    verbose=False, twd=None, **kwargs):
    """
    Sythesize a stellar spectrum given the model photosphere and transitions
    provided. 

    Optional keyword arguments include:
    dispersion_min
    dispersion_max
    dispersion_delta
    opacity_contribution

    """

    if abundances is None:
        abundances = {}

    # TODO Put the abundances on the right scale????
    formatted_abundances = "\n".join(["{0} {1:.2f}".format(species, abundance) \
        for species, abundance in abundances.items()])

    if isotopes is None:
        isotopes = {}

    formatted_isotopes = "\n".join(["{0} {1:.2f}".format(mass_number, relative_abundance) \
        for mass_number, relative_abundance in isotopes.items()])

    kwds = _turbospectrum_defaults.copy()
    kwds.update(kwargs)

    # Set default regions to synthesize.
    kwds.setdefault("dispersion_min", 
        min(transitions["wavelength"]) - kwds["opacity_contribution"])
    kwds.setdefault("dispersion_max",
        max(transitions["wavelength"]) + kwds["opacity_contribution"] \
                                       + kwds["dispersion_delta"])

    # TODO put in abundances.
    if kwds["dispersion_max"] - kwds["dispersion_min"] > 1000.0:
        raise NotImplementedError("turbospectrum fails for large synth, "
                                  "and chunking not implemented yet")

    # Update keywords.
    xi = photosphere.meta.get("microturbulence", None) or 1.0
    kwds.update(
        microturbulence=xi,
        is_spherical=["F", "T"][photosphere.meta["radius"] > 0],
        num_isotopes=len(isotopes),
        formatted_isotopes=formatted_isotopes,
        num_abundances=len(abundances),
        formatted_abundances=formatted_abundances,
    )

    path = utils.twd_path(twd=twd, **kwargs)

    # Write atmosphere file.
    photosphere.write(path(kwds["photosphere_path"]), format="turbospectrum")

    # Requires environment variable for the Turbospectrum data path.
    os.symlink(os.environ.get("TURBODATA"), path("DATA"))

    # Calculate opacities.
    op_proc = subprocess.Popen(["babsma_lu"], stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        cwd=path(""))

    with resource_stream(__name__, "babsma_lu.in") as fp:
        babsma_contents = fp.read()

    op_out, op_err = op_proc.communicate(input=babsma_contents.format(**kwds))

    if op_proc.returncode:
        logging.exception(
            "Exception when calculating opacities in Turbospectrum: {}".format(
                op_err))
        raise ValueError("exception when calculating Turbospectrum opacities")


    # Calculate the spectrum.
    transitions.write(path(kwds["transitions_path"]), format="turbospectrum")
    
    # TODO: Molecules
    n_files = ["DATA/Hlinedata", kwds["transitions_path"]]
    kwds.update(formatted_n_files="\n".join(n_files), n_files=len(n_files))

    # Calculate the spectrum.
    proc_synth = subprocess.Popen(["bsyn_lu"], stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=path(""))

    with resource_stream(__name__, "bsyn_lu.in") as fp:
        bsyn_contents = fp.read()

    synth_out, synth_err = proc_synth.communicate(input=bsyn_contents.format(**kwds))
    if proc_synth.returncode:
        logging.exception("Exception when calculating spectrum in Turbospectrum:"\
            "\n{}".format(synth_err))
        raise ValueError("exception when calculating spectrum with Turbospectrum")

    dispersion, flux = np.loadtxt(path(kwds["result_path"]), usecols=(0, 1)).T

    meta = {}

    return [(dispersion, flux, meta)]
