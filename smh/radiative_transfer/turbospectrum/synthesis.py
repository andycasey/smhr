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
    _turbospectrum_defaults = yaml.load(fp)





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

    path = utils.twd_path(twd=twd, **kwargs)


    # Calculate opacities.
    babsma_lu_kwds = kwds.copy()
    babsma_lu_kwds["opacities_path"] = path("opacities.in")
    babsma_lu_kwds["photosphere_path"] = path("photosphere.in")

    # Write atmosphere file.
    photosphere.write(babsma_lu_kwds["photosphere_path"], format="turbospectrum")

    # Get babsma_lu template and format it.
    with resource_stream(__name__, "babsma_lu.in") as fp:
        babsma_lu_template = fp.read()

    if abundances is None:
        abundances = {}

    # TODO Put the abundances on the right scale????

    formatted_abundances = "\n".join(["{0} {1:.2f}".format(species, abundance) \
        for species, abundance in abundances.items()])

    babsma_lu_kwds.update(
        num_abundances=len(abundances),
        formatted_abundances=formatted_abundances)



    # XI
    # num abundances
    # formatted abundances.

    # If no XI value, then estimate it.
    xi = photosphere.meta.get("microturbulence", None)
    if xi is None:
        xi = 1.0 # [km/s]
        logger.warn("No microturbulence given in photosphere. Assuming 1 km/s.")

    # TODO microturbulence may have been specified by kwargs....
    if "microturbulence" in babsma_lu_kwds:
        raise a

    babsma_lu_kwds["microturbulence"] = xi

    # Set some default values???

    babsma_lu_contents = babsma_lu_template.format(**babsma_lu_kwds)


    # Calcualte opacities.
    command = kwds["turbospectrum_babsma_lu"]

    os.symlink(kwds["turbospectrum_data_dir"], path("DATA"))
    proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        cwd=path(""))

    out, err = proc.communicate(input=babsma_lu_contents)
    errcode = proc.returncode


    # Calculate the spectrum.
    os.symlink(kwds["turbospectrum_molecules_dir"], path("molecules"))
    command = kwds["turbospectrum_bsyn_lu"]

    with resource_stream(__name__, "bsyn_lu.in") as fp:
        bsyn_lu_template = fp.read()

    if isotopes is None:
        isotopes = {}
    formatted_isotopes = "\n".join(["{0} {1:.2f}".format(mass_number, relative_abundance) \
        for mass_number, relative_abundance in isotopes.items()])


    babsma_lu_kwds.update(spectrum_path=path("spectrum.out"),
        num_isotopes=len(isotopes), formatted_isotopes=formatted_isotopes,
        is_spherical=["F", "T"][photosphere.meta["radius"] > 0])

    transitions_basename = "transitions.in"
    transitions.write(path(transitions_basename), format="turbospectrum")
        
    # TODO: MOlecule files.

    n_files = ["DATA/Hlinedata", transitions_basename]

    babsma_lu_kwds.update(formatted_n_files="\n".join(n_files),
        n_files=len(n_files))


    bsyn_lu_contents = bsyn_lu_template.format(**babsma_lu_kwds)
    command = kwds["turbospectrum_bsyn_lu"].split()


    proc_synth = subprocess.Popen(command, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out_synth, err_synth = proc_synth.communicate(input=bsyn_lu_contents)
    errcode_synth = proc_synth.returncode




    raise a
    # Get regions to synthesise.



    # Update abundances

    # Calculate opacities

    # Synthesise



