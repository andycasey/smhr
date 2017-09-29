#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Functionality to derive equivalent widths from a transition's curve-of-growth using the
Turbospectrum radiative transfer code.
"""

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


with resource_stream(__name__, "defaults.yaml") as fp:
    _turbospectrum_defaults = yaml.load(fp)


def equivalent_widths(photosphere, transitions, full_output=False, verbose=False,
    twd=None, **kwargs):
    """
    Calculate equivalent widths from a curve-of-growth.

    :param photosphere:
        A formatted photosphere.

    :param transitions:
        A list of atomic transitions with measured equivalent widths.

    :param verbose: [optional]
        Specify verbose flags.
    """

    kwds = _turbospectrum_defaults.copy()
    kwds.update(kwargs)

    # Set default regions to synthesize.
    kwds.setdefault("dispersion_min", 
        min(transitions["wavelength"]) - kwds["opacity_contribution"])
    kwds.setdefault("dispersion_max",
        max(transitions["wavelength"]) + kwds["opacity_contribution"] \
                                       + kwds["dispersion_delta"])

    # Update keywords.
    xi = photosphere.meta.get("microturbulence", None) or 1.0
    kwds.update(
        microturbulence=xi,
        is_spherical=["F", "T"][photosphere.meta["radius"] > 0],
        num_isotopes=0, formatted_isotopes="",
        num_abundances=0, formatted_abundances="",
        abfind="false"
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

    if verbose:
        logger.info("Output from opacitices:\n{}".format(op_out))
        logger.warn("Error from opacities:\n{}".format(op_err))

    if op_proc.returncode:
        logging.exception(
            "Exception when calculating opacities in Turbospectrum: {}".format(
                op_err))
        raise ValueError("exception when calculating Turbospectrum opacities")

    transitions.write(path(kwds["transitions_path"]), format="turbospectrum")

    # TODO: Molecules
    n_files = [kwds["transitions_path"]]
    kwds.update(formatted_n_files="\n".join(n_files), n_files=len(n_files))

    # Calculate the abundances.
    proc = subprocess.Popen(["eqwidt_lu"], stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=path(""))

    with resource_stream(__name__, "eqwidt_lu.in") as fp:
        bsyn_contents = fp.read()

    out, err = proc.communicate(input=bsyn_contents.format(**kwds))
    if proc.returncode:
        logging.exception("Exception when calculating EWs in Turbospectrum:"\
            "\n{}".format(err))
        raise ValueError("exception when calculating EWs with Turbospectrum")

    transitions = transitions.copy()
    transitions.sort(["species", "wavelength"])
    transitions["equivalent_width"] = np.nan * np.ones(len(transitions))

    wavelengths, equivalent_widths = np.loadtxt(
        path(kwds["result_path"]), usecols=(2, 5, )).T

    indices = np.where(np.in1d(transitions["wavelength"], wavelengths))[0]
    transitions["equivalent_width"][indices] = equivalent_widths

    return transitions


