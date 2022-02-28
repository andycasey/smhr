#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Synthesis functionality with the MOOG(SILENT) radiative transfer code. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)


import logging
import numpy as np
import yaml
from pkg_resources import resource_stream

from . import utils

logger = logging.getLogger(__name__)

# Load the MOOG defaults.
with resource_stream(__name__, "defaults.yaml") as fp:
    try:
        _moog_defaults = yaml.load(fp, yaml.FullLoader)
    except AttributeError:
        _moog_defaults = yaml.load(fp)


def synthesize(photosphere, transitions, abundances=None, isotopes=None,
    verbose=False, twd=None, **kwargs):
    """
    Sythesize a stellar spectrum given the model photosphere and list of
    transitions provided. This wraps the MOOG `synth` driver.

    :param photosphere:
        A formatted photosphere.

    :param transitions:
        A list of atomic transitions with measured equivalent widths.

    :param abundances: [optional]
        The non-scaled-Solar abundances to use in the synthesis.

    :param isotopes: [optional]
        The isotopic fractions to use in the synthesis.

    :param verbose: [optional]
        Specify verbose flags to MOOG. This is primarily used for debugging.
    """

    # Create a temporary directory and write out the photoshere and transitions.
    path = utils.twd_path(twd=twd,**kwargs)
    model_in, lines_in = path("model.in"), path("lines.in")
    photosphere.write(model_in, format="moog")
    transitions.write(lines_in, format="moog")
    
    # Load the synth driver template.
    with resource_stream(__name__, "synth.in") as fp:
        template = fp.read().decode("utf-8")

    # Not these are SMH defaults, not MOOG defaults.
    kwds = _moog_defaults.copy()
    if verbose:
        kwds.update({
            "atmosphere": 2,
            "molecules": 2,
            "lines": 3, # 4 is max verbosity, but MOOG falls over.
        })

    # Abundances.
    # These are given to us as log_epsilon but MOOG wants them relative to the
    # photospheric metallicity.
    mh = photosphere.meta["stellar_parameters"]["metallicity"]
    abundances_formatted, num_synth = utils._format_abundances(
        abundances, subtract_solar=True, subtract_metallicity=mh)

    kwds["abundances_formatted"] = abundances_formatted

    # Isotopes.
    kwds["isotopes_formatted"] = utils._format_isotopes(isotopes, 
        kwargs.pop("isotope_ionisation_states", (0, 1)), num_synth=num_synth)

    # Parse keyword arguments.
    kwds.update(kwargs)

    # Update with default synthesis limits.
    # Note that opacity_contribution comes from defaults.yaml and can be
    # overwritten by the kwargs to this function.
    kwds.setdefault("dispersion_min", min(transitions["wavelength"])
        - kwds["opacity_contribution"])
    kwds.setdefault("dispersion_max", max(transitions["wavelength"]) \
        + kwds["opacity_contribution"] + kwds["dispersion_delta"])

    # Parse I/O files (these must be overwritten for us to run things.)
    kwds.update({
        "standard_out": path("synth.std.out"),
        "summary_out": path("synth.sum.out"),
        "model_in": model_in,
        "lines_in": lines_in,
    })

    # Put this into a while loop only in case we have to iteratively check for
    # edge effects due to syn_contribute
    while True:
        contents = template.format(**kwds)

        # Write this to batch.par and execute.
        moog_in = path("batch.par")
        with open(moog_in, "w") as fp:
            fp.write(contents)

        # Execute MOOG in the TWD.
        code, out, err = utils.moogsilent(moog_in, **kwargs)

        # Returned normally?
        assert code == 0 # HACK # TODO

        # Parse the output.
        spectra = _parse_synth_summary(kwds["summary_out"])

        # TODO: Check for physically unrealistic intensity jumps due to the 
        # value of opacity_contribution being too low.
        """
        for meta, dispersion, intensity in spectra:
            _check_for_unphysical_intensity_jumps(transitions, dispersion,
                intensity, kwds["opacity_contribution"])
        """

        break

    return spectra


def _parse_single_spectrum(lines):
    """
    Parse the header, dispersion and depth information for a single spectrum 
    that has been output to a summary synthesis file.

    :param lines:
        A list of string lines partitioned from a MOOG summary output file.

    :returns:
        A three-length tuple containing the (1) header rows, (2) an array of
        dispersion values, and (3) an array of intensity values.
    """

    # Skip over the header information
    for i, line in enumerate(lines):
        if line.startswith("MODEL:"): break

    else:
        raise ValueError("could not find model information for spectrum")

    # Get the dispersion information.
    start, end, delta, _ = np.array(lines[i + 1].strip().split(), dtype=float)

    # If MOOG doesn't have opacity contributions at a given wavelength, it will
    # just spit out ****** entries.
    #_pre_format = lambda l: l.strip().replace("*******", " 0.0000").replace("-0.0000"," 0.0000").split()
    def _pre_format(l):
        l = l.replace("*******", " 0.0000").rstrip()
        assert len(l) % 7 == 0, len(l)
        return [l[7*i:7*(i+1)] for i in range(len(l)//7)]

    depths = np.array(
        sum([_pre_format(line) for line in lines[i + 2:]], []),
        dtype=float)

    dispersion = np.arange(start, end + delta, delta)[:depths.size]
    intensity = 1.0 - depths
    
    # Parse the headers into metadata
    meta = {
        "raw": lines[:i + 2]
    }
    return (dispersion, intensity, meta)


def _parse_synth_summary(summary_out_path):
    """
    Parse the summary output from a MOOG synth operation.

    :param summary_out:
        The path of the summary output file produced by MOOG.
    """

    with open(summary_out_path, "r") as fp:
        summary = fp.readlines()

    # There could be multiple spectra in this output file, each separated by 
    # headers. Scan to find the headers, then extract the information.
    partition_indices = [i for i, line in enumerate(summary) \
        if line.lower().startswith("all abundances not listed below differ")] \
        + [None]

    spectra = []
    for i, start in enumerate(partition_indices[:-1]):
        end = partition_indices[i + 1]
        spectra.append(_parse_single_spectrum(summary[start:end]))

    return spectra


def _check_for_unphysical_intensity_jumps(transitions, dispersion, intensity,
    opacity_contribution):
    """
    Check for unphysical jumps in synthesized intensities that result from
    having too small a value of `opacity_contribution`.

    :param transitions:
        The transitions that contributed to this spectrum.

    :param dispersion:
        The dispersion values of the synthesized spectrum.

    :param intensity:
        The intensity values of the synthesized spectrum.

    :param opacity_contribution:
        The value of opacity contribution provided to MOOG.
    """


    # Identify pixels where the jumps could occur.
    pixels = np.hstack([
        dispersion.searchsorted(transitions["wavelength"]-opacity_contribution),
        dispersion.searchsorted(transitions["wavelength"]+opacity_contribution)
    ])

    l = lambda value: np.clip(value, 0, intensity.size - 1)
    derivatives = [np.diff(_) for _ in [intensity[l(px - 1): l(px + 1)] for px in pixels]]

    # Need some good reasoning for limiting what this should be -- perhaps from
    # H-alpha
    raise NotImplementedError
