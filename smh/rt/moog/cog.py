#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Functionality to derive abundances from a transition's curve-of-growth using the
MOOG(SILENT) radiative transfer code.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import numpy as np
import re
from pkg_resources import resource_stream

from . import utils
from smh.utils import element_to_species

logger = logging.getLogger(__name__)

# TODO: Put these into a defaults.yaml file?
_moog_defaults = {
    # atmosphere:
    # 0     do not output the atmosphere
    #*1     output standard information about the atmosphere
    # 2     output more details (continuous opacity factors, etc)
    "atmosphere": 0,

    # molecules:
    #*0     do not do molecular equilibrium
    # 1     do molecular equilibrium but do not output results
    # 2     do molecular equilibrium and output results
    "molecules": 1,

    # trudamp:
    # 0     no, stick with the standard damping formulae
    #*1     sure, why not? It's a black art anyway!
    "trudamp": 0,

    # lines:
    # 0     output nothing about the input lines
    #*1     output standard information about the input line list
    # 2     output line opacities at line centers
    # 3     output information about mean depth of line formation
    # 4     output the partition functions
    "lines": 0,

    # terminal:
    #*0     MOOG will query the user for the desired terminal type
    #       (not compatible with wrapper)
    # 7     Sunview (LickMONGO)
    # 11    Sun OpenWindows, or any X11 (LickMONGO)
    # 13    Graph-on 200 series, or any VT/Retrographics (LickMONGO)
    # X11   Sun OpenWindos, or any X11 (sm)
    # xterm xterm tektronix window (sm)
    # sunview SunView window(sm)
    # graphon graphon GO-250 (sm)
    "terminal": "x11",

    # flux/int:
    #*0     perform integrated flux calculations
    # 1     perform central intensity calculations
    "flux_int": 0,

    # damping:
    #*0     use the Unsold approximation, BUT: if a factor is read from the
    #       line list for an individual line, then if the factor is greater
    #       than 10^(-10), multiply the Unsold value by the factor,
    #       otherwise replace the Unsold value by the factor
    # 1     use the Unsold approximation multiplied by 6.3
    # 2     use the Unsold approximation multiplied by a factor recommended
    #       by the Blackwell group
    "damping": 2,

    # units:
    #*0     Angstroms
    # 1     microns
    "units": 0,

    # opacit:
    #*0     no fudge factor is to be applied
    # a     multiply the nominal continuous opacity by a factor: 10000a/T
    "opacit": 0

    # Others, requiring special treatment:
    # isotopes, abundances

    # Others, not specified here:
    # synlimits, fluxlimits, coglimits, blenlimits, iraf

    # Others, not allowed to be passed through this:
    # obspectrum, iraf, plot, freeform, strong, plotpars
    # TODO: revisit strong
}

def abundance_cog(photosphere, transitions, verbose=False, **kwargs):
    """
    Calculate atomic line abundances by interpolating their position from the
    curve-of-growth. This wraps the MOOG `abfind` driver.

    :param photosphere:
        A formatted photosphere.

    :param transitions:
        A list of atomic transitions with measured equivalent widths.

    :param verbose: [optional]
        Specify verbose flags to MOOG. This is primarily used for debugging.
    """

    # Create a temporary directory.
    path = utils.twd_path(**kwargs)

    # Write out the photosphere.
    moog_in, model_in, lines_in \
        = path("batch.par"), path("model.in"), path("lines.in")
    photosphere.write(model_in, format="moog")

    # Write out the transitions.
    transitions.write(lines_in, format="moog")
    
    # Load the abfind driver template.
    with resource_stream(__name__, "abfind.in") as fp:
        template = fp.read()

    # Not these are SMH defaults, not MOOG defaults.
    kwds = _moog_defaults.copy()
    if verbose:
        kwds.update({
            "atmosphere": 2,
            "molecules": 2,
            "lines": 3, # 4 is max verbosity, but MOOG falls over.
        })

    # Parse keyword arguments.
    kwds.update(kwargs)

    # Parse I/O files:
    kwds.update({
        "standard_out": path("abfind.std.out"),
        "summary_out": path("abfind.sum.out"),
        "model_in": model_in,
        "lines_in": lines_in,
    })
    contents = template.format(**kwds)

    # Write this to batch.par
    with open(moog_in, "w") as fp:
        fp.write(contents)

    # Execute MOOG in the TWD.
    code, out, err = utils.moogsilent(moog_in, **kwargs)

    # Returned normally?
    assert code == 0 # HACK # TODO

    # Parse the output.
    transitions_array, linear_fits = _parse_abfind_summary(kwds["summary_out"])

    # Match transitions. Check for anything missing.

    # Return abundances w.r.t. the inputs.

    raise NotImplementedError



def _parse_abfind_summary(summary_out_path):
    """
    Parse output from the summary output find from the MOOGSILENT `abfind`
    driver.

    :param summary_out_path:
        The path of the summary output file from the `abfind` driver.
    """

    with open(summary_out_path, "r") as fp:
        summary = fp.readlines()

    # Just load in all the abundance values, since everything else can be
    # calculated from them
    species = None
    abundances = []
    moog_slopes = {}
    
    # Map (characters are cheap)
    name_map = {
        'ep': 'excitation_potential',
        'rw': 'reduced_ew',
        'wav': 'wavelength'
    }
    
    for line in summary:
        if line.startswith('Abundance Results for Species'):
            
            _ = line.split('Species')[1].split('(')[0].strip()
            species = element_to_species(_)
            moog_slopes[species] = {}
            
            if len(abundances) > 0:
                # Check to see if we already have abundances for this species
                exists = np.where(np.array(abundances)[:, 1] == species)
                
                if len(exists[0]) > 0:
                    logger.debug("Detecting more than one iteration from MOOG")
                    abundances = list(delete(abundances, exists, axis=0))
            
        elif re.match("^\s{2,3}[0-9]", line):
            if species is None:
                raise IOError("Could not find the species!")
                
            line = map(float, line.split())

            # July 2014 version of MOOG has an additional column w/ the species
            if len(line) == 7:
                line.insert(1, species)

            elif len(abundances) == 0:
                logger.debug("Detecting MOOG version >= July 2014; "
                             "not inserting species")
            abundances.append(line)
        
        elif 'corr. coeff.' in line:
            line = line.split()
            moog_slopes[species][name_map[line[0].replace('.', '').lower()]] \
                = map(float, [value.replace('D-', 'E-') for value in \
                    [line[4], line[7], line[11]]])
            
    transitions_array = np.array(abundances, dtype=np.float)

    return (transitions_array, moog_slopes)
