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
import yaml
from pkg_resources import resource_stream

from . import utils
from smh.utils import element_to_species

logger = logging.getLogger(__name__)

# Load the MOOG defaults.
with resource_stream(__name__, "defaults.yaml") as fp:
    _moog_defaults = yaml.load(fp)

class MOOGError(BaseException):
    pass

def abundance_cog(photosphere, transitions, verbose=False, **kwargs):
    """
    Calculate atomic line abundances by interpolating the measured 
    equivalent width from the curve-of-growth. 
    This wraps the MOOG `abfind` driver.

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
    # Note that this must write out the EW too
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

    # No isotopes here

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
    if code != 0:
        logger.error("MOOG returned the following standard output:")
        logger.error(out)
        logger.error("MOOG returned the following errors (code: {0:d}):".format(code))
        logger.error(err)
        logger.exception(MOOGError(err))
    else:
        logger.info("MOOG executed {0} successfully".format(moog_in))
        #logger.info("Standard output:")
        #logger.info(strip_control_characters(out))
        #logger.info("Standard error:")
        #logger.info(err.rstrip())

    # Parse the output.
    transitions_array, linear_fits = _parse_abfind_summary(kwds["summary_out"])

    # Match transitions. Check for anything missing.
    print(transitions_array)
    print(linear_fits)

    # Return abundances w.r.t. the inputs.

    raise NotImplementedError


def strip_control_characters(out):
    for x in np.unique(re.findall(r"\x1b\[K|\x1b\[\d+;1H",out)):
        out = out.replace(x,'')
    return out

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
                    abundances = list(np.delete(abundances, exists, axis=0))
            
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
