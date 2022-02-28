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
from .utils import RTError
from smh.utils import element_to_species
from smh import linelists

logger = logging.getLogger(__name__)

# Load the MOOG defaults.
with resource_stream(__name__, "defaults.yaml") as fp:
    try:
        _moog_defaults = yaml.load(fp, yaml.FullLoader)
    except AttributeError:
        _moog_defaults = yaml.load(fp)

def abundance_cog(photosphere, transitions, full_output=False, verbose=False,
    twd=None, **kwargs):
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
    path = utils.twd_path(twd=twd,**kwargs)

    # Write out the photosphere.
    moog_in, model_in, lines_in \
        = path("batch.par"), path("model.in"), path("lines.in")
    photosphere.write(model_in, format="moog")


    # Write out the transitions.
    # Note that this must write out the EW too
    # Versions of MOOG < 2017 (e.g. 2014 and before) take log10 
    # of linelists with all positive loggf. Add a fake line to compensate.
    all_positive_loggf = np.all(transitions['loggf'] >= 0)
    if all_positive_loggf:
        logger.debug("Adding fake line with negative loggf!!!")
        # Add a fake line with negative loggf
        fakeline = linelists.LineList.create_basic_linelist([5006.126],[26.0],[2.833],[-3])
        fakeline[0]["equivalent_width"] = 50.
        transitions = linelists.table.vstack([fakeline, transitions])
    transitions.write(lines_in, format="moog")
    
    # Load the abfind driver template.
    with resource_stream(__name__, "abfind.in") as fp:
        template = fp.read().decode("utf-8")

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
    if code != 0:
        logger.error("MOOG returned the following standard output:")
        logger.error(out)
        logger.error("MOOG returned the following errors (code: {0:d}):".format(code))
        logger.error(err)
        logger.exception(RTError(err))
    else:
        logger.info("MOOG executed {0} successfully".format(moog_in))
        #logger.debug("Standard output:")
        #logger.debug(strip_control_characters(out))
        #logger.debug("Standard error:")
        #logger.debug(err.rstrip())

    # Parse the output.
    transitions_array, linear_fits = _parse_abfind_summary(kwds["summary_out"])

    if len(transitions_array)==0:
        logger.debug("Standard output:")
        logger.debug(strip_control_characters(out))
        logger.debug("Standard error:")
        logger.debug(err.rstrip())
        raise RTError("No measurements returned!")
    if len(transitions_array)!=len(transitions):
        logger.debug("Standard output:")
        logger.debug(strip_control_characters(out))
        logger.debug("Standard error:")
        logger.debug(err.rstrip())
        raise RTError("Num lines returned {} != {} Num lines input".format(len(transitions_array),len(transitions)))
    
    if all_positive_loggf:
        # Remove the fakeline
        transitions = transitions[1:]
        transitions_array = transitions_array[1:,:]

    # Match transitions. Check for anything missing.
    col_wl, col_species, col_ep, col_loggf, col_ew, col_logrw, col_abund, col_del_avg = range(8)
    moog_wl      = transitions_array[:,col_wl]
    moog_species = transitions_array[:,col_species]
    moog_abund   = transitions_array[:,col_abund]
    
    ii_orig = np.lexsort((transitions['species'],transitions['wavelength']))
    ii_moog = np.lexsort((moog_species,moog_wl))
    tol = .01
    matched_abund = np.zeros(len(transitions))*np.nan
    for ii1,ii2 in zip(ii_orig,ii_moog):
        assert transitions['species'][ii1]==moog_species[ii2]
        assert np.abs(transitions['wavelength'][ii1]-moog_wl[ii2]) < tol
        matched_abund[ii1] = moog_abund[ii2]
    
    # Return abundances w.r.t. the inputs.
    if full_output:
        raise NotImplementedError
        return (matched_abund, linear_fits)
    
    return matched_abund

    #raise NotImplementedError

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
            moog_slopes.setdefault(species, {})
            
            """
            if len(abundances) > 0:
                # Check to see if we already have abundances for this species
                exists = np.where(np.array(abundances)[:, 1] == species)
                
                #if len(exists[0]) > 0:
                #    raise a
                #    logger.debug("Detecting more than one iteration from MOOG")
                #    abundances = list(np.delete(abundances, exists, axis=0))
            """
            
        elif re.match("^\s{2,3}[0-9]", line):
            if species is None:
                raise IOError("Could not find the species!")
                
            line = list(map(float, line.split()))

            # July 2014 version of MOOG has an additional column w/ the species
            if len(line) == 7:
                line.insert(1, species)
                logger.debug("Detecting MOOG version < July 2014; "
                             "inserting species")
            #elif len(abundances) == 0:
            #    logger.debug("Detecting MOOG version >= July 2014; "
            #                 "not inserting species")
            abundances.append(line)
        
        elif 'corr. coeff.' in line:
            line = line.split()
            moog_slopes[species][name_map[line[0].replace('.', '').lower()]] \
                = list(map(float, [value.replace('D-', 'E-') for value in \
                    [line[4], line[7], line[11]]]))
            
    transitions_array = np.array(abundances, dtype=np.float).reshape((-1, 8))

    return (transitions_array, moog_slopes)
