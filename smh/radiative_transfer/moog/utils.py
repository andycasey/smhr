#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Utility functions related to MOOG. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import numpy as np
import os
import signal
import subprocess
#import tempfile

from smh.photospheres.abundances import asplund_2009 as solar_composition
from smh.utils import elems_isotopes_ion_to_species, element_to_atomic_number
from smh.utils import mkdtemp
from six import iteritems, string_types

logger = logging.getLogger(__name__)

# Get the path of MOOGSILENT/moogsilent.
for executable in ("MOOGSILENT", "moogsilent"):
    try:
        moogsilent_path = subprocess.check_output(
            "which {}".format(executable), shell=True)
    except subprocess.CalledProcessError:
        continue
    else:
        moogsilent_path = moogsilent_path.strip()
        acceptable_moog_return_codes = (0, )

try:
    moogsilent_path
except NameError:    
    logger.exception("Failed to find moogsilent executable")
    raise IOError("cannot find MOOGSILENT")

class RTError(BaseException):
    pass

def twd_path(twd=None,**kwargs):
    """
    Create a temporary working directory and return a function that will format
    basenames from that temporary working directory.
    """

    if twd is None:
        kwds = {}
        kwds["dir"] = kwargs.get("dir", "/tmp/")
        kwds["prefix"] = kwargs.get("prefix", "smh-")
        kwds["suffix"] = kwargs.get("suffix", "")
        #kwds = kwargs.copy()
        #kwds.setdefault("dir", "/tmp/")
        #kwds.setdefault("prefix", "smh-")
        twd = mkdtemp(**kwds)
    if len(twd) > 30:
        logger.warn(
            "Temporary working directory should be as short as possible to "\
            "prevent MOOG(SILENT) from falling over. Current length ({0}): "\
            "{1}".format(len(twd), twd))
    return lambda filename: os.path.join(twd, filename)


class MOOGError(BaseException):
    pass



def moogsilent(input_filename, cwd=None, timeout=30, shell=False, env=None,
    **kwargs):
    """ 
    Execute a MOOGSILENT-compatible input file with a timeout after which it
    will be forcibly killed.

    :param input_filename:
        The full path of the MOOGSILENT-compatible input file.

    :param cwd: [optional]
        The current working directory to specify. If this is not specified, the
        directory of MOOGSILENT will be used.

    :param timeout: [optional]
        The number of seconds to wait before killing the process.

    :param shell: [optional]
        Whether to execute MOOGSILENT in a shell.

    :param env: [optional]
        A dictionary of environment variables to supply.
    """

    logger.debug("Executing MOOG input file: {}".format(input_filename))

    class Alarm(Exception):
        pass

    def alarm_handler(signum, frame):
        raise Alarm

    if cwd is None:
        cwd = os.path.dirname(input_filename)
    
    if env is None and len(os.path.dirname(moogsilent_path)) > 0:
        env = {"PATH": os.path.dirname(moogsilent_path)}

    p = subprocess.Popen([os.path.basename(moogsilent_path)], 
        shell=shell, bufsize=2056, cwd=cwd, 
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, env=env, close_fds=True)

    if timeout != -1:
        signal.signal(signal.SIGALRM, alarm_handler)
        signal.alarm(timeout)

    try:
        # Stromlo clusters may need a "\n" prefixed to the input for p.communicate
        pipe_input = "\n" if -6 in acceptable_moog_return_codes else ""
        pipe_input += os.path.basename(input_filename) + "\n"*100

        stdout, stderr = p.communicate(input=pipe_input.encode())

        # Parse the version of MOOG.
        #index = stdout.find("VERSION")
        #version = (" ".join(stdout[index:].split(" ")[1:3])).strip("()")
        #logger.debug("MOOG at {} is version {}".format(moogsilent_path, version))

        if timeout != -1:
            signal.alarm(0)

    except Alarm:

        # process might have died before getting to this line
        # so wrap to avoid OSError: no such process
        try:
            os.kill(p.pid, signal.SIGKILL)
        except OSError:
            pass
        return (-9, '', '')

    return (p.returncode, stdout, stderr)


def _format_abundances(elemental_abundances=None, subtract_solar=False,
    subtract_metallicity=0):

    if elemental_abundances is None:
        return ("0 1", 1)

    # First-pass to check the abundances and get the number of syntheses.
    elemental_abundances = elemental_abundances.copy()
    
    # Make sure that the abundances are specified as (atomic_number: abundance)
    for key in list(elemental_abundances.keys()):
        if isinstance(key, string_types):
            # It's an element. Convert to atomic number.
            atomic_number = element_to_atomic_number(key)
            if atomic_number in elemental_abundances.keys():
                raise ValueError(
                    "element {} has abundances specified by element and species"\
                    .format(key))
            elemental_abundances[atomic_number] = elemental_abundances.pop(key)

    sorted_atomic_numbers, max_synth = sorted(elemental_abundances.keys()), None
    for atomic_number in sorted_atomic_numbers:
        abundance = elemental_abundances[atomic_number]
        try:
            abundance[0]
        except (IndexError, TypeError):
            abundance = [abundance]
        abundance = np.array(abundance).flatten().copy()

        # Subtract any compositions.
        abundance -= subtract_metallicity
        if subtract_solar:
            abundance -= solar_composition(atomic_number)

        elemental_abundances[atomic_number] = abundance
        if max_synth is None or len(abundance) > max_synth:
            max_synth = len(abundance)

        if len(abundance) > 1 and len(abundance) != max_synth:
            raise ValueError("abundance entries must be fully described "
                             "or have one abundance per element (Z = {})"\
                             .format(atomic_number))

    assert 5 >= max_synth

    str_format = ["{0} {1}".format(len(sorted_atomic_numbers), max_synth)]
    for atomic_number in sorted_atomic_numbers:
        abundance = elemental_abundances[atomic_number]

        if len(abundance) == 1 and max_synth > 1:
            abundance = list(abundance) * max_synth

        _ = "      {0:3.0f} ".format(atomic_number)
        str_format.append(_ + " ".join(
            ["{0:8.3f}".format(a) for a in np.array(abundance).flatten()]))

    return ("\n".join(str_format), max_synth)




def _format_isotopes(isotopes=None, ionisation_states=(0, 1), num_synth=1):
    """
    Format isotopic information for a MOOG input file.

    :param isotopes:
        A dictionary containing isotopic information. The keys of the dictionary
        should specify elements/molecules and the value should have another
        dictionary that contains the mass number as keys, and the isotope 
        fraction as a value. The sum of all values in a sub-dictionary should 
        total 1.
    """

    if isotopes is None:
        return "0 {0:.0f}".format(num_synth) # No isotopes, N synthesis.

    if not isinstance(isotopes, dict):
        raise TypeError("isotopes must be provided as a dictionary")

    # TODO
    # For now, just output both ionization states in the isotopes
    # The easiest way to avoid this is to pass in the species as well
    N = 0
    outstr = []

    #fmt = "  {0:} {1:.3f} {1:.3f} {1:.3f}"
    fmt = "  "+" ".join(["{0:}"]+["{1:.3f}" for x in range(num_synth)])
    for elem in isotopes:
        for A,frac in iteritems(isotopes[elem]):
            if frac==0:
                invfrac = 99999.9
            else:
                invfrac = 1./frac

            try:
                e1,e2 = elem.split('-')
            except ValueError: #Single element
                e1 = elem
                A1 = A
                e2 = ''
                A2 = 0
            else: #Molecule
                A1 = int(A/100.)
                A2 = A-A1*100
            ## TODO temporary kludge is to just put both ionization states in for every isotope
            for ion in ionisation_states:
                species = elems_isotopes_ion_to_species(e1,e2,A1,A2,ion+1,as_str=True)
                outstr.append(fmt.format(species,invfrac))
                N += 1
    outstr.insert(0,"{} {:.0f}".format(N,num_synth))
    return "\n".join(outstr)
