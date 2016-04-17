#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Utility functions related to MOOG. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import os
import signal
import subprocess
import tempfile
#from signal import alarm, signal, SIGALRM, SIGKILL


logger = logging.getLogger(__name__)

# Get the path of MOOGSILENT.
try:
    moogsilent_path = subprocess.check_output("which MOOGSILENT", shell=True)
except subprocess.CalledProcessError:
    logger.exception("Failed to find MOOGSILENT")
    raise IOError("cannot find MOOGSILENT")
    
else:
    moogsilent_path = moogsilent_path.strip()
    acceptable_moog_return_codes = (0, )

def twd_path(**kwargs):
    """
    Create a temporary working directory and return a function that will format
    basenames from that temporary working directory.
    """

    kwds = kwargs.copy()
    kwds.setdefault("dir", "/tmp/")
    kwds.setdefault("prefix", "smh-")
    twd = tempfile.mkdtemp(**kwds)
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

    logger.info("Executing MOOG input file: {0}".format(input_filename))

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

        stdout, stderr = p.communicate(input=pipe_input)
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
