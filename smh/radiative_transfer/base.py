#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" An abstract base class to manage different radiative transfer codes. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)


import os
import logging
import subprocess
import tempfile
from shutil import rmtree


logger = logging.getLogger(__name__)


class RadiativeTransferBaseClass(object):
    pass



class RadiativeTransfer(RadiativeTransferBaseClass):

    _executables = ()

    def __init__(self, executable_path=None):
        """
        Initiate a radiative transfer code object.

        :param executable_path: [optional]
            The path of the executable.
        """

        self.executable_path = executable_path or self._find_executable()
        if not os.path.exists(self.executable_path):
            raise IOError(
                "executable {} does not exist".format(self.executable_path))
        return None


    def __repr__(self):
        return "<{0}.{1} object ({2} {3}) at {4}>".format(self.__module__, 
            type(self).__name__, self.name, self.version, hex(id(self)))

    def _find_executable(self):
        """ Find the executable for this radiative transfer code. """

        for executable in self._executables:
            try:
                executable_path = subprocess.check_output(
                    "which {}".format(executable), shell=True)

            except subprocess.CalledProcessError:
                continue

            else:
                return executable_path.strip()

        raise RuntimeError("cannot find executables: {}".format(
            ", ".join(self._executables)))

    @property
    def version(self):
        """ Return the version of the radiative transfer code. """
        return self._version


    @property
    def twd(self):
        """ The temporary working directory for the context manager. """
        return getattr(self, "_twd", None)


    def __enter__(self, remove_when_finished=True, **kwargs):
        """
        A context manager for this class.

        :param twd_path: [optional]
            The path to create a temporary working directory.

        :param remove_when_finished: [optional]
            Remove the temporary working directory (`twd_path`) when the context
            manager exits.
        """

        kwds = kwargs.copy()
        kwds.setdefault("dir", "/tmp/")
        kwds.setdefault("prefix", "smh-")

        self._twd = tempfile.mkdtemp(**kwds)
        self._remove_when_finished = True

        return self


    def __exit__(self, exc_type, exc_value, exc_traceback):
        """
        Exit method for the context manager.
        """

        if self._remove_when_finished and self.twd is not None:
            rmtree(self.twd, ignore_errors=True)
            delattr(self, "_twd")

        logger.info("{} {} {}".format(exc_type, exc_value, exc_traceback))
        
        return None


    def __call__(self, stdin=None, timeout=-1, **kwargs):
        """
        Open a subprocess call to execute the radiative transfer code.

        :param stdin: [optional]
            A string to pass as standard input when the code is executed.

        :param timeout: [optional]
            The number of seconds to wait before killing the subprocess.

        :returns:
            A three-length tuple containing the return code, standard output,
            and standard error.
        """

        class Alarm(Exception):
            pass

        def alarm_handler(signum, frame):
            raise Alarm


        kwds = kwargs.copy()
        kwds.setdefault("shell", True)
        kwds.setdefault("bufsize", 2056)
        kwds.setdefault("close_fds", True)
        #kwds.setdefault("cwd", self.twd or os.path.dirname(self.executable_path))
        kwds.setdefault("env", {
            "PATH": os.path.dirname(self.executable_path)
        })
        kwds.update({
            "stdin": subprocess.PIPE,
            "stdout": subprocess.PIPE,
            "stderr": subprocess.PIPE
        })

        p = subprocess.Popen([self.executable_path],  **kwds)

        if timeout != -1:
            signal.signal(signal.SIGALRM, alarm_handler)
            signal.alarm(timeout)

        try:
            stdout, stderr = p.communicate(input=stdin)

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


    def cog_abundance(self, *args, **kwargs):
        """
        Calculate abundances of atomic transitions from their position along the
        curve of growth.
        """
        raise NotImplementedError("this must be implemented by the base class")


    def synthesize(self, *args, **kwargs):
        """
        Synthesize a continuum-normalized stellar spectrum.
        """
        raise NotImplementedError("this must be implemented by the base class")



