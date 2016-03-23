#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" An object to manage SMH sessions. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["Session"]


from six import string_types

import specutils


class BaseSession(object):
    """
    An abstract class for a SMH session.
    """
    pass


class Session(BaseSession):
    """
    An object to manage a session in SMH.
    """

    def __init__(self, spectrum_paths, **kwargs):
        """
        Create a new session from input spectra.

        :param spectrum_paths:
            Filename of a single spectrum, or an iterable of spectrum filenames.
        """

        if isinstance(spectrum_paths, string_types):
            spectrum_paths = (spectrum_paths, )

        # Load the spectra and flatten all orders into a single list.
        input_spectra = list(map(specutils.Spectrum1D.read, spectrum_paths))
        input_spectra = sum([], input_spectra)

        # Sort orders from blue to red.

        # Store the path names internally for provenance.

        # Extract basic metadata information from the spectrum headers if
        # possible: RA, DEC, OBJECT

        # Load default session settings.

        raise NotImplementedError


    @classmethod
    def from_filename(cls, session_path, **kwargs):
        """
        Create a Session from a path saved to disk.
        """

        raise NotImplementedError