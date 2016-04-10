#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Models for fitting spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["BaseSpectralModel"]

from smh.session import BaseSession


class BaseSpectralModel(object):
    """
    A base class for modelling spectra.
    """

    def __init__(self, transitions, session, *args, **kwargs):
        """
        Initialize a spectral model class.

        :param transitions:
            Row(s) from a line list of transitions to associate with this model.

        :param session:
            The parent session that this model will be associated with.
        """

        # TODO: validate the transitions

        if not isinstance(session, BaseSession):
            raise TypeError("session must be a sub-class from BaseSession")

        self._session = session
        self._transitions = transitions

        self.metadata = {}
        return None


    @property
    def transitions(self):
        """ Return the transitions associateed with this class. """
        return self._transitions


    @property
    def session(self):
        """ Return the parent session that this model is associated with. """
        return self._session


    @property
    def parameters(self):
        """
        Return the model parameters.
        This must be implemented by the sub-classes.
        """
        raise NotImplementedError(
            "parameters property must be implemented by the sub-classes")


    def __call__(self, dispersion, *args, **kwargs):
        """ The data-generating function. """
        raise NotImplementedError(
            "the data-generating function must be implemented by sub-classes")

