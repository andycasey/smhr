
#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Model for fitting synthetic spectra to spectral data. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["SpectralSynthesisModel"]

import logging

from .base import BaseSpectralModel

logger = logging.getLogger(__name__)


class SpectralSynthesisModel(BaseSpectralModel):
    """
    A class to fit spectral data by synthesis, given a photospheric structure.
    """
    pass
