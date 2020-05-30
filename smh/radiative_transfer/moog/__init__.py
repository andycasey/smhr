#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A standardized interface to the MOOG(SILENT) radiative transfer code. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["abundance_cog", "synthesize", "RTError"]

# See stackoverflow.com/questions/19913653/no-unicode-in-all-for-a-packages-init
# __all__ = [_.encode("ascii") for _ in __all__]

from .cog import abundance_cog
from .synthesis import synthesize
from .utils import RTError
