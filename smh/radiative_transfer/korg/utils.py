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

class RTError(BaseException):
    pass

