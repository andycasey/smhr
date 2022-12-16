#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Spectroscopy Made Hard """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import logging
import os
from shutil import copyfile
from subprocess import check_output

# Software version.
__version__ = "0.2"

# Python 2/3 compatibility:
safe_check_output = lambda x, shell=True: check_output(x, shell=shell).decode(
    "ascii", "ignore")

try:
    git_hash = safe_check_output('git log -1 --date=short --format="%h"')
    unstaged_changes = safe_check_output('git status -s -uno')

except:
    __git_status__ = None

else:
    unstaged_changes = "*" if len(unstaged_changes) > 0 else ""
    __git_status__ = "".join([git_hash.strip(), unstaged_changes])
    del unstaged_changes, git_hash
    
# Set up logging.
logger = logging.getLogger("smh")
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))

#handler.setFormatter(logging.Formatter(
#    "%(asctime)s [%(levelname)-8s] (%(name)s/%(lineno)d): %(message)s"))

logger.addHandler(handler)


# Get the location for the Session defaults file.
from .session import Session
from . import (photospheres, radiative_transfer, spectral_models, robust_polyfit)
from .optimize_stellar_params import optimize_stellar_parameters, optimize_stellar_parameters_2, optimize_feh
from .linelists import LineList

# If there isn't a local copy of the default Session settings file, create one.
if not os.path.exists(Session._default_settings_path):
    copyfile(
        os.path.join(
            os.path.dirname(os.path.join(__file__)),
            "default_session.yaml"
        ),
        Session._default_settings_path)
