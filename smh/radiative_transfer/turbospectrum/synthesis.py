#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Synthesis functionality using the Turbospectrum radiative transfer code. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging
import numpy as np
import yaml
from pkg_resources import resource_stream

from . import utils

logger = logging.getLogger(__name__)

# Load defaults.
with resource_stream(__name__, "defaults.yaml") as fp:
    _turbospectrum_defaults = yaml.load(fp)





def synthesize(photosphere, transitions, abundances=None, isotopes=None,
    verbose=False, twd=None, **kwargs):
    """
    Sythesize a stellar spectrum given the model photosphere and transitions
    provided. 

    Optional keyword arguments include:
    dispersion_min
    dispersion_max
    dispersion_delta
    opacity_contribution

    """


    kwds = _turbospectrum_defaults.copy()
    kwds.update(kwargs)

    # Set default regions to synthesize.
    kwds.setdefault("dispersion_min", 
        min(transitions["wavelength"]) - kwds["opacity_contribution"])
    kwds.setdefault("dispersion_max",
        max(transitions["wavelength"]) + kwds["opacity_contribution"] \
                                       + kwds["dispersion_delta"])


    # TODO put in abundances.
    if kwds["dispersion_max"] - kwds["dispersion_min"] > 1000.0:
        raise NotImplementedError("turbospectrum fails for large synth, "
                                  "and chunking not implemented yet")

    path = utils.twd_path(twd=twd, **kwargs)


    # Calculate opacities.
    babsma_lu_kwds = kwds.copy()
    babsma_lu_kwds["opacities_path"] = path("opacities.in")
    babsma_lu_kwds["photosphere_path"] = path("photosphere.in")

    # Write atmosphere file.
    photosphere.write(babsma_lu_kwds["photosphere_path"], format="turbospectrum")

    # Get babsma_lu template and format it.
    with resource_stream(__name__, "babsma_lu.in") as fp:
        babsma_lu_template = fp.read()

    # XI
    # num abundances
    # formatted abundances.

    # If no XI value, then estimate it.
    xi = photosphere.meta.get("microturbulence", None)
    if xi is None:
        xi = 1.0 # [km/s]
        logger.warn("No microturbulence given in photosphere. Assuming 1 km/s.")

    # TODO microturbulence may have been specified by kwargs....
    if "microturbulence" in babsma_lu_kwds:
        raise a

    babsma_lu_kwds["microturbulence"] = xi

    # Set some default values???

    babsma_lu_contents = babsma_lu_template.format(**babsma_lu_kwds)






    raise a
    # Get regions to synthesise.



    # Update abundances

    # Calculate opacities

    # Synthesise



