#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A photosphere class. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import logging
from textwrap import dedent

import astropy.io
import astropy.table

# Create logger.
logger = logging.getLogger(__name__)

class Photosphere(astropy.table.Table):
    """
    A model photosphere object.
    """
    pass


# MOOG writer and identifier.
def _moog_writer(photosphere, filename, **kwargs):
    """
    Writes an :class:`oracle.photospheres.photosphere` to file in a MOOG-friendly
    format.
    """

    def _get_xi():
        xi = photosphere.meta["stellar_parameters"].get("microturbulence", 0.0)
        if 0 >= xi:
            logger.warn("Invalid microturbulence value: {:.3f} km/s".format(xi))
        return xi

    if photosphere.meta["kind"] == "marcs":

        xi = _get_xi()

        output = dedent("""
            WEBMARCS
             ORACLE 1D MARCS (2011) TEFF/LOGG/[M/H]/XI {1:.0f}/{2:.3f}/{3:.3f}/{4:.3f}
            NTAU       {0:.0f}
            5000.0
            """.format(len(photosphere),
                photosphere.meta["stellar_parameters"]["effective_temperature"],
                photosphere.meta["stellar_parameters"]["surface_gravity"],
                photosphere.meta["stellar_parameters"]["metallicity"],
                xi)).lstrip()

        for i, line in enumerate(photosphere):
            output += " {0:>3.0f} {0:>3.0f} {1:10.3e} {0:>3.0f} {2:10.3e} "\
                "{3:10.3e} {4:10.3e}\n".format(i + 1, line["lgTau5"], line["T"],
                    line["Pe"], line["Pg"])

        output += "         {0:.3f}\n".format(xi)
        output += "NATOMS        0     {0:.3f}\n".format(
            photosphere.meta["stellar_parameters"]["metallicity"])
        output += "NMOL          0\n"


    elif photosphere.meta["kind"] == "castelli/kurucz":

        xi = _get_xi()
        
        output = dedent("""
            KURUCZ
             ORACLE 1D CASTELLI/KURUCZ (2004) TEFF/LOGG/[M/H]/[alpha/M]/XI {1:.0f}/{2:.3f}/{3:.3f}/{4:.3f}/{5:.3f}
            NTAU       {0:.0f}
            """.format(len(photosphere),
                photosphere.meta["stellar_parameters"]["effective_temperature"],
                photosphere.meta["stellar_parameters"]["surface_gravity"],
                photosphere.meta["stellar_parameters"]["metallicity"],
                photosphere.meta["stellar_parameters"]["alpha_enhancement"],
                xi)).lstrip()

        for line in photosphere:
            output += " {0:.8e} {1:10.3e}{2:10.3e}{3:10.3e}{4:10.3e}\n".format(
                line["RHOX"], line["T"], line["P"], line["XNE"], line["ABROSS"])

        output += "         {0:.3f}\n".format(xi)
        output += "NATOMS        0     {0:.3f}\n".format(
            photosphere.meta["stellar_parameters"]["metallicity"])
        output += "NMOL          0\n"

    else:
        raise ValueError("photosphere kind '{}' cannot be written to a MOOG-"\
            "compatible format".format(photosphere.meta["kind"]))

    with open(filename, "w") as fp:
        fp.write(output)

    return None


def _moog_identifier(*args, **kwargs):
    return isinstance(args[0], basestring) and args[0].lower().endswith(".moog")

# Register the MOOG writer and identifier
astropy.io.registry.register_writer("moog", Photosphere, _moog_writer)
astropy.io.registry.register_identifier("moog", Photosphere, _moog_identifier)
