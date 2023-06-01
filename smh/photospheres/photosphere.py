
#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" A photosphere class. """

from __future__ import division, absolute_import, print_function

import logging
from textwrap import dedent

import astropy.io
import astropy.table

# Create logger.
logger = logging.getLogger(__name__)

class Photosphere(astropy.table.Table):
    """ A model photosphere object. """
    pass


def _turbospectrum_writer(photosphere, filename, **kwargs):
    """
    Writes a :class:`photospheres.Photosphere` to file in a format that can be
    read by Turbospectrum.
    
    :param photosphere:
        The photosphere.

    :param filename:
        The filename to write the photosphere to.
    """

    if photosphere.meta["kind"] != "marcs":
        raise ValueError("only marcs photospheres supported with turbospectrum")

    output = ""
    radius = photosphere.meta["radius"]
    if radius > 0:
        output += \
            "spherical model\n"\
            "  1.0        Mass [Msun]\n"\
            "  {:.4e} Radius [cm] At Tau(Rosseland)=1.0\n".format(radius)

    else:
        output += \
            "plane-parallel model\n"\
            "  0.0        No mass for plane-parallel models\n"\
            "  1.0000E+00 1 cm radius for plane-parallel models\n"

    output += "  {:.0f} Number of depth points\n".format(len(photosphere))
    output += "Model structure\n"
    output += " k lgTauR  lgTau5    Depth     T        Pe          Pg         Prad       Pturb\n"

    lgTauR = -5
    for i, layer in enumerate(photosphere):
        lgTauR = -5 + i * 0.20
        output += "{:3.0f} {:5.2f} {:7.4f} {:10.3e} {:7.1f} {:10.3e} {:10.3e} {:10.3e} {:10.3e}\n"\
            .format(i + 1, lgTauR, layer["lgTau5"], layer["Depth"], layer["T"],
                layer["Pe"], layer["Pg"], layer["Prad"], 0)

    with open(filename, "w") as fp:
        fp.write(output)

    return None


def _moog_writer(photosphere, filename, **kwargs):
    """
    Writes an :class:`photospheres.photosphere` to file in a MOOG-friendly
    format.

    :param photosphere:
        The photosphere.

    :path filename:
        The filename to write the photosphere to.
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
             MARCS (2011) TEFF/LOGG/[M/H]/XI {1:.0f}/{2:.3f}/{3:.3f}/{4:.3f}
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

        output += "        {0:.3f}\n".format(xi)
        output += "NATOMS        0     {0:.3f}\n".format(
            photosphere.meta["stellar_parameters"]["metallicity"])
        output += "NMOL          0\n"


    elif photosphere.meta["kind"] == "castelli/kurucz":

        xi = _get_xi()
        
        output = dedent("""
            KURUCZ
             CASTELLI/KURUCZ (2004) {1:.0f}/{2:.3f}/{3:.3f}/{4:.3f}/{5:.3f}
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
            
        output += "        {0:.3f}\n".format(xi)
        output += "NATOMS        0     {0:.3f}\n".format(
            photosphere.meta["stellar_parameters"]["metallicity"])
        output += "NMOL          0\n"
        # MOOG11 fails to read if you don't add an extra line
        output += "\n"

    else:
        raise ValueError("photosphere kind '{}' cannot be written to a MOOG-"\
            "compatible format".format(photosphere.meta["kind"]))

    with open(filename, "w") as fp:
        fp.write(output)

    return None


# Register writers.
astropy.io.registry.register_writer("moog", Photosphere, _moog_writer)
astropy.io.registry.register_writer("turbospectrum", Photosphere, _turbospectrum_writer)

# Register identifiers.
astropy.io.registry.register_identifier("moog", Photosphere, 
    lambda *a, **k: "{0}".format(a[0]).lower().endswith(".moog"))
astropy.io.registry.register_identifier("turbospectrum", Photosphere,
    lambda *a, **k: "{0}".format(a[0]).lower().endswith(".turbospectrum"))
