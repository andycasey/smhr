# coding: utf-8

""" Utility functions for Spectroscopy Made Hard """

__author__ = "Andy Casey <andy@astrowizici.st>"

# Standard library
import os
import logging
import platform
import sys
import traceback

from commands import getstatusoutput
from hashlib import sha1 as sha
from socket import gethostname, gethostbyname

# Third party imports
import numpy as np

__all__ = ["element_to_species", "species_to_element", "get_common_letters", \
    "find_common_start", "extend_limits", "get_version", \
    "approximate_stellar_jacobian", "approximate_sun_hermes_jacobian",\
    "hashed_id"]

logger = logging.getLogger(__name__)


def hashed_id():
    try:
        salt = getstatusoutput("git config --get user.name")[1]
    except:
        import uuid
        salt = uuid.uuid3(uuid.NAMESPACE_DNS, "")
    return sha(salt).hexdigest()
hashed_id = hashed_id()


def approximate_stellar_jacobian(stellar_parameters, *args):
    """ Approximate the Jacobian of the stellar parameters and
    minimisation parameters, based on calculations from the Sun """

    logger.info("Updated approximation of the Jacobian")

    teff, vt, logg, feh = stellar_parameters[:4]

    # This is the black magic.
    full_jacobian = np.array([
        [ 5.4393e-08*teff - 4.8623e-04, -7.2560e-02*vt + 1.2853e-01,  1.6258e-02*logg - 8.2654e-02,  1.0897e-02*feh - 2.3837e-02],
        [ 4.2613e-08*teff - 4.2039e-04, -4.3985e-01*vt + 8.0592e-02, -5.7948e-02*logg - 1.2402e-01, -1.1533e-01*feh - 9.2341e-02],
        [-3.2710e-08*teff + 2.8178e-04,  3.8185e-03*vt - 1.6601e-02, -1.2006e-02*logg - 3.5816e-03, -2.8592e-05*feh + 1.4257e-03],
        [-1.7822e-08*teff + 1.8250e-04,  3.5564e-02*vt - 1.1024e-01, -1.2114e-02*logg + 4.1779e-02, -1.8847e-02*feh - 1.0949e-01]
    ])
    return full_jacobian.T


def approximate_sun_hermes_jacobian(stellar_parameters, *args):
    """
    Approximate the Jacobian of the stellar parameters and
    minimisation parameters, based on calculations using the Sun
    and the HERMES atomic line list, after equivalent widths
    were carefully inspected.
    """

    logger.info("Updated approximation of the Jacobian")

    teff, vt, logg, feh = stellar_parameters[:4]

    full_jacobian = np.array([
        [ 4.4973e-08*teff - 4.2747e-04, -1.2404e-03*vt + 2.4748e-02,  1.6481e-02*logg - 5.1979e-02,  1.0470e-02*feh - 8.5645e-03],
        [-9.3371e-08*teff + 6.9953e-04,  5.0115e-02*vt - 3.0106e-01, -6.0800e-02*logg + 6.7056e-02, -4.1281e-02*feh - 6.2085e-02],
        [-2.1326e-08*teff + 1.9121e-04,  1.0508e-03*vt + 1.1099e-03, -6.1479e-03*logg - 1.7401e-02,  3.4172e-03*feh + 3.7851e-03],
        [-9.4547e-09*teff + 1.1280e-04,  1.0033e-02*vt - 3.6439e-02, -9.5015e-03*logg + 3.2700e-02, -1.7947e-02*feh - 1.0383e-01]
    ])

    # After culling abundance outliers,..
    full_jacobian = np.array([
        [ 4.5143e-08*teff - 4.3018e-04, -6.4264e-04*vt + 2.4581e-02,  1.7168e-02*logg - 5.3255e-02,  1.1205e-02*feh - 7.3342e-03],
        [-1.0055e-07*teff + 7.5583e-04,  5.0811e-02*vt - 3.1919e-01, -6.7963e-02*logg + 7.3189e-02, -4.1335e-02*feh - 6.0225e-02],
        [-1.9097e-08*teff + 1.8040e-04, -3.8736e-03*vt + 7.6987e-03, -6.4754e-03*logg - 2.0095e-02, -4.1837e-03*feh - 4.1084e-03],
        [-7.3958e-09*teff + 1.0175e-04,  6.5783e-03*vt - 3.6509e-02, -9.7692e-03*logg + 3.2322e-02, -1.7391e-02*feh - 1.0502e-01]
    ])
    return full_jacobian.T


def element_to_species(element_repr):
    """ Converts a string representation of an element and its ionization state
    to a floating point """
    
    periodic_table = """H                                                  He
                        Li Be                               B  C  N  O  F  Ne
                        Na Mg                               Al Si P  S  Cl Ar
                        K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                        Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
                        Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                        Fr Ra Lr""" 
    
    lanthanoids    =   "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"
    actinoids      =   "Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"
    
    periodic_table = periodic_table.replace(" Ba ", " Ba " + lanthanoids + " ") \
        .replace(" Ra ", " Ra " + actinoids + " ").split()
    del actinoids, lanthanoids
    
    if not isinstance(element_repr, (unicode, str)):
        raise TypeError("element must be represented by a string-type")
        
    if element_repr.count(" ") > 0:
        element, ionization = element_repr.split()[:2]
    else:
        element, ionization = element_repr, "I"
    
    if element not in periodic_table:
        # Don"t know what this element is
        return float(element_repr)
    
    ionization = max([0, ionization.upper().count("I") - 1]) /10.
    transition = periodic_table.index(element) + 1 + ionization
    return transition


def species_to_element(species):
    """ Converts a floating point representation of a species to a string
    representation of the element and its ionization state """
    
    periodic_table = """H                                                  He
                        Li Be                               B  C  N  O  F  Ne
                        Na Mg                               Al Si P  S  Cl Ar
                        K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                        Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
                        Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                        Fr Ra Lr Rf"""
    
    lanthanoids    =   "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"
    actinoids      =   "Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"
    
    periodic_table = periodic_table.replace(" Ba ", " Ba " + lanthanoids + " ") \
        .replace(" Ra ", " Ra " + actinoids + " ").split()
    del actinoids, lanthanoids
    
    if not isinstance(species, (float, int)):
        raise TypeError("species must be represented by a floating point-type")
    
    if species + 1 >= len(periodic_table) or 1 > species:
        # Don"t know what this element is. It"s probably a molecule.
        common_molecules = {
            112: ["Mg", "H"],
            606: ["C", "C"],
            607: ["C", "N"],
            106: ["C", "H"],
            108: ["O", "H"],
            126: ["Fe", "H"],
            107: ["N", "H"],
            114: ["Si", "H"],
            822: ["Ti", "O"],
            823: ["V", "O"],
            840: ["Zr", "O"]
        }
        if species in common_molecules.keys():
            elements_in_molecule = common_molecules[species]
            if len(list(set(elements_in_molecule))): return "{0}_{1}".format(elements_in_molecule[0], len(elements_in_molecule))

            return "-".join(elements_in_molecule)

        else:
            # No idea
            return str(species)
        
    atomic_number = int(species)
    element = periodic_table[int(species) - 1]
    ionization = int(round(10 * (species - int(species)) + 1))

    # The special cases
    if element in ("C", "H", "He"): return element
    return "%s %s" % (element, "I" * ionization)


def get_common_letters(strlist):
    return "".join([x[0] for x in zip(*strlist) \
        if reduce(lambda a,b:(a == b) and a or None,x)])


def find_common_start(strlist):
    strlist = strlist[:]
    prev = None
    while True:
        common = get_common_letters(strlist)
        if common == prev:
            break
        strlist.append(common)
        prev = common

    return get_common_letters(strlist)


def extend_limits(values, fraction=0.10, tolerance=1e-2):
    """ Extend the values of a list by a fractional amount """

    values = np.array(values)
    finite_indices = np.isfinite(values)

    if np.sum(finite_indices) == 0:
        raise ValueError("no finite values provided")

    lower_limit, upper_limit = np.min(values[finite_indices]), np.max(values[finite_indices])
    ptp_value = np.ptp([lower_limit, upper_limit])

    new_limits = lower_limit - fraction * ptp_value, ptp_value * fraction + upper_limit

    if np.abs(new_limits[0] - new_limits[1]) < tolerance:
        if np.abs(new_limits[0]) < tolerance:
            # Arbitrary limits, since we"ve just been passed zeros
            offset = 1

        else:
            offset = np.abs(new_limits[0]) * fraction
            
        new_limits = new_limits[0] - offset, offset + new_limits[0]

    return np.array(new_limits)


def get_version():
    """ Retrieves the version of Spectroscopy Made Hard based on the
    git version """
  
    if getstatusoutput("which git")[0] == 0:
        git_commands = ("git rev-parse --abbrev-ref HEAD", "git log --pretty=format:'%h' -n 1")
        return "0.1dev:" + ":".join([getstatusoutput(command)[1] for command in git_commands])
    else:
        return "Unknown"
    
