# coding: utf-8

""" Utility functions for Spectroscopy Made Hard """

__author__ = "Andy Casey <andy@astrowizici.st>"

# Standard library
import os
import logging
import platform
import string
import sys
import traceback

from collections import Counter

from commands import getstatusoutput
from hashlib import sha1 as sha
from random import choice
from socket import gethostname, gethostbyname

# Third party imports
import numpy as np

__all__ = ["element_to_species", "species_to_element", "get_common_letters", \
    "elems_isotopes_ion_to_species", "species_to_elems_isotopes_ion", \
    "find_common_start", "extend_limits", "get_version", \
    "approximate_stellar_jacobian", "approximate_sun_hermes_jacobian",\
    "hashed_id"]

logger = logging.getLogger(__name__)


def random_string(N=10):
    return ''.join(choice(string.ascii_uppercase + string.digits) for _ in range(N))


def equilibrium_state(transitions, columns=("expot", "rew"), group_by="species",
    ycolumn="abundance", yerr_column=None):
    """
    Perform linear fits to the abundances provided in the transitions table
    with respect to x-columns.

    :param transitions:
        A table of atomic transitions with measured equivalent widths and
        abundances.

    :param x: [optional]
        The names of the columns to make fits against.

    :param group_by: [optional]
        The name of the column in `transitions` to calculate states.
    """

    lines = {}
    transitions = transitions.group_by(group_by)
    for i, start_index in enumerate(transitions.groups.indices[:-1]):
        end_index = transitions.groups.indices[i + 1]

        # Do excitation potential first.
        group_lines = {}
        for x_column in columns:
            x = transitions[x_column][start_index:end_index]
            y = transitions["abundance"][start_index:end_index]

            if yerr_column is not None:
                try:
                    yerr = transitions[yerr_column][start_index:end_index]

                except KeyError:
                    logger.exception("Cannot find yerr column '{}':".format(
                        yerr_column))
                    
                    yerr = np.ones(len(y))
            
            else:
                yerr = np.ones(len(y))

            # Only use finite values.
            finite = np.isfinite(x * y * yerr)
            if not np.any(finite):
                #group_lines[x_column] = (np.nan, np.nan, np.nan, np.nan, 0)
                continue

            x, y, yerr = x[finite], y[finite], yerr[finite]

            A = np.vstack((np.ones_like(x), x)).T
            C = np.diag(yerr**2)
            try:
                cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
                b, m = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))

            except np.linalg.LinAlgError:
                #group_lines[x_column] \
                #    = (np.nan, np.nan, np.median(y), np.std(y), len(x))
                None

            else:
                group_lines[x_column] = (m, b, np.median(y), np.std(y), len(x))

        identifier = transitions[group_by][start_index]
        if group_lines:
            lines[identifier] = group_lines

    return lines


def spectral_model_conflicts(spectral_models, line_list):
    """
    Identify abundance conflicts in a list of spectral models.

    :param spectral_models:
        A list of spectral models to check for conflicts.

    :param line_list:
        A table of energy transitions.

    :returns:
        A list containing tuples of spectral model indices where there is a
        conflict about which spectral model to use for the determination of
        stellar parameters and/or composition.
    """

    transition_hashes = {}
    for i, spectral_model in enumerate(spectral_models):
        for transition_hash in spectral_model._transition_hashes:
            transition_hashes.setdefault(transition_hash, [])
            transition_hashes[transition_hash].append(i)

    # Which of the transition_hashes appear more than once?
    conflicts = []
    for transition_hash, indices in transition_hashes.iteritems():
        if len(indices) < 2: continue

        # OK, what element is this transition?
        match = (line_list["hash"] == transition_hash)
        element = line_list["element"][match][0].split()[0]

        # Of the spectral models that use this spectral hash, what are they
        # measuring?
        conflict_indices = []
        for index in indices:
            if element not in spectral_models[index].metadata["elements"]:
                # This transition is not being measured in this spectral model.
                continue

            else:
                # This spectral model is modeling this transition. 
                # Does it say this should be used for the determination of
                # stellar parameters or composition?
                if spectral_models[index].use_for_stellar_parameter_inference \
                or spectral_models[index].use_for_stellar_composition_inference:
                    conflict_indices.append(index)

        if len(conflict_indices) > 1:
            conflicts.append(conflict_indices)
    
    return conflicts







# List the periodic table here so that we can use it outside of a single
# function scope (e.g., 'element in utils.periodic_table')

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


def element_to_atomic_number(element_repr):
    """
    Converts a string representation of an element and its ionization state
    to a floating point.

    :param element_repr:
        A string representation of the element. Typical examples might be 'Fe',
        'Ti I', 'si'.
    """
    
    if not isinstance(element_repr, (unicode, str)):
        raise TypeError("element must be represented by a string-type")
    
    element = element_repr.title().strip().split()[0]
    try:
        index = periodic_table.index(element)

    except IndexError:
        raise ValueError("unrecognized element '{}'".format(element_repr))

    return 1 + index
    





def species_to_element(species):
    """ Converts a floating point representation of a species to a string
    representation of the element and its ionization state """
    
    if not isinstance(species, (float, int)):
        raise TypeError("species must be represented by a floating point-type")
    
    if round(species,1) != species:
        # Then you have isotopes, but we will ignore that
        species = int(species*10)/10.

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
            if len(list(set(elements_in_molecule))) == 1: return "{0}_{1}".format(elements_in_molecule[0], len(elements_in_molecule))

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


def elems_isotopes_ion_to_species(elem1,elem2,isotope1,isotope2,ion):
    Z1 = int(element_to_species(elem1.strip()))
    if isotope1==0: isotope1=''
    else: isotope1 = str(isotope1).zfill(2)

    if elem2.strip()=='': # Atom
        mystr = "{}.{}{}".format(Z1,int(ion-1),isotope1)
    else: # Molecule
        #assert ion==1,ion
        Z2 = int(element_to_species(elem2.strip()))

        # If one isotope is specified but the other isn't, use a default mass
        # These masses are taken from MOOG for Z=1 to 95
        amu = [1.008,4.003,6.941,9.012,10.81,12.01,14.01,16.00,19.00,20.18,
               22.99,24.31,26.98,28.08,30.97,32.06,35.45,39.95,39.10,40.08,
               44.96,47.90,50.94,52.00,54.94,55.85,58.93,58.71,63.55,65.37,
               69.72,72.59,74.92,78.96,79.90,83.80,85.47,87.62,88.91,91.22,
               92.91,95.94,98.91,101.1,102.9,106.4,107.9,112.4,114.8,118.7,
               121.8,127.6,126.9,131.3,132.9,137.3,138.9,140.1,140.9,144.2,
               145.0,150.4,152.0,157.3,158.9,162.5,164.9,167.3,168.9,173.0,
               175.0,178.5,181.0,183.9,186.2,190.2,192.2,195.1,197.0,200.6,
               204.4,207.2,209.0,210.0,210.0,222.0,223.0,226.0,227.0,232.0,
               231.0,238.0,237.0,244.0,243.0]
        amu = [int(round(x,0)) for x in amu]
        if isotope1 == '':
            if isotope2 == 0:
                isotope2 = ''
            else:
                isotope1 = str(amu[Z1-1]).zfill(2)
        else:
            if isotope2 == 0:
                isotope2 = str(amu[Z2-1]).zfill(2)
            else:
                isotope2 = str(isotope2).zfill(2)
        # Swap if needed
        if Z1 < Z2:
            mystr = "{}{:02}.{}{}{}".format(Z1,Z2,int(ion-1),isotope1,isotope2)
        else:
            mystr = "{}{:02}.{}{}{}".format(Z2,Z1,int(ion-1),isotope2,isotope1)

    return float(mystr)

def species_to_elems_isotopes_ion(species):
    element = species_to_element(species)
    if species >= 100:
        # Molecule
        Z1 = int(species/100)
        Z2 = int(species - Z1*100)
        elem1 = species_to_element(Z1).split()[0]
        elem2 = species_to_element(Z2).split()[0]
        # All molecules that we use are unionized
        ion = 1
        if species == round(species,1):
            # No isotope specified
            isotope1 = 0
            isotope2 = 0
        else: #Both isotopes need to be specified!
            isotope1 = int(species*1000) - int(species*10)*100
            isotope2 = int(species*100000) - int(species*1000)*100
            if isotope1 == 0 or isotope2 == 0: 
                raise ValueError("molecule species must have both isotopes specified: {} -> {} {}".format(species,isotope1,isotope2))
        # Swap if needed
    else:
        # Element
        try:
            elem1,_ion = element.split()
        except ValueError as e:
            if element == 'C':
                elem1,_ion = 'C','I'
            else:
                print(element)
                raise e
        ion = len(_ion)
        assert _ion == 'I'*ion, "{}; {}".format(_ion,ion)
        if species == round(species,1):
            isotope1 = 0
        elif species == round(species,4):
            isotope1 = int(species*10000) - int(species*10)*1000
        elif species == round(species,3):
            isotope1 = int(species*1000) - int(species*10)*100
        else:
            raise ValueError("problem determining isotope: {}".format(species))
        elem2 = ''
        isotope2 = 0
    return elem1,elem2,isotope1,isotope2,ion


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



    
