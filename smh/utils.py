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
import tempfile
from six import string_types

from collections import Counter, OrderedDict

from subprocess import getstatusoutput
from hashlib import sha1 as sha
from random import choice
from socket import gethostname, gethostbyname

# Third party imports
import numpy as np
import astropy.table
from scipy import stats, integrate, optimize

common_molecule_name2Z = {
    'Mg-H': 12,'H-Mg': 12,
    'C-C':  6,
    'C-N':  7, 'N-C':  7, #TODO
    'C-H':  6, 'H-C':  6,
    'O-H':  8, 'H-O':  8,
    'Fe-H': 26,'H-Fe': 26,
    'N-H':  7, 'H-N':  7,
    'Si-H': 14,'H-Si': 14,
    'Ti-O': 22,'O-Ti': 22,
    'V-O':  23,'O-V':  23,
    'Zr-O': 40,'O-Zr': 40
    }
common_molecule_name2species = {
    'Mg-H': 112,'H-Mg': 112,
    'C-C':  606,
    'C-N':  607,'N-C':  607,
    'C-H':  106,'H-C':  106,
    'O-H':  108,'H-O':  108,
    'Fe-H': 126,'H-Fe': 126,
    'N-H':  107,'H-N':  107,
    'Si-H': 114,'H-Si': 114,
    'Ti-O': 822,'O-Ti': 822,
    'V-O':  823,'O-V':  823,
    'Zr-O': 840,'O-Zr': 840
    }
common_molecule_species2elems = {
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

__all__ = ["element_to_species", "element_to_atomic_number", "species_to_element", "get_common_letters", \
    "elems_isotopes_ion_to_species", "species_to_elems_isotopes_ion", \
    "find_common_start", "extend_limits", "get_version", \
    "approximate_stellar_jacobian", "approximate_sun_hermes_jacobian",\
    "hashed_id"]

logger = logging.getLogger(__name__)


def mkdtemp(**kwargs):
    if not os.path.exists(os.environ["HOME"]+"/.smh"):
        logger.info("Making "+os.environ["HOME"]+"/.smh")
        os.mkdir(os.environ["HOME"]+"/.smh")
    if 'dir' not in kwargs:
        kwargs['dir'] = os.environ["HOME"]+"/.smh"
    return tempfile.mkdtemp(**kwargs)
def mkstemp(**kwargs):
    if not os.path.exists(os.environ["HOME"]+"/.smh"):
        logger.info("Making "+os.environ["HOME"]+"/.smh")
        os.mkdir(os.environ["HOME"]+"/.smh")
    if 'dir' not in kwargs:
        kwargs['dir'] = os.environ["HOME"]+"/.smh"
    return tempfile.mkstemp(**kwargs)

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

    :param columns: [optional]
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
            try: # fix for masked arrays
                finite = finite.filled(False)
            except:
                pass
            if not np.any(finite):
                #group_lines[x_column] = (np.nan, np.nan, np.nan, np.nan, 0)
                continue

            m, b, medy, stdy, stdm, N = fit_line(x, y, None)
            group_lines[x_column] = (m, b, medy, (stdy, stdm), N)
#            x, y, yerr = np.array(x[finite]), np.array(y[finite]), np.array(yerr[finite])
#
#            # Let's remove the covariance between m and b by making the mean of x = 0
#            xbar = np.mean(x)
#            x = x - xbar
#            # y = mx+b = m(x-xbar) + (b+m*xbar), so m is unchanged but b is shifted.
#
##            A = np.vstack((np.ones_like(x), x)).T
##            C = np.diag(yerr**2)
##            try:
##                cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
##                b, m = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))
##
##            except np.linalg.LinAlgError:
##                #group_lines[x_column] \
##                #    = (np.nan, np.nan, np.median(y), np.std(y), len(x))
##                None
##
##            else:
##                #group_lines[x_column] = (m, b, np.median(y), (np.std(y), np.sqrt(cov[1,1])), len(x))
##                group_lines[x_column] = (m, b+m*xbar, np.median(y), (np.std(y), np.sqrt(cov[1,1])), len(x))
#            m, b, r, p, m_stderr = stats.linregress(x, y)
#            group_lines[x_column] = (m, b-m*xbar, np.median(y), (np.std(y), m_stderr), len(x))


        identifier = transitions[group_by][start_index]
        if group_lines:
            lines[identifier] = group_lines

    return lines


def fit_line(x, y, yerr=None):
    if yerr is not None: raise NotImplementedError("Does not fit with error bars yet")
    finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum()==0:
        return np.nan, np.nan, np.nan, np.nan, np.nan, 0
    x, y = x[finite], y[finite]
    xbar = np.mean(x)
    x = x - xbar
    m, b_bar, r, p, m_stderr = stats.linregress(x, y)
    b = b_bar - m*xbar
    return m, b, np.median(y), np.std(y), m_stderr, len(x)

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

    line_list_hashes = line_list.compute_hashes()

    transition_hashes = {}
    for i, spectral_model in enumerate(spectral_models):
        for transition in spectral_model.transitions:
            transition_hash = line_list.hash(transition)
            transition_hashes.setdefault(transition_hash, [])
            transition_hashes[transition_hash].append(i)

    # Which of the transition_hashes appear more than once?
    conflicts = []
    for transition_hash, indices in transition_hashes.iteritems():
        if len(indices) < 2: continue

        # OK, what element is this transition?
        match = (line_list_hashes == transition_hash)
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
    return sha(salt.encode("utf-8")).hexdigest()
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


# E. Holmbeck added this jacobian to be used for solving only vt and feh.
# No idea if it's implemented properly...
def approximate_feh_jacobian(stellar_parameters, *args):
    """ Approximate the Jacobian of vt and feh and
    minimisation parameters, based on calculations from the Sun """

    logger.info("Updated approximation of the Jacobian")
	
    params_to_optimize = args[0]
	
    jacobian_params = [np.nan]*4
    next_param = 0
    # There must be a cleaner way to do this...
    for pi,p in enumerate(params_to_optimize):
        if p==True:
            jacobian_params[pi] = stellar_parameters[next_param]
            next_param += 1

    teff, vt, logg, feh = jacobian_params

    full_jacobian = np.array([
        [ 5.4393e-08*teff - 4.8623e-04, -7.2560e-02*vt + 1.2853e-01,  1.6258e-02*logg - 8.2654e-02,  1.0897e-02*feh - 2.3837e-02],
        [ 4.2613e-08*teff - 4.2039e-04, -4.3985e-01*vt + 8.0592e-02, -5.7948e-02*logg - 1.2402e-01, -1.1533e-01*feh - 9.2341e-02],
        [-3.2710e-08*teff + 2.8178e-04,  3.8185e-03*vt - 1.6601e-02, -1.2006e-02*logg - 3.5816e-03, -2.8592e-05*feh + 1.4257e-03],
        [-1.7822e-08*teff + 1.8250e-04,  3.5564e-02*vt - 1.1024e-01, -1.2114e-02*logg + 4.1779e-02, -1.8847e-02*feh - 1.0949e-01]
    ])

    full_jacobian = full_jacobian[params_to_optimize]
    return full_jacobian.T[params_to_optimize]


def approximate_sun_hermes_jacobian(stellar_parameters, *args):
    """
    Approximate the Jacobian of the stellar parameters and
    minimisation parameters, based on calculations using the Sun
    and the HERMES atomic line list, after equivalent widths
    were carefully inspected.
    """

#    logger.info("Updated approximation of the Jacobian")

    teff, vt, logg, feh = stellar_parameters[:4]

#    full_jacobian = np.array([
#        [ 4.4973e-08*teff - 4.2747e-04, -1.2404e-03*vt + 2.4748e-02,  1.6481e-02*logg - 5.1979e-02,  1.0470e-02*feh - 8.5645e-03],
#        [-9.3371e-08*teff + 6.9953e-04,  5.0115e-02*vt - 3.0106e-01, -6.0800e-02*logg + 6.7056e-02, -4.1281e-02*feh - 6.2085e-02],
#        [-2.1326e-08*teff + 1.9121e-04,  1.0508e-03*vt + 1.1099e-03, -6.1479e-03*logg - 1.7401e-02,  3.4172e-03*feh + 3.7851e-03],
#        [-9.4547e-09*teff + 1.1280e-04,  1.0033e-02*vt - 3.6439e-02, -9.5015e-03*logg + 3.2700e-02, -1.7947e-02*feh - 1.0383e-01]
#    ])

    # After culling abundance outliers,..
    full_jacobian = np.array([
        [ 4.5143e-08*teff - 4.3018e-04, -6.4264e-04*vt + 2.4581e-02,  1.7168e-02*logg - 5.3255e-02,  1.1205e-02*feh - 7.3342e-03],
        [-1.0055e-07*teff + 7.5583e-04,  5.0811e-02*vt - 3.1919e-01, -6.7963e-02*logg + 7.3189e-02, -4.1335e-02*feh - 6.0225e-02],
        [-1.9097e-08*teff + 1.8040e-04, -3.8736e-03*vt + 7.6987e-03, -6.4754e-03*logg - 2.0095e-02, -4.1837e-03*feh - 4.1084e-03],
        [-7.3958e-09*teff + 1.0175e-04,  6.5783e-03*vt - 3.6509e-02, -9.7692e-03*logg + 3.2322e-02, -1.7391e-02*feh - 1.0502e-01]
    ])
    return full_jacobian.T


def approximate_stellar_jacobian_2(stellar_parameters, *args):
    """ Approximate the Jacobian of the stellar parameters and
    minimisation parameters, based on calculations from the Sun """

    logger.info("Updated approximation of the Jacobian {}".format(stellar_parameters))

    teff, logg, vt, feh = stellar_parameters[:4]
    #if np.isnan(teff): teff = 5000.; logger.info("jacobian: teff=nan->5000")
    #if np.isnan(logg): logg = 2.0; logger.info("jacobian: logg=nan->2.0")
    #if np.isnan(vt): vt = 1.75; logger.info("jacobian: vt=nan->1.75")
    #if np.isnan(feh): feh = -2.0; logger.info("jacobian: feh=nan->-2.0")

    # This is the black magic.
    full_jacobian = np.array([
        [ 5.4393e-08*teff - 4.8623e-04,  1.6258e-02*logg - 8.2654e-02, -7.2560e-02*vt + 1.2853e-01,  1.0897e-02*feh - 2.3837e-02],
        [ 4.2613e-08*teff - 4.2039e-04, -5.7948e-02*logg - 1.2402e-01, -4.3985e-01*vt + 8.0592e-02, -1.1533e-01*feh - 9.2341e-02],
        [-3.2710e-08*teff + 2.8178e-04, -1.2006e-02*logg - 3.5816e-03,  3.8185e-03*vt - 1.6601e-02, -2.8592e-05*feh + 1.4257e-03],
        [-1.7822e-08*teff + 1.8250e-04, -1.2114e-02*logg + 4.1779e-02,  3.5564e-02*vt - 1.1024e-01, -1.8847e-02*feh - 1.0949e-01]
    ])
    return full_jacobian.T


def approximate_sun_hermes_jacobian_2(stellar_parameters, *args):
    """
    Approximate the Jacobian of the stellar parameters and
    minimisation parameters, based on calculations using the Sun
    and the HERMES atomic line list, after equivalent widths
    were carefully inspected.
    """

#    logger.info("Updated approximation of the Jacobian")

    teff, logg, vt, feh = stellar_parameters[:4]

#    full_jacobian = np.array([
#        [ 4.4973e-08*teff - 4.2747e-04, -1.2404e-03*vt + 2.4748e-02,  1.6481e-02*logg - 5.1979e-02,  1.0470e-02*feh - 8.5645e-03],
#        [-9.3371e-08*teff + 6.9953e-04,  5.0115e-02*vt - 3.0106e-01, -6.0800e-02*logg + 6.7056e-02, -4.1281e-02*feh - 6.2085e-02],
#        [-2.1326e-08*teff + 1.9121e-04,  1.0508e-03*vt + 1.1099e-03, -6.1479e-03*logg - 1.7401e-02,  3.4172e-03*feh + 3.7851e-03],
#        [-9.4547e-09*teff + 1.1280e-04,  1.0033e-02*vt - 3.6439e-02, -9.5015e-03*logg + 3.2700e-02, -1.7947e-02*feh - 1.0383e-01]
#    ])

    # After culling abundance outliers,..
    full_jacobian = np.array([
        [ 4.5143e-08*teff - 4.3018e-04,  1.7168e-02*logg - 5.3255e-02, -6.4264e-04*vt + 2.4581e-02,  1.1205e-02*feh - 7.3342e-03],
        [-1.0055e-07*teff + 7.5583e-04, -6.7963e-02*logg + 7.3189e-02,  5.0811e-02*vt - 3.1919e-01, -4.1335e-02*feh - 6.0225e-02],
        [-1.9097e-08*teff + 1.8040e-04, -6.4754e-03*logg - 2.0095e-02, -3.8736e-03*vt + 7.6987e-03, -4.1837e-03*feh - 4.1084e-03],
        [-7.3958e-09*teff + 1.0175e-04, -9.7692e-03*logg + 3.2322e-02,  6.5783e-03*vt - 3.6509e-02, -1.7391e-02*feh - 1.0502e-01]
    ])
    return full_jacobian.T


def _debytify(x):
    if isinstance(x, bytes):
        return x.decode("utf-8")
    return x
def _fix_bytes_dict(d):
    new_dict = {}
    for k,v in d.items():
        sk = _debytify(k)
        new_dict[sk] = v
    return new_dict

def element_to_species(element_repr):
    """ Converts a string representation of an element and its ionization state
    to a floating point """
    
    element_repr = _debytify(element_repr)
    if not isinstance(element_repr, string_types):
        raise TypeError("element must be represented by a string-type {} {}".format(element_repr, type(element_repr)))
        
    if element_repr.count(" ") > 0:
        element, ionization = element_repr.split()[:2]
    else:
        element, ionization = element_repr, "I"
    
    if element not in periodic_table:
        try:
            return common_molecule_name2species[element]
        except KeyError:
            # Don't know what this element is
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
    
    element_repr = _debytify(element_repr)
    if not isinstance(element_repr,  string_types):
        raise TypeError("element must be represented by a string-type {} {}".format(element_repr, type(element_repr)))
    
    element = element_repr.title().strip().split()[0]
    try:
        index = periodic_table.index(element)

    except IndexError:
        raise ValueError("unrecognized element '{}'".format(element_repr))
    except ValueError:
        try:
            return common_molecule_name2Z[element]
        except KeyError:
            raise ValueError("unrecognized element '{}'".format(element_repr))
        

    return 1 + index
    





def species_to_element(species):
    """ Converts a floating point representation of a species to a string
    representation of the element and its ionization state """
    
    if not isinstance(species, (float, int)):
        raise TypeError("species must be represented by a floating point-type {} {}".format(species, type(species)))
    
    if round(species,1) != species:
        # Then you have isotopes, but we will ignore that
        species = int(species*10)/10.

    if species + 1 >= len(periodic_table) or 1 > species:
        # Don"t know what this element is. It"s probably a molecule.
        try:
            elems = common_molecule_species2elems[species]
            return "-".join(elems)
        except KeyError:
            # No idea
            return str(species)
        
    atomic_number = int(species)
    element = periodic_table[int(species) - 1]
    ionization = int(round(10 * (species - int(species)) + 1))

    # The special cases
    if element in ("C", "H", "He"): return element
    return "%s %s" % (element, "I" * ionization)


def elems_isotopes_ion_to_species(elem1,elem2,isotope1,isotope2,ion,as_str=False):
    Z1 = int(element_to_species(elem1.strip()))
    if isotope1==0: isotope1=''
    else: isotope1 = str(isotope1).zfill(2)

    if elem2.strip()=='': # Atom
        mystr = "{}.{}{:03}".format(Z1,int(ion-1),int(isotope1))
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

    if as_str: return mystr
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
            elif element == 'H':
                elem1,_ion = 'H','I'
            elif element == 'He':
                elem1,_ion = 'He','I'
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



def struct2array(x):
    """ Convert numpy structured array of simple type to normal numpy array """
    Ncol = len(x.dtype)
    type = x.dtype[0].type
    assert np.all([x.dtype[i].type == type for i in range(Ncol)])
    return x.view(type).reshape((-1,Ncol))

    

def _make_rhomat(rho_Tg=0.0, rho_Tv=0.0, rho_TM=0.0, rho_gv=0.0, rho_gM=0.0, rho_vM=0.0):
    rhomat = np.array([[1.0, rho_Tg, rho_Tv, rho_TM],
                       [rho_Tg, 1.0, rho_gv, rho_gM],
                       [rho_Tv, rho_gv, 1.0, rho_vM],
                       [rho_TM, rho_gM, rho_vM, 1.0]])
    return rhomat
def process_session_uncertainties_lines(session, rhomat, minerr=0.001):
    """
    Using Sergey's estimator
    """
    from .spectral_models import ProfileFittingModel, SpectralSynthesisModel
    from .photospheres.abundances import asplund_2009 as solar_composition
    cols = ["index","wavelength","species","expot","loggf",
            "logeps","e_stat","eqw","e_eqw","fwhm",
            "e_Teff","e_logg","e_vt","e_MH","e_sys",
            "e_tot","weight"]
    data = OrderedDict(zip(cols, [[] for col in cols]))
    for i, model in enumerate(session.spectral_models):
        if not model.is_acceptable: continue
        if model.is_upper_limit: continue
        
        wavelength = model.wavelength
        species = np.ravel(model.species)[0]
        expot = model.expot
        loggf = model.loggf
        if np.isnan(expot) or np.isnan(loggf):
            print(i, species, model.expot, model.loggf)
        try:
            logeps = model.abundances[0]
            staterr = model.metadata["1_sigma_abundance_error"]
            if isinstance(model, SpectralSynthesisModel):
                (named_p_opt, cov, meta) = model.metadata["fitted_result"]
                if np.isfinite(cov[0,0]**0.5):
                    staterr = max(staterr, cov[0,0]**0.5)
                assert ~np.isnan(staterr)
            # apply minimum
            staterr = np.sqrt(staterr**2 + minerr**2)
            sperrdict = model.metadata["systematic_stellar_parameter_abundance_error"]
            e_Teff = sperrdict["effective_temperature"]
            e_logg = sperrdict["surface_gravity"]
            e_vt = sperrdict["microturbulence"]
            e_MH = sperrdict["metallicity"]
            e_all = np.array([e_Teff, e_logg, e_vt, e_MH])
            syserr_sq = e_all.T.dot(rhomat.dot(e_all))
            syserr = np.sqrt(syserr_sq)
            fwhm = model.fwhm
        except Exception as e:
            print("ERROR!!!")
            print(i, species, model.wavelength)
            print("Exception:",e)
            logeps, staterr, e_Teff, e_logg, e_vt, e_MH, syserr = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        if isinstance(model, ProfileFittingModel):
            eqw = model.equivalent_width or np.nan
            e_eqw = model.equivalent_width_uncertainty or np.nan
        else:
            eqw = -999
            e_eqw = -999
        #toterr = np.sqrt(staterr**2 + syserr**2)
        input_data = [i, wavelength, species, expot, loggf,
                      logeps, staterr, eqw, e_eqw, fwhm,
                      e_Teff, e_logg, e_vt, e_MH, syserr,
                      np.nan, np.nan]
        for col, x in zip(cols, input_data):
            data[col].append(x)
    tab = astropy.table.Table(data)

    # Calculate systematic error and effective weights for each species
    tab["e_sys"] = np.nan
    for species in np.unique(tab["species"]):
        ix = np.where(tab["species"]==species)[0]
        t = tab[ix]
        
        # Estimate systematic error s
        s = s_old = 0.
        #s_max = 2.
        s_max = 3. # Updated to a larger value because it was not always converging.
        delta = struct2array(t["e_Teff","e_logg","e_vt","e_MH"].as_array())
        ex = t["e_stat"]
        for i in range(100):
            sigma_tilde = np.diag(s**2 + ex**2) + (delta.dot(rhomat.dot(delta.T)))
            sigma_tilde_inv = np.linalg.inv(sigma_tilde)
            w = np.sum(sigma_tilde_inv, axis=1)
            
            xhat = np.sum(w*t["logeps"])/np.sum(w)
            dx = t["logeps"] - xhat
            def func(s):
                return np.sum(dx**2 / (ex**2 + s**2)**2) - np.sum(1/(ex**2 + s**2))
            if func(0) < func(s_max):
                s = 0
                break
            try:
                s = optimize.brentq(func, 0, s_max, xtol=.001)
            except ValueError as e:
                print("ERROR FOR SPECIES",species)
                print(e)
                print("s_max:",s_max)
                print("func(0)",func(0))
                print("func(s_max)",func(s_max))
                print("Figure out what you should do to s_max here:")
                import pdb; pdb.set_trace()
                raise
            
            if np.abs(s_old - s) < 0.01:
                break
            s_old = s
        else:
            print(species,"s did not converge!")
        print("Final in {} iter: {:.1f} {:.3f}".format(i+1, species, s))
        tab["e_sys"][ix] = s
        tab["e_tot"][ix] = np.sqrt(s**2 + ex**2)
        
        sigma_tilde = np.diag(tab["e_tot"][ix]**2) + (delta.dot(rhomat.dot(delta.T)))
        sigma_tilde_inv = np.linalg.inv(sigma_tilde)
        w = np.sum(sigma_tilde_inv, axis=1)
        wb = np.sum(sigma_tilde_inv, axis=0)
        assert np.allclose(w,wb,rtol=1e-6), "Problem in species {:.1f}, Nline={}, e_sys={:.2f}".format(species, len(t), s)
        tab["weight"][ix] = w
        
    for col in tab.colnames:
        if col in ["index", "wavelength", "species", "loggf", "star"]: continue
        tab[col].format = ".3f"
    
    return tab
def process_session_uncertainties_covariance(summary_tab, rhomat):
    ## Add in [X/Fe]
    # cov_XY = Cov(X,Y). Diagonal entries are Var(X). The matrix is symmetric.
    delta_XY = struct2array(np.array(summary_tab["e_Teff_w","e_logg_w","e_vt_w","e_MH_w"]))
    cov_XY = delta_XY.dot(rhomat.dot(delta_XY.T))
    assert np.all(np.abs(cov_XY - cov_XY.T) < 0.01**2), np.max(np.abs(np.abs(cov_XY - cov_XY.T)))
    # Add statistical errors to the diagonal
    #var_X = cov_XY[np.diag_indices_from(cov_XY)] + summary_tab["stderr_w"]**2
    var_X = summary_tab["e_XH"]**2 #cov_XY[np.diag_indices_from(cov_XY)] + 
    return var_X, cov_XY
def process_session_uncertainties_calc_xfe_errors(summary_tab, var_X, cov_XY):
    """
    Computes the following
    Var([X/Fe]) = Var(X) + Var(Fe) - 2 Cov(X, Fe)
    
    Does *not* compute covariances, but you can do that this way:
    Cov([X/Fe], [Fe/H]) = Cov(X,Fe) - Cov(Fe, Fe)
    """
    # [X/Fe] errors are the Fe1 and Fe2 parts of the covariance matrix
    try:
        ix1 = np.where(summary_tab["species"]==26.0)[0][0]
    except IndexError:
        print("No feh1: setting to nan")
        feh1 = np.nan
        exfe1 = np.nan
    else:
        feh1 = summary_tab["[X/H]"][ix1]
        var_fe1 = var_X[ix1]
        # Var(X/Fe1) = Var(X) + Var(Fe1) - 2*Cov(X,Fe1)
        exfe1 = np.sqrt(var_X + var_fe1 - 2*cov_XY[ix1,:])
    try:
        ix2 = np.where(summary_tab["species"]==26.1)[0][0]
    except IndexError:
        print("No feh2: setting to feh1")
        feh2 = feh1
        try:
            exfe2 = np.sqrt(var_X[ix1])
        except UnboundLocalError: # no ix1 either
            exfe2 = np.nan
    else:
        feh2 = summary_tab["[X/H]"][ix2]
        var_fe2 = var_X[ix2]
        # Var(X/Fe2) = Var(X) + Var(Fe2) - 2*Cov(X,Fe2)
        exfe2 = np.sqrt(var_X + var_fe2 - 2*cov_XY[ix2,:])
    
    return feh1, exfe1, feh2, exfe2
def process_session_uncertainties_abundancesummary(tab, rhomat):
    """
    Take a table of lines and turn them into standard abundance table
    """
    from .spectral_models import ProfileFittingModel, SpectralSynthesisModel
    from .photospheres.abundances import asplund_2009 as solar_composition
    unique_species = np.unique(tab["species"])
    cols = ["species","elem","N",
            "logeps","sigma","stderr",
            "logeps_w","sigma_w","stderr_w",
            "e_Teff","e_logg","e_vt","e_MH","e_sys",
            "e_Teff_w","e_logg_w","e_vt_w","e_MH_w","e_sys_w",
            "[X/H]","e_XH","s_X"]
    data = OrderedDict(zip(cols, [[] for col in cols]))
    for species in unique_species:
        ttab = tab[tab["species"]==species]
        elem = species_to_element(species)
        N = len(ttab)
        logeps = np.mean(ttab["logeps"])
        stdev = np.std(ttab["logeps"])
        stderr = stdev/np.sqrt(N)
        
        w = ttab["weight"]
        finite = np.isfinite(w)
        if finite.sum() != N:
            print("WARNING: species {:.1f} N={} != finite weights {}".format(species, N, finite.sum()))
        x = ttab["logeps"]
        logeps_w = np.sum(w*x)/np.sum(w)
        stdev_w = np.sqrt(np.sum(w*(x-logeps_w)**2)/np.sum(w))
        stderr_w = np.sqrt(1/np.sum(w))
        
        sperrs = []
        sperrs_w = []
        for spcol in ["Teff","logg","vt","MH"]:
            x_new = x + ttab["e_"+spcol]
            e_sp = np.mean(x_new) - logeps
            sperrs.append(e_sp)
            #e_sp_w = np.sum(w*x_new)/np.sum(w) - logeps_w
            e_sp_w = np.sum(w*ttab["e_"+spcol])/np.sum(w)
            sperrs_w.append(e_sp_w)
        sperrs = np.array(sperrs)
        sperrs_w = np.array(sperrs_w)
        sperrtot = np.sqrt(sperrs.T.dot(rhomat.dot(sperrs)))
        sperrtot_w = np.sqrt(sperrs_w.T.dot(rhomat.dot(sperrs_w)))
        
        XH = logeps_w - solar_composition(species)
        #e_XH = np.sqrt(stderr_w**2 + sperrtot_w**2)
        e_XH = stderr_w
        s_X = ttab["e_sys"][0]
        assert np.allclose(ttab["e_sys"], s_X), s_X
        input_data = [species, elem, N,
                      logeps, stdev, stderr,
                      logeps_w, stdev_w, stderr_w,
                      sperrs[0], sperrs[1], sperrs[2], sperrs[3], sperrtot,
                      sperrs_w[0], sperrs_w[1], sperrs_w[2], sperrs_w[3], sperrtot_w,
                      XH, e_XH, s_X
        ]
        assert len(cols) == len(input_data)
        for col, x in zip(cols, input_data):
            data[col].append(x)
    summary_tab = astropy.table.Table(data)
    
    ## Add in [X/Fe]
    var_X, cov_XY = process_session_uncertainties_covariance(summary_tab, rhomat)
    feh1, efe1, feh2, efe2 = process_session_uncertainties_calc_xfe_errors(summary_tab, var_X, cov_XY)
    
    if len(summary_tab["[X/H]"]) > 0:
        summary_tab["[X/Fe1]"] = summary_tab["[X/H]"] - feh1
        summary_tab["e_XFe1"] = efe1
        summary_tab["[X/Fe2]"] = summary_tab["[X/H]"] - feh2
        summary_tab["e_XFe2"] = efe2
    
        ixion = np.array([x - int(x) > .01 for x in summary_tab["species"]])
        summary_tab["[X/Fe]"] = summary_tab["[X/Fe1]"]
        summary_tab["e_XFe"] = summary_tab["e_XFe1"]
        summary_tab["[X/Fe]"][ixion] = summary_tab["[X/Fe2]"][ixion]
        summary_tab["e_XFe"][ixion] = summary_tab["e_XFe2"][ixion]
        for col in summary_tab.colnames:
            if col=="N" or col=="species" or col=="elem": continue
            summary_tab[col].format = ".3f"
    else:
        for col in ["[X/Fe]","[X/Fe1]","[X/Fe2]",
                    "e_XFe","e_XFe1","e_XFe2"]:
            summary_tab.add_column(astropy.table.Column(np.zeros(0),col))
            #summary_tab[col] = np.nan #.add_column(col)
    return summary_tab
def process_session_uncertainties_abundancesummary_totweight(tab, rhomat):
    """
    Take a table of lines and turn them into standard abundance table
    """
    from .spectral_models import ProfileFittingModel, SpectralSynthesisModel
    from .photospheres.abundances import asplund_2009 as solar_composition
    unique_species = np.unique(tab["species"])
    cols = ["species","elem","N",
            "logeps","sigma","stderr",
            "logeps_w","sigma_w","stderr_w",
            "e_Teff","e_logg","e_vt","e_MH","e_sys",
            "e_Teff_w","e_logg_w","e_vt_w","e_MH_w","e_sys_w",
            "[X/H]","e_XH","s_X"]
    data = OrderedDict(zip(cols, [[] for col in cols]))
    for species in unique_species:
        ttab = tab[tab["species"]==species]
        elem = species_to_element(species)
        N = len(ttab)
        logeps = np.mean(ttab["logeps"])
        stdev = np.std(ttab["logeps"])
        stderr = stdev/np.sqrt(N)
        
        ## NOTE: PART 1 OF 3 THAT IS CHANGING
        #w = ttab["weight"]
        w = ttab["e_tot"]**-2

        finite = np.isfinite(w)
        if finite.sum() != N:
            print("WARNING: species {:.1f} N={} != finite weights {}".format(species, N, finite.sum()))
        x = ttab["logeps"]
        logeps_w = np.sum(w*x)/np.sum(w)
        stdev_w = np.sqrt(np.sum(w*(x-logeps_w)**2)/np.sum(w))
        stderr_w = np.sqrt(1/np.sum(w))
        
        sperrs = []
        sperrs_w = []
        for spcol in ["Teff","logg","vt","MH"]:
            x_new = x + ttab["e_"+spcol]
            e_sp = np.mean(x_new) - logeps
            sperrs.append(e_sp)
            #e_sp_w = np.sum(w*x_new)/np.sum(w) - logeps_w
            e_sp_w = np.sum(w*ttab["e_"+spcol])/np.sum(w)
            sperrs_w.append(e_sp_w)
        sperrs = np.array(sperrs)
        sperrs_w = np.array(sperrs_w)
        sperrtot = np.sqrt(sperrs.T.dot(rhomat.dot(sperrs)))
        sperrtot_w = np.sqrt(sperrs_w.T.dot(rhomat.dot(sperrs_w)))
        
        XH = logeps_w - solar_composition(species)
        ## NOTE: PART 2 OF 3 THAT CHANGED
        e_XH = np.sqrt(stderr_w**2 + sperrtot**2)
        #e_XH = stderr_w
        s_X = ttab["e_sys"][0]
        assert np.allclose(ttab["e_sys"], s_X), s_X
        input_data = [species, elem, N,
                      logeps, stdev, stderr,
                      logeps_w, stdev_w, stderr_w,
                      sperrs[0], sperrs[1], sperrs[2], sperrs[3], sperrtot,
                      sperrs_w[0], sperrs_w[1], sperrs_w[2], sperrs_w[3], sperrtot_w,
                      XH, e_XH, s_X
        ]
        assert len(cols) == len(input_data)
        for col, x in zip(cols, input_data):
            data[col].append(x)
    summary_tab = astropy.table.Table(data)
    
    ## Add in [X/Fe]
    var_X, cov_XY = process_session_uncertainties_covariance(summary_tab, rhomat)
    ## NOTE: PART 3 OF 3 THAT IS ACTUALLY CHANGING
    #feh1, efe1, feh2, efe2 = process_session_uncertainties_calc_xfe_errors(summary_tab, var_X, cov_XY)
    ## Because we didn't do it self-consistently, these will be nan-ish.
    ## Here let's just do it by adding uncertainties in quadrature...
    try:
        ix1 = np.where(summary_tab["species"]==26.0)[0][0]
    except IndexError:
        print("No feh1: setting to nan")
        feh1 = np.nan
        efe1 = var_X * np.nan
    else:
        feh1 = summary_tab["[X/H]"][ix1]
        delta_XY = struct2array(np.array(summary_tab["e_Teff","e_logg","e_vt","e_MH"]))
        delta_XY = delta_XY - delta_XY[ix1,:]
        e2_sys = np.sqrt(np.sum(delta_XY**2, axis=1))
        efe1 = np.sqrt(var_X + e2_sys)
    try:
        ix2 = np.where(summary_tab["species"]==26.1)[0][0]
    except IndexError:
        print("No feh2: setting to nan")
        feh2 = np.nan
        efe2 = var_X * np.nan
    else:
        feh2 = summary_tab["[X/H]"][ix2]
        delta_XY = struct2array(np.array(summary_tab["e_Teff","e_logg","e_vt","e_MH"]))
        delta_XY = delta_XY - delta_XY[ix2,:]
        e2_sys = np.sum(delta_XY**2, axis=1)
        efe2 = np.sqrt(var_X + e2_sys)
    
    if len(summary_tab["[X/H]"]) > 0:
        summary_tab["[X/Fe1]"] = summary_tab["[X/H]"] - feh1
        summary_tab["e_XFe1"] = efe1
        summary_tab["[X/Fe2]"] = summary_tab["[X/H]"] - feh2
        summary_tab["e_XFe2"] = efe2
    
        ixion = np.array([x - int(x) > .01 for x in summary_tab["species"]])
        summary_tab["[X/Fe]"] = summary_tab["[X/Fe1]"]
        summary_tab["e_XFe"] = summary_tab["e_XFe1"]
        summary_tab["[X/Fe]"][ixion] = summary_tab["[X/Fe2]"][ixion]
        summary_tab["e_XFe"][ixion] = summary_tab["e_XFe2"][ixion]
        for col in summary_tab.colnames:
            if col=="N" or col=="species" or col=="elem": continue
            summary_tab[col].format = ".3f"
    else:
        for col in ["[X/Fe]","[X/Fe1]","[X/Fe2]",
                    "e_XFe","e_XFe1","e_XFe2"]:
            summary_tab.add_column(astropy.table.Column(np.zeros(0),col))
            #summary_tab[col] = np.nan #.add_column(col)
    return summary_tab
def process_session_uncertainties_limits(session, tab, summary_tab, rhomat):
    from .spectral_models import ProfileFittingModel, SpectralSynthesisModel
    from .photospheres.abundances import asplund_2009 as solar_composition
    ## Add in upper limits to line data
    cols = ["index","wavelength","species","expot","loggf",
            "logeps","e_stat","eqw","e_eqw","fwhm",
            "e_Teff","e_logg","e_vt","e_MH","e_sys",
            "e_tot","weight"]
    var_X, cov_XY = process_session_uncertainties_covariance(summary_tab, rhomat)
    feh1, efe1, feh2, efe2 = process_session_uncertainties_calc_xfe_errors(summary_tab, var_X, cov_XY)

    assert len(cols)==len(tab.colnames)
    data = OrderedDict(zip(cols, [[] for col in cols]))
    for i, model in enumerate(session.spectral_models):
        if not model.is_upper_limit: continue
        if not model.is_acceptable: continue
        
        wavelength = model.wavelength
        species = np.ravel(model.species)[0]
        expot = model.expot or np.nan
        loggf = model.loggf or np.nan
        try:
            logeps = model.abundances[0]
        except:
            logeps = np.nan

        input_data = [i, wavelength, species, expot, loggf,
                      logeps, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                      np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        for col, x in zip(cols, input_data):
            data[col].append(x)
    tab_ul = astropy.table.Table(data)
    tab_ul["logeps"].format = ".3f"
    tab = astropy.table.vstack([tab, tab_ul])
    
    ## Add in upper limits to summary table
    ul_species = np.unique(tab_ul["species"])
    cols = ["species","elem","N",
            "logeps","sigma","stderr",
            "logeps_w","sigma_w","stderr_w",
            "e_Teff","e_logg","e_vt","e_MH","e_sys",
            "e_Teff_w","e_logg_w","e_vt_w","e_MH_w","e_sys_w",
            "[X/H]","e_XH","s_X"] + ["[X/Fe1]","e_XFe1","[X/Fe2]","e_XFe2","[X/Fe]","e_XFe"]
    assert len(cols)==len(summary_tab.colnames)
    data = OrderedDict(zip(cols, [[] for col in cols]))
    for species in ul_species:
        if species in summary_tab["species"]: continue
        ttab_ul = tab_ul[tab_ul["species"]==species]
        elem = species_to_element(species)
        N = len(ttab_ul)
        limit_logeps = np.min(ttab_ul["logeps"])
        limit_XH = limit_logeps - solar_composition(species)
        limit_XFe1 = limit_XH - feh1
        limit_XFe2 = limit_XH - feh2
        limit_XFe = limit_XFe2 if (species - int(species) > .01) else limit_XFe1
        input_data = [species, elem, N,
                      limit_logeps, np.nan, np.nan,
                      limit_logeps, np.nan, np.nan,
                      np.nan, np.nan, np.nan, np.nan, np.nan,
                      np.nan, np.nan, np.nan, np.nan, np.nan,
                      limit_XH, np.nan, np.nan, limit_XFe1, np.nan, limit_XFe2, np.nan,
                      limit_XFe, np.nan
        ]
        for col, x in zip(cols, input_data):
            data[col].append(x)
    summary_tab_ul = astropy.table.Table(data)
    if len(summary_tab_ul) > 0:
        if len(summary_tab) > 0:
            summary_tab = astropy.table.vstack([summary_tab, summary_tab_ul])
        else:
            summary_tab = summary_tab_ul
    
    return tab, summary_tab
def process_session_uncertainties(session,
                                  rho_Tg=0.0, rho_Tv=0.0, rho_TM=0.0, rho_gv=0.0, rho_gM=0.0, rho_vM=0.0):
    """
    After you have run session.compute_all_abundance_uncertainties(),
    this pulls out a big array of line data
    and computes the final abundance table and errors
    
    By default assumes no correlations in stellar parameters. If you specify rho_XY
    it will include that correlated error.
    (X,Y) in [T, g, v, M]
    """
    ## Correlation matrix. This is multiplied by the errors to get the covariance matrix.
    # rho order = [T, g, v, M]
    rhomat = _make_rhomat(rho_Tg, rho_Tv, rho_TM, rho_gv, rho_gM, rho_vM)
    ## Make line measurement table (no upper limits yet)
    tab = process_session_uncertainties_lines(session, rhomat)
    ## Summarize measurements
    summary_tab = process_session_uncertainties_abundancesummary(tab, rhomat)
    ## Add upper limits
    tab, summary_tab = process_session_uncertainties_limits(session, tab, summary_tab, rhomat)
    return tab, summary_tab

def get_synth_eqw(model, window=1.0, wavelength=None,
                  get_spec=False):
    """
    Calculate the equivalent width associated with the synthetic line.
    This is done by synthesizing the line in absence of any other elements,
    then integrating the synthetic spectrum in a window around the central wavelength.
    
    The user can specify the size of the window (default +/-1A) 
    and the central wavelength (default None -> model.wavelength)
    """
    from .spectral_models import ProfileFittingModel, SpectralSynthesisModel
    assert isinstance(model, SpectralSynthesisModel)
    assert len(model.elements)==1, model.elements
    
    abundances = model.metadata["rt_abundances"].copy()
    for key in abundances:
        if key != model.elements[0]: abundances[key] = -9.0
    abundances[model.elements[0]] = model.metadata["fitted_result"][0].values()[0]
    print(abundances)
    
    synth_dispersion, intensities, meta = model.session.rt.synthesize(
        model.session.stellar_photosphere, model.transitions,
        abundances,
        isotopes=model.session.metadata["isotopes"], twd=model.session.twd)[0]
    if wavelength is None: wavelength = model.wavelength
    ii = (synth_dispersion > wavelength - window) & (synth_dispersion < wavelength + window)
    # integrate with the trapezoid rule, get milliangstroms
    eqw = 1000.*integrate.trapz(1.0-intensities[ii], synth_dispersion[ii])
    # integrate everything with the trapezoid rule, get milliangstroms
    eqw_all = 1000.*integrate.trapz(1.0-intensities, synth_dispersion)
    
    for key in abundances:
        abundances[key] = -9.0
    blank_dispersion, blank_flux, blank_meta = model.session.rt.synthesize(
        model.session.stellar_photosphere, model.transitions,
        abundances,
        isotopes=model.session.metadata["isotopes"], twd=model.session.twd)[0]
    blank_eqw = 1000.*integrate.trapz(1.0-blank_flux[ii], blank_dispersion[ii])
    # integrate everything with the trapezoid rule, get milliangstroms
    blank_eqw_all = 1000.*integrate.trapz(1.0-blank_flux, blank_dispersion)
    
    if get_spec:
        return eqw, eqw_all, blank_eqw, blank_eqw_all, synth_dispersion, intensities
    return eqw, eqw_all, blank_eqw, blank_eqw_all
