#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Atmosphere utilities. """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"


def estimate_microturbulence(effective_temperature, surface_gravity):
    """
    Estimate microtubulence from relations between effective temperature and
    surface gravity. For giants (logg < 3.5) the relationship employed is from
    Kirby et al. (2008, ) and for dwarfs (logg >= 3.5) the Reddy et al. (2003)
    relation is used.

    :param effective_temperature:
        The effective temperature of the star in Kelvin.

    :type effective_temperature:
        float

    :param surface_gravity:
        The surface gravity of the star.

    :type surface_gravity:
        float

    :returns:
        The estimated microturbulence (km/s) from the given stellar parameters.

    :rtype:
        float
    """

    if surface_gravity >= 3.5:
        return 1.28 + 3.3e-4 * (effective_temperature - 6000) \
            - 0.64 * (surface_gravity - 4.5)
    else:
        return 2.70 - 0.509 * surface_gravity


def element(atomic_number):
    """
    Return the element of a given atomic number.

    :param atomic_number:
        The atomic number for the element in question (e.g., 26).

    :type atomic_number:
        int-like

    :returns:
        The short-hand element for a given atomic number.

    :rtype:
        str
    """

    atomic_number = int(atomic_number)
    periodic_table = """H                                                  He
                        Li Be                               B  C  N  O  F  Ne
                        Na Mg                               Al Si P  S  Cl Ar
                        K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
                        Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
                        Cs Ba Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
                        Fr Ra Lr Rf Db Sg Bh Hs Mt Ds Rg Cn UUt"""
    
    lanthanoids    =   "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"
    actinoids      =   "Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No"
    
    periodic_table = periodic_table.replace(" Ba ", " Ba " + lanthanoids + " ") \
        .replace(" Ra ", " Ra " + actinoids + " ").split()
    del actinoids, lanthanoids
    return periodic_table[atomic_number - 1]