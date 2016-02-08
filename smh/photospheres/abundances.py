#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Solar abundances """

from __future__ import division, absolute_import, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

from numpy import array
from oracle.utils import element


def asplund_2009(elements):
    """
    Return the Asplund 2009 solar abundance for the given element(s).

    :param elements:
        Input elements provided as integer-like or string representations.

    :type element:
        int, str, or list/array-like of integer-like objects

    :returns:
        The abundances of the input elements given the Asplund (2009) solar
        composition.
    """

    asplund_2009 = {
        "Pr": 0.72, 
        "Ni": 6.22, 
        "Gd": 1.07, 
        "Pd": 1.57, 
        "Pt": 1.62, 
        "Ru": 1.75, 
        "S": 7.12, 
        "Na": 6.24, 
        "Nb": 1.46, 
        "Nd": 1.42, 
        "Mg": 7.6, 
        "Li": 1.05, 
        "Pb": 1.75, 
        "Re": 0.26, 
        "Tl": 0.9, 
        "Tm": 0.1, 
        "Rb": 2.52, 
        "Ti": 4.95, 
        "As": 2.3, 
        "Te": 2.18, 
        "Rh": 0.91, 
        "Ta": -0.12, 
        "Be": 1.38, 
        "Xe": 2.24, 
        "Ba": 2.18, 
        "Tb": 0.3, 
        "H": 12.0, 
        "Yb": 0.84, 
        "Bi": 0.65, 
        "W": 0.85, 
        "Ar": 6.4, 
        "Fe": 7.5, 
        "Br": 2.54, 
        "Dy": 1.1, 
        "Hf": 0.85, 
        "Mo": 1.88, 
        "He": 10.93, 
        "Cl": 5.5, 
        "C": 8.43, 
        "B": 2.7, 
        "F": 4.56, 
        "I": 1.55, 
        "Sr": 2.87, 
        "K": 5.03, 
        "Mn": 5.43, 
        "O": 8.69, 
        "Ne": 7.93, 
        "P": 5.41, 
        "Si": 7.51, 
        "Th": 0.02, 
        "U": -0.54, 
        "Sn": 2.04, 
        "Sm": 0.96, 
        "V": 3.93, 
        "Y": 2.21, 
        "Sb": 1.01, 
        "N": 7.83, 
        "Os": 1.4, 
        "Se": 3.34, 
        "Sc": 3.15, 
        "Hg": 1.17, 
        "Zn": 4.56, 
        "La": 1.1, 
        "Ag": 0.94, 
        "Kr": 3.25, 
        "Co": 4.99, 
        "Ca": 6.34, 
        "Ir": 1.38, 
        "Eu": 0.52, 
        "Al": 6.45, 
        "Ce": 1.58, 
        "Cd": 1.71, 
        "Ho": 0.48, 
        "Ge": 3.65, 
        "Lu": 0.1, 
        "Au": 0.92, 
        "Zr": 2.58, 
        "Ga": 3.04, 
        "In": 0.8, 
        "Cs": 1.08, 
        "Cr": 5.64, 
        "Cu": 4.19, 
        "Er": 0.92
    }

    def parse(x):

        if isinstance(x, (str, unicode)):
            return (asplund_2009[x], True)

        elif isinstance(x, (int, float)):
            el = element(x)
            return (asplund_2009[el], True)

        else:
            # Assume list-type
            return ([parse(el)[0] for el in x], False)

    abundances, is_scalar = parse(elements)

    if is_scalar:
        return abundances
        
    return array(abundances)