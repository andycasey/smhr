from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from six import iteritems

import numpy as np
import cPickle as pickle
import os

from .linelists import LineList
from .utils import element_to_species

__all__ = ['get_needed_isotopes','identify_isotopes','identify_needed_isotopes',
           'validate_isotopes','load_isotope_data']

_datadir = os.path.dirname(__file__)+'/data/isotopes'

def get_needed_isotopes(ll,isotopes):
    needed_isotopes = identify_needed_isotopes(ll)
    if len(needed_isotopes)==0: return {}
    bad_isotopes = []
    good_isotopes = {}
    for Z,isos in iteritems(needed_isotopes):
        if Z not in isotopes:
            bad_isotopes.append((Z,0))
        else:
            good_isos = {}
            for mass in isos:
                if mass not in isotopes[Z]:
                    bad_isotopes.append((Z,mass))
                else:
                    good_isos[mass] = isotopes[Z][mass]
            good_isotopes[Z] = good_isos
    if len(bad_isotopes) > 0:
        raise ValueError("Missing isotopes (Z,A): {}".format(bad_isotopes))
    return good_isotopes

def identify_isotopes(ll):
    """
    Find all isotopes in this linelist
    """
    assert isinstance(ll, LineList)
    cols = ['numelems','elem1','elem2','isotope1','isotope2']
    for col in cols:
        assert col in ll.colnames
    isotopes = {}
    for line in ll:
        nelem,e1,e2,i1,i2 = (line[col] for col in cols)
        #print(line['numelems','elem1','elem2','isotope1','isotope2'])
        #nelem,e1,e2,i1,i2 = line['numelems','elem1','elem2','isotope1','isotope2']
        Z1 = int(element_to_species(e1))
        if Z1 in isotopes:
            if i1 not in isotopes[Z1]:
                isotopes[Z1][i1] = 1.0
        else:
            isotopes[Z1] = {i1:1.0}
        if nelem==2:
            if Z2 in isotopes:
                if i2 not in isotopes[Z2]:
                    isotopes[Z2][i2] = 1.0
            else:
                isotopes[Z2] = {i2:1.0}
    return isotopes
def identify_needed_isotopes(ll):
    """
    Find elements with more than one isotope in this linelist
    """
    isotopes = identify_isotopes(ll)
    needed_isotopes = {}
    for Z,isos in iteritems(isotopes):
        if len(isos) > 1:
            needed_isotopes[Z] = isos
    return needed_isotopes

def validate_isotopes(isotopes,tol=1e-4):
    numbad = 0
    bad_isos = ""
    for Z in isotopes:
        total = 0.0
        for mass in isotopes[Z]:
            total += isotopes[Z][mass]
        if np.abs(total - 1.0) > tol: 
            bad_isos += "{} {};".format(Z,total-1)
            numbad += 1
    if numbad > 0: raise RuntimeError(bad_isos)

def load_isotope_data(whichdata):
    assert whichdata in ['rproc','sproc','sneden','asplund']
    datamap = {'rproc':'sneden08_rproc_isotopes.pkl',
               'sproc':'sneden08_sproc_isotopes.pkl',
               'sneden':'sneden08_all_isotopes.pkl',
               'asplund':'asplund09_isotopes.pkl'}
    with open(_datadir+'/'+datamap[whichdata],'r') as f:
        isotopes = pickle.load(f)
    validate_isotopes(isotopes)
    return isotopes
