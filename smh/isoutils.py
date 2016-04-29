from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from six import iteritems

import numpy as np
import cPickle as pickle
import os

from .linelists import LineList
from .utils import element_to_species

__all__ = ['IsotopeError',
           'get_needed_isotopes','identify_isotopes','identify_needed_isotopes',
           'validate_isotopes','load_isotope_data']

_datadir = os.path.dirname(__file__)+'/data/isotopes'

class IsotopeError(Exception):
    """Exception raised for missing isotopes"""
    def __init__(self,bad_isotopes):
        self.bad_isotopes = bad_isotopes
    def __str__(self):
        return repr(self.bad_isotopes)

def get_needed_isotopes(ll,new_isotopes):
    """
    Search the linelist for needed isotopes,
    tries to replace them with new_isotopes.
    Returns found_isotopes, missing_isotopes

    found_isotopes: needed isotopes that are found
    missing_isotopes: needed isotopes that are not found
    """
    needed_isotopes = identify_needed_isotopes(ll)
    if len(needed_isotopes)==0: return {}
    return find_isotopes_to_replace(needed_isotopes,new_isotopes)

def find_isotopes_to_replace(needed_isotopes, new_isotopes):
    """
    Given needed_isotopes, and new_isotopes which might the ratios you need
    Returns found_isotopes, missing_isotopes

    found_isotopes: needed isotopes that are found
    missing_isotopes: needed isotopes that are not found
    """
    missing_isotopes = {}
    found_isotopes = {}
    for elem,isos in iteritems(needed_isotopes):
        if elem not in new_isotopes:
            missing_isotopes[elem] = {0:np.nan}
        else:
            good_isos = {}
            for mass in isos:
                if mass not in new_isotopes[elem]:
                    missing_isotopes.append((elem,mass))
                else:
                    good_isos[mass] = new_isotopes[elem][mass]
            found_isotopes[elem] = good_isos
    return found_isotopes, missing_isotopes

def identify_isotopes(ll):
    """
    Find all isotopes in this linelist (even if there is only one isotope for an element)
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
        if e1 in isotopes:
            if i1 not in isotopes[e1]:
                isotopes[e1][i1] = 1.0
        else:
            isotopes[e1] = {i1:1.0}
        if nelem==2:
            if e2 in isotopes:
                if i2 not in isotopes[e2]:
                    isotopes[e2][i2] = 1.0
            else:
                isotopes[e2] = {i2:1.0}
    return isotopes
def identify_needed_isotopes(ll):
    """
    Find elements with more than one isotope in this linelist (i.e., the isotope ratio needs to be specified)
    """
    isotopes = identify_isotopes(ll)
    needed_isotopes = {}
    for elem,isos in iteritems(isotopes):
        if len(isos) > 1:
            needed_isotopes[elem] = isos
    return needed_isotopes

def validate_isotopes(isotopes,tol=1e-4):
    """
    Ensure isotope fractions for all isotopes sums to 1 (within tol=1e-4 by default)
    """
    numbad = 0
    bad_isos = {}
    for elem in isotopes:
        total = 0.0
        for mass in isotopes[elem]:
            total += isotopes[elem][mass]
        if np.abs(total - 1.0) > tol: 
            bad_isos[elem] = isotopes[elem]
            numbad += 1
    if numbad > 0: raise IsotopeError(bad_isos)

def load_isotope_data(whichdata):
    """
    Load isotope ratios: 'rproc', 'sproc', 'sneden', 'asplund'.
    'rproc': Sneden et al. 2008 r-process
    'sproc': Sneden et al. 2008 s-process
    'sneden': Sneden et al. 2008 total solar (only neutron-capture elements)
    'asplund': Asplund et al. 2009 solar (doesn't get the heaviest elements)
    """
    assert whichdata in ['rproc','sproc','sneden','asplund']
    datamap = {'rproc':'sneden08_rproc_isotopes.pkl',
               'sproc':'sneden08_sproc_isotopes.pkl',
               'sneden':'sneden08_all_isotopes.pkl',
               'asplund':'asplund09_isotopes.pkl'}
    with open(_datadir+'/'+datamap[whichdata],'r') as f:
        isotopes = pickle.load(f)
    validate_isotopes(isotopes)
    return isotopes
