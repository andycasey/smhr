from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from six import iteritems

import numpy as np
import pickle
import os

from .linelists import LineList
from .utils import element_to_species, element_to_atomic_number

__all__ = ['IsotopeError',
           'get_needed_isotopes','identify_isotopes','identify_needed_isotopes',
           'add_molecules',
           'validate_isotopes','load_isotope_data']

_datadir = os.path.dirname(__file__)+'/data/isotopes'

common_molecules = ['Mg-H','C-C','C-N','C-H','O-H','Fe-H','N-H','Si-H','Ti-O','V-O','Zr-O']
    

class IsotopeError(Exception):
    """Exception raised for missing isotopes"""
    def __init__(self,bad_isotopes):
        self.bad_isotopes = bad_isotopes
    def __str__(self):
        return repr(self.bad_isotopes)

def pretty_print_isotopes(isotopes):
    outstr = ""
    for elem in isotopes:
        outstr += str(elem)+': '
        for A,frac in iteritems(isotopes[elem]):
            outstr += '({},{:.3f}) '.format(A,frac)
        outstr += '\n'
    return outstr
def convert_isodict_to_array(isotopes,sort_by_Z=True):
    tab = []
    for elem in isotopes:
        for A in isotopes[elem]:
            # elem, A, frac
            tab.append([elem,A,isotopes[elem][A]])
    tab = np.array(tab)
    if sort_by_Z:
        def _sorter(x):
            elem,A = x
            A = int(A)
            # Works with molecules now
            Z = element_to_atomic_number(elem)
            if "-" in elem: # molecule
                A1 = int(A/100.)
                A2 = A-100*A1
                A = A1+A2
            return Z,A
        out = list(map(_sorter,zip(tab[:,0],tab[:,1])))
        out = np.array(out)
        ii = np.lexsort([out[:,1],out[:,0]])
        #Z = [int(element_to_species(elem)) for elem in tab[:,0]]
        #A = tab[:,1]
        #ii = np.lexsort([A,Z])
        tab = tab[ii]
    return tab
def convert_array_to_isodict(tab):
    isotopes = {}
    elems = np.unique(tab[:,0])
    for elem in elems:
        ttab = tab[tab[:,0]==elem]
        isotopes[elem] = dict(zip(ttab[:,1].astype(int),ttab[:,2].astype(float)))
    return isotopes

def get_needed_isotopes(ll,new_isotopes,include_molecules=False):
    """
    Search the linelist for needed isotopes,
    tries to replace them with new_isotopes.
    Returns found_isotopes, missing_isotopes

    found_isotopes: needed isotopes that are found
    missing_isotopes: needed isotopes that are not found
    """
    needed_isotopes = identify_needed_isotopes(ll,include_molecules)
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

def identify_isotopes(ll,include_molecules=False):
    """
    Find all isotopes in this linelist (even if there is only one isotope for an element).
    """
    assert isinstance(ll, LineList)
    cols = ['numelems','elem1','elem2','isotope1','isotope2']
    for col in cols:
        assert col in ll.colnames
    isotopes = {}
    for line in ll:
        nelem,e1,e2,i1,i2 = (line[col] for col in cols)
        if e1 in isotopes:
            if i1 not in isotopes[e1]:
                isotopes[e1][i1] = 1.0
        else:
            isotopes[e1] = {i1:1.0}
        if nelem==2:
            if e2 in isotopes:
                if i2 not in isotopes[e2]:
                    isotopes[e2][i2] = 1.0
                    if include_molecules:
                        if i1 < i2:
                            isotopes[e1+'-'+e2][i1*100+i2] = 1.0
                        else:
                            isotopes[e2+'-'+e1][i2*100+i1] = 1.0
            else:
                isotopes[e2] = {i2:1.0}
                if include_molecules:
                    if i1 < i2:
                        isotopes[e1+'-'+e2] = {i1*100+i2: 1.0}
                    else:
                        isotopes[e2+'-'+e1] = {i2*100+i1: 1.0}
    return isotopes
def identify_needed_isotopes(ll,include_molecules=False):
    """
    Find elements with more than one isotope in this linelist (i.e., the isotope ratio needs to be specified)
    """
    isotopes = identify_isotopes(ll,include_molecules)
    needed_isotopes = {}
    for elem,isos in iteritems(isotopes):
        if len(isos) > 1:
            needed_isotopes[elem] = isos
    return needed_isotopes

def add_molecules(isotopes,molecules=common_molecules):
    """
    Given a list of molecules to add to the isotopes dictionary, construct them from the
    element-by-element isotopes given.

    The key for isotope masses is int(A1*100+A2) (Z1 < Z2).
    This seems OK for now.
    """
    new_isotopes = isotopes.copy()
    for molecule in molecules:
        e1,e2 = molecule.split('-')
        Z1 = int(element_to_species(e1))
        Z2 = int(element_to_species(e2))
        if Z1 >= Z2: #swap e1,e2
            _e = e1; e1 = e2; e2 = _e
        elem = e1+'-'+e2

        if e1 in isotopes:
            iso1 = isotopes[e1]
        else:
            iso1 = {0:1.0}
        if e2 in isotopes:
            iso2 = isotopes[e2]
        else:
            iso2 = {0:1.0}

        iso = {}
        for A1 in iso1:
            for A2 in iso2:
                iso[A1*100+A2] = iso1[A1]*iso2[A2]
        new_isotopes[elem] = iso
    return new_isotopes

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
        try:
            e1,e2 = elem.split('-')
        except ValueError:
            pass
        else:
            Z1 = int(element_to_species(e1))
            Z2 = int(element_to_species(e2))
            assert Z1 <= Z2, "{} {}".format(Z1,Z2)
    if numbad > 0: raise IsotopeError(bad_isos)

def load_isotope_data(whichdata,include_molecules=False):
    """
    Load isotope ratios: 'rproc', 'sproc', 'sneden', 'asplund'.
    'rproc': Sneden et al. 2008 r-process
    'sproc': Sneden et al. 2008 s-process
    'sneden': Sneden et al. 2008 total solar (only neutron-capture elements)
    'asplund': Asplund et al. 2009 solar (doesn't get the heaviest elements)
    """
    assert whichdata in ['rproc','sproc','sneden','asplund']
    datamap = {'rproc':'sneden08_rproc_isotopes_unicode.pkl',
               'sproc':'sneden08_sproc_isotopes_unicode.pkl',
               'sneden':'sneden08_all_isotopes_unicode.pkl',
               'asplund':'asplund09_isotopes_unicode.pkl'}
    with open(_datadir+'/'+datamap[whichdata],'rb') as f:
        isotopes = pickle.load(f, encoding="latin1")
    if include_molecules:
        isotopes = add_molecules(isotopes)
    validate_isotopes(isotopes)
    return isotopes
