from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import json

from linelists import LineList
from utils import element_to_species

_datadir = './data/isotopes'

def get_needed_isotopes(ll,isotopes):
    needed_isotopes = identify_needed_isotopes(ll)
    if len(needed_isotopes)==0: return {}
    bad_isotopes = []
    for Z,nisos in needed_isotopes.iteritems():
        if Z not in isotopes:
            pass #TODO add bad isotopes
        else:
            for mass,frac in isos.iteritems():
                if mass not in isotopes[Z]:
                    pass #TODO add bad isotopes
                else:
                    needed_isotopes[Z][mass] = frac
    if len(bad_isotopes) > 0:
        pass #TODO raise error
    return needed_isotopes

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
    for Z,isos in isotopes.iteritems():
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

def load_isotopes(whichdata):
    assert whichdata in ['rproc','sproc','sneden','asplund']
    datamap = {'rproc':'sneden08_rproc_isotopes.json',
               'sproc':'sneden08_sproc_isotopes.json',
               'sneden':'sneden08_all_isotopes.json',
               'asplund':'asplund09_isotopes.json'}
    with open(_datadir+'/'+datamap[whichdata],'r') as f:
        isotopes = json.load(f)
    validate_isotopes(isotopes)
    return isotopes

def test_load_isotopes():
    for whichdata in ['rproc','sproc','sneden','asplund']:
        load_isotopes(whichdata)
def test_identify_isotopes():
    ll = LineList.read('lin4554new')
    isotopes = identify_isotopes(ll)
    for Z in isotopes:
        for mass in isotopes[Z]:
            print("{:2} {:3} {:.3f}".format(Z,mass,isotopes[Z][mass]))
            pass
    needed_isotopes = identify_needed_isotopes(ll)
    print()
    num_isotopes = 0
    for Z in needed_isotopes:
        for mass in needed_isotopes[Z]:
            print("{:2} {:3} {:.3f}".format(Z,mass,isotopes[Z][mass]))
            num_isotopes += 1
    assert num_isotopes == 5

if __name__=="__main__":
    test_identify_isotopes()
    test_load_isotopes()
    
