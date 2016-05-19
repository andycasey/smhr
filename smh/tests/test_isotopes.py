from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from smh.isoutils import *
from smh.linelists import LineList
from smh.radiative_transfer import moog

import os
datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'

def test_get_needed_isotopes():
    ll = LineList.read(datadir+'/linelists/lin4554new')
    rproc = load_isotope_data('rproc')
    isotopes, missing_isotopes = get_needed_isotopes(ll,rproc)
    assert len(missing_isotopes) == 0, missing_isotopes
    validate_isotopes(isotopes)
    for elem in isotopes:
        for mass in isotopes[elem]:
            assert isotopes[elem][mass] == rproc[elem][mass]
def test_get_needed_isotopes2():
    ll = LineList.read(datadir+'/linelists/masseron_linch.txt')
    asplund = load_isotope_data('asplund')
    asplund = add_molecules(asplund)
    isotopes, missing_isotopes = get_needed_isotopes(ll,asplund,include_molecules=True)
    assert len(missing_isotopes) == 0, missing_isotopes
    validate_isotopes(isotopes)
    for elem in isotopes:
        for mass in isotopes[elem]:
            assert isotopes[elem][mass] == asplund[elem][mass]

def test_load_isotope_data():
    for whichdata in ['rproc','sproc','sneden','asplund']:
        load_isotope_data(whichdata)
        load_isotope_data(whichdata,include_molecules=True)
def test_identify_isotopes():
    ll = LineList.read(datadir+'/linelists/lin4554new')
    isotopes = identify_isotopes(ll)
    for elem in isotopes:
        for mass in isotopes[elem]:
            if elem != 'Ba': 
                assert mass == 0
            assert isotopes[elem][mass] == 1
    needed_isotopes = identify_needed_isotopes(ll)
    num_isotopes = 0
    for elem in needed_isotopes:
        for mass in needed_isotopes[elem]:
            num_isotopes += 1
    assert num_isotopes == 5
def test_identify_isotopes2():
    ll = LineList.read(datadir+'/linelists/masseron_linch.txt')
    isotopes = identify_isotopes(ll,include_molecules=True)
    for elem in isotopes:
        for mass in isotopes[elem]:
            #print(elem,mass,isotopes[elem][mass])
            if elem not in ['C','H-C','H']: 
                assert mass == 0
            elif elem == 'H':
                assert mass == 1
            elif elem == 'C':
                assert mass in [12,13]
            elif elem == 'H-C':
                assert mass in [112,113]
            assert isotopes[elem][mass] == 1
    needed_isotopes = identify_needed_isotopes(ll,include_molecules=True)
    num_isotopes = 0
    for elem in needed_isotopes:
        for mass in needed_isotopes[elem]:
            #print(elem,mass,isotopes[elem][mass])
            num_isotopes += 1
    assert num_isotopes == 4, num_isotopes

def test_validate_isotopes():
    isos = load_isotope_data('rproc')
    isos['Ba'][136] = 0.1
    try:
        validate_isotopes(isos)
    except IsotopeError:
        pass
    else:
        raise Error

def test_not_needed_isotopes():
    ll = LineList.read(datadir+'/linelists/lin4554new')
    isos = load_isotope_data('asplund') #Does not have Ba!
    found_isos, missing_isos = get_needed_isotopes(ll,isos)
    if len(missing_isos) != 1:
        raise Error("{}   {}".format(found_isos,missing_isos))

def test_moog_isotopes():
    output = moog.utils._format_isotopes(isotopes={'H-C':{112:.8,113:.2},
                                                   'C-N':{1214:.8,1314:.2},
                                                   'Ba':{134:.0,135:.370,136:.0,137:.351,138:.279}
                                                   },
                                         num_synth=3)
    #print(output)
    assert len(output.split('\n')) == 1+int(output.split('\n')[0].split()[0]),output
    for line in output.split('\n')[1:]:
        species,ifrac,a,b = line.split()
        

if __name__=="__main__":
    test_get_needed_isotopes()
    test_get_needed_isotopes2()
    test_identify_isotopes()
    test_identify_isotopes2()
    test_load_isotope_data()
    test_validate_isotopes()
    test_not_needed_isotopes()
    test_moog_isotopes()
