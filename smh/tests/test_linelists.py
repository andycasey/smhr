from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os

from smh import linelists, utils
from smh.linelists import LineList
from nose.tools import assert_equals, assert_almost_equals, ok_
from astropy import table
import numpy as np

ll_filenames = ['masseron_linch.txt']
ll_filenames = ['masseron_linch.txt','lin4077new','lin4554new']

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'
lls = [LineList.read(datadir+'/linelists/'+filename) for filename in ll_filenames]
moog_lls = [LineList.read(datadir+'/linelists/'+filename,moog_columns=True) for filename in ll_filenames]

def test_species_converting():
    # Species, elem1, isotope1, elem2, isotope2, ion
    tests = [(56.1, 'Ba', 0, '', 0, 2),
             (56.0, 'Ba', 0, '', 0, 1),
             (56.1152, 'Ba', 152, '', 0, 2),
             (106.0,     'H', 0, 'C', 0, 1),
             (106.00112, 'H', 1, 'C', 12, 1),
             (106.00113, 'H', 1, 'C', 13, 1)
            ]
    failed = False
    for test in tests:
        species, elem1, isotope1, elem2, isotope2, ion = test
        input = (elem1, elem2, isotope1, isotope2, ion)
        output = utils.species_to_elems_isotopes_ion(species)
        for x,y in zip(input,output):
            if x != y: 
                print(test,x,y)
                failed = True
    assert not failed
    print("Passed species converting!")


def test_colnames():
    for ll in lls:
        required_colnames = ll.full_colnames
        for col in required_colnames:
            ok_(col in ll.columns)
def test_moog_colnames():
    for ll in moog_lls:
        required_colnames = ll.moog_colnames
        for col in required_colnames:
            ok_(col in ll.columns)
def test_merge():
    ll1 = lls[0]
    ll2 = lls[1]
    ll1.merge(ll2,raise_exception=False)
    lls[0] = LineList.read_moog(datadir+'/linelists/masseron_linch.txt')
def test_duplicates():
    fname = datadir+'/linelists/lin4077new'
    N = 50
    ll1 = LineList.read_moog(fname)
    ll2 = LineList.read_moog(fname)
    ll = table.vstack([ll1,ll2[0:N]])
    duplicate_indices,duplicate_lines = ll.find_duplicates(thresh=.01)
    assert_equals(2*N,len(duplicate_indices))

def test_readwrite_moog(fname=datadir+'/linelists/masseron_linch.txt'):
    ll = LineList.read(fname)
    ll = LineList.read(fname,format='moog')
    ll.write('_test.moog',format='moog')
    ll2 = LineList.read('_test.moog')
    for i in range(len(ll)):
        for col in ll[i].colnames:
            if isinstance(ll[i][col], str): continue
            diff = np.sum(np.abs(ll[i][col] - ll2[i][col]))
            assert np.isnan(diff) or diff < .001, "{} {}".format(ll[i][col],ll2[i][col])
    os.remove('_test.moog')
def test_writeread(fname=datadir+'/linelists/masseron_linch.txt'):
    ll = LineList.read(fname)
    ll.write('line_list.fits',format='fits')
    ll2 = LineList.read('line_list.fits',format='fits')
    for i in range(len(ll)):
        for col in ll[i].colnames:
            if isinstance(ll[i][col], str): continue
            diff = np.sum(np.abs(ll[i][col] - ll2[i][col]))
            assert np.isnan(diff) or diff < .001, "{} {}".format(ll[i][col],ll2[i][col])
    os.remove('line_list.fits')

def test_exception():
    ll = LineList.read(datadir+'/linelists/complete.list')
    N = len(ll)
    ll2 = LineList.read(datadir+'/linelists/complete.list')
    ll.merge(ll2,skip_equal_loggf=True)
    assert len(ll)==N,"{} {}".format(len(ll),N)
    
    ll2[5]['loggf'] = -3.14159
    ll2[200]['loggf'] = -3.14159
    try:
        ll.merge(ll2,skip_equal_loggf=True)
    except linelists.LineListConflict as e:
        assert len(e.conflicts1) == 2
        assert len(e.conflicts2) == 2
    else:
        raise RuntimeError("Conflicts didn't work")

if __name__=="__main__":
    test_species_converting()
    test_colnames()
    test_moog_colnames()
    test_merge()
    test_duplicates()
    for fname in ll_filenames:
        test_readwrite_moog(datadir+'/linelists/'+fname)
        test_writeread(datadir+'/linelists/'+fname)
    test_exception()
    
    ll = LineList.read(datadir+'/linelists/complete.list')
    ll.verbose = True
    ll_ti = LineList.read(datadir+'/linelists/tiII.moog')
    
    print("Merging to new object")
    new_ll = ll.merge(ll_ti, raise_exception=False, in_place=False)
    
    print("Merging in place")
    ll.merge(ll_ti, raise_exception=False)
    
    print("Sorting")
    ll.sort('wavelength')

#print("Reading GES list")
#ges = LineList.read('ges_master_v5.fits')
