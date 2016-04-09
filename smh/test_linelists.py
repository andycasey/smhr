from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import linelists
import utils
from linelists import LineList
from nose.tools import assert_equals, assert_almost_equals, ok_
from astropy import table

ll_filenames = ['masseron_linch.txt']
ll_filenames = ['masseron_linch.txt','lin4077new','lin4554new','AJ_4398_Y.lin']

lls = [LineList.read(filename) for filename in ll_filenames]
moog_lls = [LineList.read(filename,moog_columns=True) for filename in ll_filenames]

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
    #required_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E',
    #                     'E_hi','lande_hi','lande_lo','damp_stark','damp_rad','references']
    for ll in lls:
        required_colnames = ll.full_colnames
        for col in required_colnames:
            ok_(col in ll.columns)
def test_moog_colnames():
    #required_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','references']
    for ll in moog_lls:
        required_colnames = ll.moog_colnames
        for col in required_colnames:
            ok_(col in ll.columns)
def test_merge():
    ll1 = lls[0]
    ll2 = lls[1]
    ll1.merge(ll2)
    lls[0] = LineList.read_moog('masseron_linch.txt')
def test_duplicates():
    fname = 'lin4077new'
    N = 100
    ll1 = LineList.read_moog(fname)
    ll2 = LineList.read_moog(fname)
    ll = table.vstack([ll1,ll2[0:N]])
    duplicate_indices,duplicate_lines = ll.find_duplicates()
    assert_equals(2*N,len(duplicate_indices))
    
test_species_converting()
test_colnames()
test_moog_colnames()
test_merge()
test_duplicates()

ll = LineList.read('complete.list')
ll.verbose = True
ll_ncap = LineList.read('linelist_he1523_ncap_2016_03_smh_ohne.dat')

print("Merging to new object")
new_ll = ll.merge(ll_ncap, in_place=False)

print("Merging in place")
ll.merge(ll_ncap)

#print("Reading GES list")
#ges = LineList.read('ges_master_v5.fits')
