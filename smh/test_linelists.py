from linelists import *
from nose.tools import assert_equals, assert_almost_equals, ok_

ll_filenames = ['masseron_linch.txt']
ll_filenames = ['masseron_linch.txt','lin4077new','lin4554new','AJ_4398_Y.lin']

lls = [LineList.read_moog(filename) for filename in ll_filenames]
moog_lls = [LineList.read_moog(filename,full_columns=False) for filename in ll_filenames]

def test_colnames():
    required_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E',
                         'E_hi','E_lo','lande_hi','lande_lo','damp_stark','damp_rad','references']
    for ll in lls:
        for col in required_colnames:
            ok_(col in ll.data.columns)
def test_moog_colnames():
    required_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','references']
    for ll in moog_lls:
        for col in required_colnames:
            ok_(col in ll.data.columns)
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
    ll1.data = table.vstack([ll1.data,ll2.data[0:100]])
    duplicate_indices,duplicate_lines = ll1.find_duplicates()
    assert_equals(2*N,len(duplicate_indices))
    
test_colnames()
test_moog_colnames()
test_merge()
test_duplicates()
#ll = LineList.read('masseron_linch.txt')
#ll0 = LineList.read('lin4077new')
#ll.merge(ll0)
#ll.verbose=True
#ll.merge(ll0)
