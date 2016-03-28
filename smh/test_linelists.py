from linelists import *
from nose.tools import assert_equals, assert_almost_equals, ok_

ll_filenames = ['masseron_linch.txt']
#ll_filenames = ['masseron_linch.txt','lin4077new','lin4554new','AJ_4398_Y.lin']

lls = [read_moog_linelist(filename) for filename in ll_filenames]
moog_lls = [read_moog_linelist(filename,full_columns=False) for filename in ll_filenames]

def test_colnames():
    required_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E',
                         'E_hi','E_lo','lande_hi','lande_lo','damp_stark','damp_rad','references']
    for ll in lls:
        for col in required_colnames:
            ok_(col in ll.columns)
def test_moog_colnames():
    required_colnames = ['wavelength','species','expot','loggf','damp_vdw','dissoc_E','references']
    for ll in moog_lls:
        for col in required_colnames:
            ok_(col in ll.columns)
    
test_colnames()
test_moog_colnames()
