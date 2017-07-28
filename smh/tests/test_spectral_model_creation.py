from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
from nose.tools import assert_equals, assert_almost_equals, ok_

from smh import Session, LineList
import smh.spectral_models as sm

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'
ll_filenames = ['masseron_linch.txt','lin4077new','lin4554new']
synth_lls = [LineList.read(datadir+'/linelists/'+filename) for filename in ll_filenames]
synth_elems = ["C", "Sr", "Ba"]
eqw_ll = LineList.read(datadir+'/linelists/complete.list')

def test_create_profile():
    session = Session([datadir+'/spectra/hd122563.fits'])
    model = sm.ProfileFittingModel(session, eqw_ll[0])
    model = sm.ProfileFittingModel(session, eqw_ll[0])
    print("Loaded profile model!")

def test_create_synthesis():
    session = Session([datadir+'/spectra/hd122563.fits'])
    for ll,elem in zip(synth_lls,synth_elems):
        model = sm.SpectralSynthesisModel(session, ll, elem)
    print("Loaded synthesis models!")

if __name__=="__main__":
    test_create_profile()
    test_create_synthesis()
