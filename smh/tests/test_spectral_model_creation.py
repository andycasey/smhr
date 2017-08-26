from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
from nose.tools import assert_equals, assert_almost_equals, ok_

from smh import Session, LineList
import smh.spectral_models as sm

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'
synth_ll_filenames = [datadir+'/linelists/'+fname for fname in ['masseron_linch.txt','lin4077new','lin4554new']]
synth_lls = [LineList.read(filename) for filename in synth_ll_filenames]
synth_elems = ["C", "Sr", "Ba"]
eqw_ll_fname = datadir+'/linelists/complete.list'
eqw_ll = LineList.read(eqw_ll_fname)

def test_create_profile():
    session = Session([datadir+'/spectra/hd122563.fits'])
    model = sm.ProfileFittingModel(session, eqw_ll[0])
    model = sm.ProfileFittingModel(session, eqw_ll[0])
    print("Created profile model!")

def test_create_synthesis():
    session = Session([datadir+'/spectra/hd122563.fits'])
    for ll,elem in zip(synth_lls,synth_elems):
        model = sm.SpectralSynthesisModel(session, ll, elem)
    print("Created synthesis models!")

def test_import_models_into_session():
    session = Session([datadir+'/spectra/hd122563.fits'])
    session.import_linelist_as_profile_models(eqw_ll_fname)
    print ("Loaded profiles into session")
    for fname, elem in zip(synth_ll_filenames, synth_elems):
        session.import_linelist_as_synthesis_model(fname, elem)
    print ("Loaded syntheses into session")

def test_import_eqw_into_session():
    session = Session([datadir+'/spectra/hd122563.fits'])
    session.import_linelist_as_profile_models(datadir+"/linelists/frebel13_HD122563.moog2",
                                              import_equivalent_widths=True)
    print ("Loaded eqws into session")

if __name__=="__main__":
    test_create_profile()
    test_create_synthesis()
    test_import_models_into_session()
    test_import_eqw_into_session()
