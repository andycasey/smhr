from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os, time
from nose.tools import assert_equals, assert_almost_equals, ok_

from smh import Session, LineList
import smh.spectral_models as sm

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'

def test_create_session_and_analyze():
    file_to_write_to = datadir+"/.test_create_session.smh" 
    #file_to_write_to = datadir+"/test_create_session.smh" 

    session = Session([datadir+'/spectra/hd122563.fits'])
    session.import_linelist_as_profile_models(datadir+'/linelists/complete.list')
    #rv, e_rv = session.rv_measure()
    #print( rv, e_rv )
    session.rv_correct(0.0)

    # TODO put this in somewhere nice, as well as the rest of the normalization tab
    session.metadata["normalization"] = {
        "continuum": [1.0],
        "normalization_kwargs": [{}]
    }
    session.stitch_and_stack()
    
    start = time.time()
    num_fit = 0
    for model in session.spectral_models:
        try:
            res = model.fit()
        except (ValueError, RuntimeError, TypeError) as e:
            print("Fitting error", model)
            print(e)
        else:
            num_fit += 1
    print("Time to fit {}/{} models: {:.1f}".format(num_fit, len(session.spectral_models), time.time()-start))
    
    session.measure_abundances()

    session.save(file_to_write_to, overwrite=True)

def test_load_session():
    file_to_load = datadir+"/test_create_session.smh" 
    session = Session.load(file_to_load)

if __name__=="__main__":
    test_create_session_and_analyze()
    test_load_session()
    
