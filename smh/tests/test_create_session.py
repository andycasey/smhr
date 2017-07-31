from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os, time
from nose.tools import assert_equals, assert_almost_equals, ok_

from smh import Session, LineList
import smh.spectral_models as sm

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'

def test_create_session_and_analyze():
    file_to_write_to = datadir+"/.test_create_session.smh" 

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
    
    start = time.time()
    session.save(file_to_write_to, overwrite=True)
    print("Time to save session: {:.1f} (file is {}M)".format(time.time()-start, os.path.getsize(file_to_write_to)/1.e6))

    start = time.time()
    session2 = Session.load(file_to_write_to)
    print("Time to load session: {:.1f}".format(time.time()-start))
    
if __name__=="__main__":
    test_create_session_and_analyze()
