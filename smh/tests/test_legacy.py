from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os, time
from nose.tools import assert_equals, assert_almost_equals, ok_

from smh import Session, LineList, legacy
import smh.spectral_models as sm

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'

def test_convert_v01():
    fname_in = datadir+"/hd122563_v0.1.smh"
    fname_out= datadir+"/.converted.hd122563_v0.1.smh"
    legacy.convert_v0_1_to_v0_2(fname_in, fname_out, overwrite=True)
    session = Session.load(fname_out)
    
    session.measure_abundances()
    
