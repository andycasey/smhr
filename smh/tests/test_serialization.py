from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os, time
from nose.tools import assert_equals, assert_almost_equals, ok_

from smh import Session, LineList
import smh.spectral_models as sm
import cPickle as pickle

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'

def test_serialize_linelist():
    ll = LineList.read(datadir+"/linelists/complete.list")
    file_to_write_to = datadir+"/.test.serialize.linelist.pkl" 
    with open(file_to_write_to, "w") as fp:
        pickle.dump(ll, fp)
    with open(file_to_write_to, "r") as fp:
        loaded_ll = pickle.load(fp)
    
def test_load_serialized_linelist():
    ## previously dumped complete.list to this file
    file_to_load = datadir+"/test_load_serialized_linelist.pkl" 
    with open(file_to_load, "r") as fp:
        ll = pickle.load(fp)
    
if __name__=="__main__":
    test_serialize_linelist()
    test_load_serialized_linelist()
    
    # Time loading/writing
    
