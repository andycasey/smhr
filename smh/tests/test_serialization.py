from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os, time
from nose.tools import assert_equals, assert_almost_equals, ok_

from smh import Session, LineList
import smh.spectral_models as sm
from six.moves import cPickle as pickle

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
    
    fmt = "{} {}: {:.3f}/{} = {:.3f}"

    # Time loading/writing
    Niter = 20
    ll = LineList.read(datadir+"/linelists/complete.list")
    ll = LineList.read(datadir+"/linelists/complete.list")[0:1]
    print ("N lines = {}".format(len(ll)))

    # Time pkl save
    start = time.time()
    fname = "./complete.pkl"
    for i in range(Niter):
        with open(fname, "w") as fp:
            pickle.dump(ll, fp)
    end = time.time()-start
    print( fmt.format("pkl", "save", end, Niter, end/Niter) )
    # Time pkl load
    start = time.time()
    fname = "./complete.pkl"
    for i in range(Niter):
        with open(fname, "r") as fp:
            _ll = pickle.load(fp)
    end = time.time()-start
    print( fmt.format("pkl", "load", end, Niter, end/Niter) )
    # Time fits save
    start = time.time()
    fname = "./complete.fits"
    for i in range(Niter):
        ll.write(fname, format='fits', overwrite=True)
    end = time.time()-start
    print( fmt.format("fits", "save", end, Niter, end/Niter) )
    # Time fits save
    start = time.time()
    fname = "./complete.fits"
    for i in range(Niter):
        _ll = LineList.read(fname, format='fits')
    end = time.time()-start
    print( fmt.format("fits", "load", end, Niter, end/Niter) )

    print ( "File size: pkl={}, fits={}".format(os.path.getsize("complete.pkl"), os.path.getsize("complete.fits")))
