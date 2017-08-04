from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os

from smh import specutils

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'

def test_load_and_write():
    spec = specutils.Spectrum1D.read(datadir+"/spectra/hd122563.fits")
    spec.write(datadir+"/spectra/.test_load_and_write.txt")
    spec.write(datadir+"/spectra/.test_load_and_write.fits")
    
