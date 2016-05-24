from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
from smh.session import Session
import yaml

datadir = os.path.dirname(os.path.abspath(__file__))+'/test_data'

def test_rv():
    session = Session([datadir+'/spectra/hd122563.fits'])
    with open(Session._default_settings_path, "rb") as fp:
        defaults = yaml.load(fp)
    session.metadata.update(defaults)
    session.metadata["rv"]["template_spectrum"] = datadir+'/spectra/hd122563.fits'
    rv, rv_uncertainty = session.rv_measure()
    session.rv_correct(rv)

if __name__=="__main__":
    test_rv()

