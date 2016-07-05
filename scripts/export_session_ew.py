from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import smh
from smh import Session
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)
from smh.linelists import LineList

import sys, os

if __name__=="__main__":
    assert len(sys.argv)==2,"Must pass in path to session object"
    path = sys.argv[1]
    if not os.path.exists(path):
        print("{} does not exist".format(path))
        sys.exit()
    print("Extracting equivalent widths from "+path)
    session = Session.load(path)
    
    metadata = session.metadata
    line_list = metadata["line_list"]
    spectral_models = metadata["spectral_models"]

    wavelengths = []
    expots = []
    speciess = []
    loggfs = []
    equivalent_widths = []
    
    for spectral_model in spectral_models:
        # Skip unacceptable spectral models
        if not spectral_model.is_acceptable: continue
        
        # Skip all synthesis measurements
        if not isinstance(spectral_model, ProfileFittingModel): continue

        # Skip spectral models without fitted results
        if "fitted_result" not in spectral_model.metadata: continue
        
        transition = spectral_model.transitions[0]
        wavelengths.append(transition["wavelength"])
        expots.append(transition["expot"])
        speciess.append(transition["species"])
        loggfs.append(transition["loggf"])
        fitted_result = spectral_model.metadata["fitted_result"]
        equivalent_widths.append(fitted_result[2]["equivalent_width"][0]*1000.)
    
    N = len(wavelengths)
    dummy = np.ones(N)*np.nan
    dummy2= ['' for x in range(N)]
    ll = LineList([wavelengths,expots,speciess,loggfs,dummy,dummy,equivalent_widths,dummy2],
             names=['wavelength','expot','species','loggf','damp_vdw','dissoc_E','equivalent_width','comments'])
    ll.sort('wavelength')
    ll.write("out.moog",format="moog")
