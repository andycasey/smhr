from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import smh
from smh import Session
from smh.spectral_models import (ProfileFittingModel, SpectralSynthesisModel)
from smh.linelists import LineList

if __name__=="__main__":
    path = "/Users/alexji/lib/smhr/smh/tests/test_data/test2.smh"
    print("Extracting equivalent widths from "+path)
    
    session = Session.load(path)
    
    metadata = session.metadata
    line_list = metadata["line_list"]
    spectral_models = metadata["spectral_models"]

    wavelengths = []
    speciess = []
    equivalent_widths = []
    equivalent_widths_err1 = []
    equivalent_widths_err2 = []
    abundances = []
    abundance_uncertainties = []

    for spectral_model in spectral_models:
        # Skip unacceptable spectral models
        if not spectral_model.is_acceptable: continue
        
        # Skip all synthesis measurements
        if not isinstance(model, ProfileFittingModel): continue

        # Skip spectral models without fitted results
        if "fitted_result" not in spectral_model.metadata: continue
        
        # Extract information
        transition = spectral_model.transitions[0]
        wavelengths.append(transition["wavelength"])
        speciess.append(transition["species"])

        # Equivalent width is saved in units of Angstroms
        # Turn it into milli-Angstroms by multiplying by 1000.
        equivalent_widths.append(fitted_result[2]["equivalent_width"][0]*1000.)
        equivalent_widths_err1.append(fitted_result[2]["equivalent_width"][1]*1000.)
        equivalent_widths_err2.append(fitted_result[2]["equivalent_width"][2]*1000.)
        
        # Extract abundances
        # You can have a fit without abundances
        # Use a try-except block to either add the abundance or NaN
        try:
            abundances.append(fitted_result[2]["abundances"][0])
        except KeyError:
            abundances.append(np.nan)
        
        try:
            abundance_uncertainties.append(fitted_result[2]["abundance_uncertainties"][0])
        except KeyError:
            abundance_uncertainties.append(np.nan)

    wavelengths = np.array(wavelengths)
    speciess = np.array(speciess)
    equivalent_widths = np.array(equivalent_widths)
    abundances = np.array(abundances)
    abundance_uncertainties = np.array(abundance_uncertainties)
