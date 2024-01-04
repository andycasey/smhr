import os
import numpy as np
from typing import Union, Sequence
from functools import cached_property


from specutils import Spectrum1D


# Move to utilities:
        
def expand_path(path):
    return os.path.abspath(os.path.expanduser(path))

def expand_input_paths(input_paths):
    if isinstance(input_paths, (str, bytes)):
        input_paths = [input_paths]
    return tuple(map(expand_path, input_paths))


def get_closest_spectrum_index(session, wavelength):
    """Return the spectrum in the session closest to the wavelength."""
    mean_wavelengths = [np.mean(wl) for wl, *_ in session.spectra]
    index = np.argmin(np.abs(np.array(mean_wavelengths) - wavelength))
    return index


class Session:
    
    def __init__(
        self,
        input_paths: Sequence[Union[str, bytes, os.PathLike]]
    ):
        self.input_paths = expand_input_paths(input_paths)        
        return None
        
    
    @cached_property
    def spectra(self):
        print(self, hash(self), "loading spectra")
        spectra = []
        for input_path in self.input_paths:            
            wavelength, flux, ivar, meta = Spectrum1D.read_fits_multispec(input_path, flux_ext=6, ivar_ext=3)
            for i in range(len(wavelength)):
                spectra.append([wavelength[i], flux[i], ivar[i], meta])
        return spectra
    
    
    def _get_closest_spectrum_index(self, wavelength):
        return get_closest_spectrum_index(self, wavelength)
        
        

