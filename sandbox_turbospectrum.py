
from smh.radiative_transfer import turbospectrum as ts
from smh.photospheres.marcs import Interpolator

marcs = Interpolator()

solar_photosphere = marcs(5777, 4.4, 0)
