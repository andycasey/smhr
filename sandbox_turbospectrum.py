
from smh.radiative_transfer import turbospectrum as ts
from smh.photospheres.marcs import Interpolator
from smh.linelists import LineList




marcs = Interpolator()
solar_photosphere = marcs(5777, 4.4, 0)
linelist = LineList.read("ni-test.list")

# Just restrict to a few lines.
#linelist = linelist[500:505]

ts.synthesize(solar_photosphere, linelist, None, dispersion_min=5150, dispersion_max=5165)
