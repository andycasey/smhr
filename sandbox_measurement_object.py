

import smh
from smh import linelists

session = smh.Session([
    "/Users/arc/codes/smh/hd44007red_multi.fits",
    "/Users/arc/codes/smh/hd44007blue_multi.fits",
])

spectrum = smh.specutils.Spectrum1D.read("hd140283.fits")


# Check against best GES node using MOOG.
transitions = smh.linelists.LineList.read("/Users/arc/research/ges/linelist/vilnius.ew")

# Generate photosphere.
import smh.photospheres
interpolator = smh.photospheres.interpolator(kind="marcs")
photosphere = interpolator(5737, 4.395, 0.03)
photosphere.meta["stellar_parameters"]["microturbulence"] = 0.87

raise a
import smh.r
foo = smh.rt.moog.abundance_cog(photosphere, transitions)